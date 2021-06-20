// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include <iostream>

#include "DfLocalize.h"

// #define OCCU_THRESHOLD (1.0e-7)
// #define UNOCCU_THRESHOLD (1.0e-4)
// #define DELTAG_THRESHOLD (1e-14)

DfLocalize::DfLocalize(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      log_(TlLogging::getInstance()),
      isRestart_(false),
      CMatrixPath_(""),
      orbInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
    // setup log
    this->log_.setFilePath("lo-output.txt");

    this->maxIteration_ = 100;
    if ((*pPdfParam)["lo/max_iteration"].getStr().empty() != true) {
        this->maxIteration_ = (*pPdfParam)["lo/max_iteration"].getInt();
    }

    this->lo_iteration_ = 0;
    if ((*pPdfParam_)["lo/num_of_iterations"].getStr().empty() != true) {
        this->lo_iteration_ = (*pPdfParam_)["lo/num_of_iterations"].getInt();
    }

    this->threshold_ = 1.0E-4;
    if ((*pPdfParam)["lo/threshold"].getStr().empty() != true) {
        this->threshold_ = (*pPdfParam)["lo/threshold"].getDouble();
    }

    this->numOfOcc_ = (this->m_nNumOfElectrons + 1) / 2;

    // for OCC
    this->startOrb_ = 0;
    this->endOrb_ = this->numOfOcc_;
}

DfLocalize::~DfLocalize() {
}

void DfLocalize::setRestart(bool yn) {
    this->isRestart_ = yn;
}

void DfLocalize::setCMatrixPath(const std::string& path) {
    this->CMatrixPath_ = path;
}

void DfLocalize::initialize() {
    this->G_ = 0.0;

    // output information
    this->log_.info(TlUtils::format("number of Atoms = %d", this->m_nNumOfAtoms));
    this->log_.info(TlUtils::format("number of AOs   = %d", this->m_nNumOfAOs));
    this->log_.info(TlUtils::format("number of MOs   = %d", this->m_nNumOfMOs));
    this->log_.info(TlUtils::format("max iteration = %d", this->maxIteration_));
    this->log_.info(TlUtils::format("threshold = %10.5e", this->threshold_));

    this->log_.info("make group");
    this->makeGroup();

    // load S matrix
    this->log_.info("load S matrix");
    const bool isLoaded = this->getSMatrix(&(this->S_));
    if (!isLoaded) {
        const std::string path = this->getSpqMatrixPath();
        this->log_.critical(TlUtils::format("could not load: %s", path.c_str()));
        abort();
    }
}

bool DfLocalize::getSMatrix(TlDenseSymmetricMatrix_Lapack* pS) {
    const std::string path = this->getSpqMatrixPath();
    const bool isLoaded = pS->load(path);

    return isLoaded;
}

std::string DfLocalize::getCMatrixPath() {
    std::string path = "";
    if (this->CMatrixPath_.empty()) {
        if (this->isRestart_) {
            path = DfObject::getCloMatrixPath(DfObject::RUN_RKS, this->lo_iteration_);
        } else {
            const int scf_iteration = (*(this->pPdfParam_))["num_of_iterations"].getInt();
            path = DfObject::getCMatrixPath(DfObject::RUN_RKS, scf_iteration);
            this->lo_iteration_ = 0;
        }
    } else {
        path = this->CMatrixPath_;
    }

    return path;
}

void DfLocalize::getCMatrix(TlDenseGeneralMatrix_Lapack* pC) {
    this->log_.info("load C matrix");
    const std::string path = this->getCMatrixPath();

    this->log_.info(TlUtils::format("load C matrix: %s", path.c_str()));
    pC->load(path);
}

void DfLocalize::exec() {
    this->initialize();

    TlDenseGeneralMatrix_Lapack C;
    this->getCMatrix(&C);

    const int maxIteration = this->maxIteration_;
    for (int i = this->lo_iteration_ + 1; i <= maxIteration; ++i) {
        this->log_.info(TlUtils::format("localize: %d", i));
        const double sumDeltaG = this->localize(&C);

        const std::string msg = TlUtils::format("%d th: G=%10.5e, delta_G=%10.5e", i, this->G_, sumDeltaG);
        this->log_.info(msg);
        std::cout << msg << std::endl;

        DfObject::saveCloMatrix(RUN_RKS, i, C);

        if (sumDeltaG < this->threshold_) {
            (*(this->pPdfParam_))["lo/num_of_iterations"] = i;
            (*(this->pPdfParam_))["lo/satisfied"] = "yes";
            break;
        }
    }
}

double DfLocalize::localize(TlDenseGeneralMatrix_Lapack* pC) {
    const index_type numOfAOs = this->m_nNumOfAOs;

    this->makeQATable(*pC);
    this->makeJobList();

    this->log_.info("localize: start");
    double sumDeltaG = 0.0;
    const double rotatingThreshold = this->threshold_ * 0.01;

    JobItem jobItem;
    (void)this->getJobItem(&jobItem, true);
    while (this->getJobItem(&jobItem) != 0) {
        const index_type orb_i = jobItem.orb_i;
        const index_type orb_j = jobItem.orb_j;
        // this->log_.info(TlUtils::format("calc: %d, %d", orb_i, orb_j));

        double A_ij, B_ij;
        TlDenseVector_Lapack Cpi = pC->getColVector<TlDenseVector_Lapack>(orb_i);
        TlDenseVector_Lapack Cpj = pC->getColVector<TlDenseVector_Lapack>(orb_j);
        this->calcQA_ij(Cpi, Cpj, &A_ij, &B_ij);
        const double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
        const double deltaG = A_ij + normAB;
        // std::cout << TlUtils::format("%f: A=%f B=%f", deltaG, A_ij, B_ij) << std::endl;

        if (std::fabs(deltaG) > rotatingThreshold) {
            sumDeltaG += deltaG;
            TlDenseGeneralMatrix_Lapack rot(2, 2);
            this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
            this->rotateVectors(&Cpi, &Cpj, rot);

            pC->setColVector(orb_i, Cpi.data(), numOfAOs);
            pC->setColVector(orb_j, Cpj.data(), numOfAOs);
        }
    }

    this->log_.info("localize: end");
    return sumDeltaG;
}

void DfLocalize::makeGroup() {
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfAtoms = this->m_nNumOfAtoms;

    const int numOfGroups = numOfAtoms;

    this->groupV_.clear();
    this->groupV_.resize(numOfGroups);
    for (int i = 0; i < numOfGroups; ++i) {
        this->groupV_[i].resize(numOfAOs);
    }

    for (index_type orb = 0; orb < numOfAOs; ++orb) {
        const index_type atomID = this->orbInfo_.getAtomIndex(orb);
        this->groupV_[atomID].set(orb, 1.0);
    }
}

void DfLocalize::makeQATable(const TlDenseGeneralMatrix_Lapack& C) {
    this->log_.info("make QA table");

    const index_type startOrb = this->startOrb_;
    const index_type endOrb = this->endOrb_;

    double sumQAii2 = 0.0;
    // sort
    index_type size = endOrb - startOrb;
    this->orb_QA_table_.resize(size);
#pragma omp parallel for schedule(runtime) reduction(+ : sumQAii2)
    for (index_type i = 0; i < size; ++i) {
        const index_type orb = startOrb + i;
        const double QAii2 = this->calcQA_ii(C, orb);
        sumQAii2 += QAii2;
        Orb_QA_Item item(orb, QAii2);
        this->orb_QA_table_[i] = item;
    }
    this->G_ = sumQAii2;

    // QAが大きい順にソート
    this->log_.info("sort QA table");
    std::sort(this->orb_QA_table_.begin(), this->orb_QA_table_.end(), do_OrbQAItem_sort_functor_cmp());

    this->log_.info("end QA table");
}

void DfLocalize::rotateVectors(TlDenseVector_Lapack* pCpi, TlDenseVector_Lapack* pCpj,
                               const TlDenseGeneralMatrix_Lapack& rot) {
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(pCpi->getSize() == numOfAOs);
    assert(pCpj->getSize() == numOfAOs);

    TlDenseGeneralMatrix_Lapack A(2, numOfAOs);
    A.setRowVector(0, pCpi->data(), numOfAOs);
    A.setRowVector(1, pCpj->data(), numOfAOs);

    const TlDenseGeneralMatrix_Lapack B = rot * A;  // rotate A to form B

    *pCpi = B.getRowVector<TlDenseVector_Lapack>(0);
    *pCpj = B.getRowVector<TlDenseVector_Lapack>(1);
    assert(pCpi->getSize() == numOfAOs);
    assert(pCpj->getSize() == numOfAOs);
}

void DfLocalize::getRotatingMatrix(const double A_ij, const double B_ij, const double normAB,
                                   TlDenseGeneralMatrix_Lapack* pRot) {
    assert(pRot != NULL);
    pRot->resize(2, 2);

    const double cos4a = -A_ij / normAB;

    // cos(2\theta) = 2*(cos(\theta))^2 -1
    const double cos2a = std::sqrt(0.5 * (cos4a + 1.0));
    const double cos_a_2 = 0.5 * (cos2a + 1.0);  // (cos_a)^2
    const double cos_a = std::sqrt(cos_a_2);
    double sin_a = std::sqrt(1.0 - cos_a_2);

    // (gamma = alpha) means sin_g = sin_a, cos_g = cos_a
    if (B_ij < 0.0) {
        // it means (sin4a < 0.0)
        sin_a = -sin_a;
    }

    pRot->set(0, 0, cos_a);
    pRot->set(0, 1, sin_a);
    pRot->set(1, 0, -sin_a);
    pRot->set(1, 1, cos_a);
}

void DfLocalize::calcQA_ij(const TlDenseVector_Lapack& Cpi, const TlDenseVector_Lapack& Cpj, double* pA_ij,
                           double* pB_ij) {
    assert(pA_ij != NULL);
    assert(pB_ij != NULL);

    const index_type numOfGrps = this->groupV_.size();

    const TlDenseVector_Lapack SCpi = this->S_ * Cpi;
    const TlDenseVector_Lapack SCpj = this->S_ * Cpj;

    double sumAij = 0.0;
    double sumBij = 0.0;

#pragma omp parallel for schedule(runtime) reduction(+ : sumAij, sumBij)
    for (index_type grp = 0; grp < numOfGrps; ++grp) {
        TlDenseVector_Lapack tCqi = Cpi;
        TlDenseVector_Lapack tCqj = Cpj;
        tCqi.dotInPlace(this->groupV_[grp]);
        tCqj.dotInPlace(this->groupV_[grp]);

        const double QAii = tCqi * SCpi;
        const double QAjj = tCqj * SCpj;
        const double QAij = 0.5 * (tCqj * SCpi + tCqi * SCpj);

        const double QAii_QAjj = QAii - QAjj;
        sumAij += QAij * QAij - 0.25 * QAii_QAjj * QAii_QAjj;
        sumBij += QAij * QAii_QAjj;
    }

    *pA_ij = sumAij;
    *pB_ij = sumBij;
}

double DfLocalize::calcQA_ii(const TlDenseGeneralMatrix_Lapack& C, const index_type orb_i) {
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfGrps = this->groupV_.size();
    double sumQAii2 = 0.0;

    const TlDenseVector_Lapack Cpi = C.getColVector<TlDenseVector_Lapack>(orb_i);
    const TlDenseVector_Lapack SCpi = this->S_ * Cpi;

    // #pragma omp parallel for reduction(+ : sumQAii2)
    for (index_type grp = 0; grp < numOfGrps; ++grp) {
        TlDenseVector_Lapack tCqi = Cpi;
        tCqi.dotInPlace(this->groupV_[grp]);

        const double QAii = tCqi * SCpi;
        sumQAii2 += QAii * QAii;
    }

    return sumQAii2;
}

void DfLocalize::makeJobList() {
    const std::size_t orbQATableSize = this->orb_QA_table_.size();
    this->jobList_.clear();

    JobItem item;
    for (std::size_t i = 0; i < orbQATableSize; ++i) {
        item.orb_i = this->orb_QA_table_[i].orb;
        for (std::size_t j = orbQATableSize - 1; j > i; --j) {
            item.orb_j = this->orb_QA_table_[j].orb;
            this->jobList_.push_back(item);
        }
    }
}

bool DfLocalize::getJobItem(DfLocalize::JobItem* pJob, bool isInitialized) {
    assert(pJob != NULL);

    static std::list<JobItem>::iterator it;
    if (isInitialized == true) {
        it = this->jobList_.begin();

        return true;
    }

    bool answer = false;
    if (it != this->jobList_.end()) {
        *pJob = *it;
        ++it;
        answer = true;
    }

    return answer;
}
