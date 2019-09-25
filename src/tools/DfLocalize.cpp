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

#include "DfLocalize.h"
#include <iostream>

// #define OCCU_THRESHOLD (1.0e-7)
// #define UNOCCU_THRESHOLD (1.0e-4)
// #define DELTAG_THRESHOLD (1e-14)

DfLocalize::DfLocalize(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
    this->maxIteration_ = 100;
    if ((*pPdfParam)["lo/max_iteration"].getStr().empty() != true) {
        this->maxIteration_ = (*pPdfParam)["lo/max_iteration"].getInt();
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

DfLocalize::~DfLocalize() {}

void DfLocalize::makeGroup() {
    const std::size_t numOfAOs = this->m_nNumOfAOs;
    const std::size_t numOfAtoms = this->m_nNumOfAtoms;

    this->group_.clear();
    this->group_.resize(numOfAtoms);
    for (std::size_t orb = 0; orb < numOfAOs; ++orb) {
        const std::size_t atomID = this->orbInfo_.getAtomIndex(orb);
        this->group_[atomID].push_back(orb);
    }
}

void DfLocalize::localize(const std::string& inputCMatrixPath) {
    std::cout << TlUtils::format("number of Atoms = %ld", this->m_nNumOfAtoms)
              << std::endl;
    std::cout << TlUtils::format("number of AOs   = %ld", this->m_nNumOfAOs)
              << std::endl;
    std::cout << TlUtils::format("number of MOs   = %ld", this->m_nNumOfMOs)
              << std::endl;
    this->makeGroup();

    const std::string S_path = this->getSpqMatrixPath();
    if (this->S_.load(S_path) == false) {
        std::cerr << "could not load: " << S_path << std::endl;
        abort();
    }

    std::string CMatrixPath = inputCMatrixPath;
    if (CMatrixPath.empty() == true) {
        CMatrixPath =
            this->getCMatrixPath(DfObject::RUN_RKS, this->m_nIteration);
    }
    if (this->C_.load(CMatrixPath) == false) {
        std::cerr << "could not load: " << CMatrixPath << std::endl;
        abort();
    }

    const int maxIteration = this->maxIteration_;
    std::cout << TlUtils::format("max iteration = %d", maxIteration)
              << std::endl;
    std::cout << TlUtils::format("threshold = %10.5e", this->threshold_)
              << std::endl;
    for (int num_iteration = 0; num_iteration < maxIteration; ++num_iteration) {
        this->makeQATable();

        this->makeJobList();

        double sumDeltaG = 0.0;
        JobItem jobItem;
        (void)this->getJobItem(&jobItem, true);
        while (this->getJobItem(&jobItem) != 0) {
            const std::size_t orb_i = jobItem.orb_i;
            const std::size_t orb_j = jobItem.orb_j;

            double A_ij, B_ij;
            this->calcQA(orb_i, orb_j, &A_ij, &B_ij);
            double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
            double deltaG = A_ij + normAB;

            if (std::fabs(deltaG) > 1.0E-16) {
                sumDeltaG += deltaG;
                TlDenseGeneralMatrix_Lapack rot(2, 2);
                this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
                this->rotateCmatrix(orb_i, orb_j, rot);
            }
        }

        std::cout << TlUtils::format("%d th: sum of delta_g: %10.5e",
                                     num_iteration + 1, sumDeltaG)
                  << std::endl;
        DfObject::saveCloMatrix(RUN_RKS, num_iteration + 1, this->C_);

        if (sumDeltaG < this->threshold_) {
            std::cout << "number of iteration: " << num_iteration + 1
                      << std::endl;
            (*(this->pPdfParam_))["lo/num_of_iterations"] = num_iteration + 1;
            (*(this->pPdfParam_))["lo/satisfied"] = "yes";
            break;
        }
    }
}

void DfLocalize::makeQATable() {
    std::size_t startOrb = this->startOrb_;
    std::size_t endOrb = this->endOrb_;

    // sort
    std::size_t size = endOrb - startOrb;
    this->orb_QA_table_.resize(size);
    for (std::size_t i = 0; i < size; ++i) {
        std::size_t orb = startOrb + i;
        double QAii = this->calcQA(orb);
        Orb_QA_Item item(orb, QAii);
        this->orb_QA_table_[i] = item;
    }
    // QAが大きい順にソート
    std::sort(this->orb_QA_table_.begin(), this->orb_QA_table_.end(),
              do_OrbQAItem_sort_functor_cmp());
}

void DfLocalize::rotateCmatrix(std::size_t orb_i, std::size_t orb_j,
                               const TlDenseGeneralMatrix_Lapack& rot) {
    const std::size_t numOfAOs = this->m_nNumOfAOs;
    TlDenseGeneralMatrix_Lapack A(2, numOfAOs);
    for (std::size_t i = 0; i < numOfAOs; ++i) {
        A.set(0, i, this->C_.get(i, orb_i));
        A.set(1, i, this->C_.get(i, orb_j));
    }

    const TlDenseGeneralMatrix_Lapack B = rot * A;  // rotate A to form B

    for (std::size_t i = 0; i < numOfAOs; ++i) {
        this->C_.set(i, orb_i, B.get(0, i));
        this->C_.set(i, orb_j, B.get(1, i));
    }
}

void DfLocalize::getRotatingMatrix(const double A_ij, const double B_ij,
                                   const double normAB,
                                   TlDenseGeneralMatrix_Lapack* pRot) {
    assert(pRot != NULL);
    pRot->resize(2, 2);

    const double cos4a = -A_ij / normAB;
    const double cos_a2 =
        0.5 * (1.0 + std::sqrt(0.5 * (cos4a + 1.0)));  // (cos_a)^2
    double sin_a = std::sqrt(1.0 - cos_a2);
    const double cos_a = std::sqrt(cos_a2);
    if (B_ij < 0.0) {
        // it means (sin4a < 0.0)
        sin_a = -sin_a;
    }

    pRot->set(0, 0, cos_a);
    pRot->set(0, 1, sin_a);
    pRot->set(1, 0, -sin_a);
    pRot->set(1, 1, cos_a);
}

void DfLocalize::calcQA(const std::size_t orb_i, const std::size_t orb_j,
                        double* pA_ij, double* pB_ij) {
    assert(pA_ij != NULL);
    assert(pB_ij != NULL);
    assert(this->C_.getNumOfRows() == this->m_nNumOfAOs);
    assert(this->C_.getNumOfCols() == this->m_nNumOfMOs);
    assert(this->S_.getNumOfRows() == this->m_nNumOfAOs);
    assert(this->S_.getNumOfCols() == this->m_nNumOfAOs);

    const std::size_t numOfGrps = this->group_.size();
    double sumAij = 0.0;
    double sumBij = 0.0;
    for (std::size_t grp = 0; grp < numOfGrps; ++grp) {
        const std::vector<std::size_t> orbList = this->group_[grp];
        const std::size_t dim = orbList.size();
        if (dim == 0) {
            continue;
        }

        TlDenseGeneralMatrix_Lapack C_i(dim, 1);
        TlDenseGeneralMatrix_Lapack C_j(dim, 1);
        for (std::size_t i = 0; i < dim; ++i) {
            const std::size_t orb = orbList[i];
            C_i.set(i, 0, this->C_.get(orb, orb_i));
            C_j.set(i, 0, this->C_.get(orb, orb_j));
        }

        TlDenseGeneralMatrix_Lapack Ci_ = C_i;
        Ci_.transposeInPlace();
        TlDenseGeneralMatrix_Lapack Cj_ = C_j;
        Cj_.transposeInPlace();

        TlDenseSymmetricMatrix_Lapack partS(dim);
        for (std::size_t i = 0; i < dim; ++i) {
            const std::size_t orb1 = orbList[i];
            partS.set(i, i, this->S_.get(orb1, orb1));
            for (std::size_t j = 0; j < i; ++j) {
                const std::size_t orb2 = orbList[j];
                partS.set(i, j, this->S_.get(orb1, orb2));
            }
        }

        const TlDenseGeneralMatrix_Lapack QAiiMat = Ci_ * partS * C_i;
        assert(QAiiMat.getNumOfRows() == 1);
        assert(QAiiMat.getNumOfCols() == 1);
        const double QAii = QAiiMat.get(0, 0);

        const TlDenseGeneralMatrix_Lapack QAjjMat = Cj_ * partS * C_j;
        assert(QAjjMat.getNumOfRows() == 1);
        assert(QAjjMat.getNumOfCols() == 1);
        const double QAjj = QAjjMat.get(0, 0);

        const TlDenseGeneralMatrix_Lapack QAijMat = Ci_ * partS * C_j;
        const TlDenseGeneralMatrix_Lapack QAjiMat = Cj_ * partS * C_i;
        assert(QAijMat.getNumOfRows() == 1);
        assert(QAijMat.getNumOfCols() == 1);
        assert(QAjiMat.getNumOfRows() == 1);
        assert(QAjiMat.getNumOfCols() == 1);
        const double QAij = (QAijMat.get(0, 0) + QAjiMat.get(0, 0)) * 0.5;

        const double QAii_QAjj = QAii - QAjj;
        sumAij += QAij * QAij - 0.25 * QAii_QAjj * QAii_QAjj;
        sumBij += QAij * QAii_QAjj;
    }

    *pA_ij = sumAij;
    *pB_ij = sumBij;
}

double DfLocalize::calcQA(const std::size_t orb_i) {
    const std::size_t numOfGrps = this->group_.size();
    double sumQAii = 0.0;
    for (std::size_t grp = 0; grp < numOfGrps; ++grp) {
        const std::vector<std::size_t> orbList = this->group_[grp];
        const std::size_t dim = orbList.size();
        if (dim == 0) {
            continue;
        }

        TlDenseGeneralMatrix_Lapack C_i(dim, 1);
        for (std::size_t i = 0; i < dim; ++i) {
            const std::size_t orb = orbList[i];
            C_i.set(i, 0, this->C_.get(orb, orb_i));
        }

        TlDenseGeneralMatrix_Lapack Ci_ = C_i;
        Ci_.transposeInPlace();

        TlDenseSymmetricMatrix_Lapack partS(dim);
        for (std::size_t i = 0; i < dim; ++i) {
            const std::size_t orb1 = orbList[i];
            partS.set(i, i, this->S_.get(orb1, orb1));
            for (std::size_t j = 0; j < i; ++j) {
                const std::size_t orb2 = orbList[j];
                partS.set(i, j, this->S_.get(orb1, orb2));
            }
        }

        const TlDenseGeneralMatrix_Lapack QAiiMat = Ci_ * partS * C_i;
        assert(QAiiMat.getNumOfRows() == 1);
        assert(QAiiMat.getNumOfCols() == 1);
        const double QAii = QAiiMat.get(0, 0);

        sumQAii += QAii;
    }

    return sumQAii;
}

void DfLocalize::makeJobList() {
    const std::size_t orbQATableSize = this->orb_QA_table_.size();
    const std::size_t maxJobIndex = orbQATableSize * (orbQATableSize - 1) /
                                    2;  // '-1'は対角項が必要無いため
    this->jobList_.resize(maxJobIndex);

    std::size_t index = 0;
    JobItem item;
    std::vector<Orb_QA_Item>::const_iterator itrEnd = this->orb_QA_table_.end();
    for (std::vector<Orb_QA_Item>::const_iterator p =
             this->orb_QA_table_.begin();
         p != itrEnd; ++p) {
        item.orb_i = p->orb;
        for (std::vector<Orb_QA_Item>::const_iterator q = p + 1; q != itrEnd;
             ++q) {
            item.orb_j = q->orb;
            this->jobList_[index] = item;
            ++index;
        }
    }
}

int DfLocalize::getJobItem(DfLocalize::JobItem* pJob, bool isInitialized) {
    assert(pJob != NULL);

    static std::size_t jobIndex = 0;
    static std::size_t maxJobIndex = 0;

    if (isInitialized == true) {
        jobIndex = 0;
        maxJobIndex = this->jobList_.size();

        return true;
    }

    int answer = 0;
    if (jobIndex < maxJobIndex) {
        *pJob = this->jobList_[jobIndex];
        ++jobIndex;
        answer = 1;
    }

    return answer;
}
