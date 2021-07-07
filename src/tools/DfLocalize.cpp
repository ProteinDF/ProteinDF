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
#include "TlTime.h"

// #define OCCU_THRESHOLD (1.0e-7)
// #define UNOCCU_THRESHOLD (1.0e-4)
// #define DELTAG_THRESHOLD (1e-14)

DfLocalize::DfLocalize(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      log_(TlLogging::getInstance()),
      isRestart_(false),
      CMatrixPath_(""),
      orbInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
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
    // for OCC
    this->startOrb_ = 0;
    this->endOrb_ = this->numOfOcc_;

    this->G_ = 0.0;

    // output information
    this->log_.info(TlUtils::format("number of Atoms: %d", this->m_nNumOfAtoms));
    this->log_.info(TlUtils::format("number of AOs: %d", this->m_nNumOfAOs));
    this->log_.info(TlUtils::format("number of MOs: %d", this->m_nNumOfMOs));
    this->log_.info(TlUtils::format("max iteration: %d", this->maxIteration_));
    this->log_.info(TlUtils::format("threshold: %10.5e", this->threshold_));
    this->log_.info(TlUtils::format("processing MO: %d -> %d", this->startOrb_, this->endOrb_));

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
    assert(this->S_.getNumOfRows() == this->m_nNumOfAOs);

    this->checkOpenMP();
}

void DfLocalize::checkOpenMP() {
#ifdef _OPENMP
    {
        this->log_.info(">>>> OpenMP info ");
        const int numOfOmpThreads = omp_get_max_threads();
        this->log_.info(TlUtils::format("OpenMP max threads: %d", numOfOmpThreads));

        omp_sched_t kind;
        int modifier = 0;
        omp_get_schedule(&kind, &modifier);
        switch (kind) {
            case omp_sched_static:
                this->log_.info("OpenMP schedule: static");
                break;
            case omp_sched_dynamic:
                this->log_.info("OpenMP schedule: dynamic");
                break;
            case omp_sched_guided:
                this->log_.info("OpenMP schedule: guided");
                break;
            case omp_sched_auto:
                this->log_.info("OpenMP schedule: auto");
                break;
            default:
                this->log_.info(TlUtils::format("OpenMP schedule: unknown: %d", static_cast<int>(kind)));
                break;
        }

        this->log_.info("<<<<");
    }
#else
    { this->log_.info("OpenMP is not enabled."); }
#endif  // _OPENMP
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
    const std::string path = this->getCMatrixPath();

    this->log_.info(TlUtils::format("load C matrix: %s", path.c_str()));
    pC->load(path);
}

void DfLocalize::exec() {
    this->initialize();

    TlDenseGeneralMatrix_Lapack C;
    {
        TlDenseGeneralMatrix_Lapack C_full;
        this->getCMatrix(&C_full);
        assert(this->startOrb_ < C_full.getNumOfCols());
        assert(this->endOrb_ < C_full.getNumOfCols());

        C_full.block(0, this->startOrb_, C_full.getNumOfRows(), this->endOrb_ - this->startOrb_ + 1, &C);
    }

    const int maxIteration = this->maxIteration_;
    for (int i = this->lo_iteration_ + 1; i <= maxIteration; ++i) {
        this->log_.info(TlUtils::format("localize: %d", i));
        // const double sumDeltaG = this->localize(&C);
        const double sumDeltaG = this->localize_byPop(&C);

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
    this->G_ = this->calcG(*pC, this->startOrb_, this->endOrb_);
    const double sumDeltaG = this->localize_core(pC, this->startOrb_, this->endOrb_, this->startOrb_, this->endOrb_);

    return sumDeltaG;
}

double DfLocalize::localize_core(TlDenseGeneralMatrix_Lapack* pC, const index_type startMO1, const index_type endMO1,
                                 const index_type startMO2, const index_type endMO2) {
    assert(startMO1 < endMO1);
    assert(startMO2 < endMO2);

    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(pC->getNumOfRows() == numOfAOs);

    this->log_.info("localize: start");

    std::vector<TaskItem> taskList;
    if ((startMO1 == startMO2) && (endMO1 == endMO2)) {
        taskList = this->getTaskList(startMO1, endMO1);
    } else {
        taskList = this->getTaskList(startMO1, endMO1, startMO2, endMO2);
    }

    const double rotatingThreshold = this->threshold_ * 0.01;
    double sumDeltaG = 0.0;

    assert(endMO1 <= endMO2);
    this->initLockMO(endMO2 + 1);

    const std::size_t numOfTasks = taskList.size();
#pragma omp parallel for reduction(+ : sumDeltaG) schedule(runtime)
    for (std::size_t i = 0; i < numOfTasks; ++i) {
        const TaskItem& task = taskList[i];
        const index_type orb_i = task.first;
        const index_type orb_j = task.second;

        double sleep = 100.0;
        while (this->isLockedMO(orb_i, orb_j)) {
            TlTime::sleep(sleep);
            sleep = std::min(4000.0, sleep * 2.0);
        }
#pragma omp critical(lock_MO)
        {
            this->lockMO(orb_i);
            this->lockMO(orb_j);
        }

        // this->log_.info(TlUtils::format("calc: %d, %d", orb_i, orb_j));
        TlDenseVector_Lapack Cpi = pC->getColVector(orb_i);
        TlDenseVector_Lapack Cpj = pC->getColVector(orb_j);

        double A_ij = 0.0;
        double B_ij = 0.0;
        this->calcQA_ij(Cpi, Cpj, &A_ij, &B_ij);

        const double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
        const double deltaG = A_ij + normAB;

        if (std::fabs(deltaG) > rotatingThreshold) {
            sumDeltaG += deltaG;
            TlDenseGeneralMatrix_Lapack rot(2, 2);
            this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
            this->rotateVectors(&Cpi, &Cpj, rot);

            pC->setColVector(orb_i, numOfAOs, Cpi.data());
            pC->setColVector(orb_j, numOfAOs, Cpj.data());
        }

#pragma omp critical(unlock_MO)
        {
            this->unlockMO(orb_i);
            this->unlockMO(orb_j);
        }
    }

    this->log_.info("localize: end");
    return sumDeltaG;
}

void DfLocalize::makeGroup() {
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfAtoms = this->m_nNumOfAtoms;

    const int numOfGroups = numOfAtoms;

    this->group_.resize(numOfGroups);
    for (int i = 0; i < numOfGroups; ++i) {
        this->group_[i].resize(numOfAOs);
    }

    for (index_type orb = 0; orb < numOfAOs; ++orb) {
        const index_type atomID = this->orbInfo_.getAtomIndex(orb);
        this->group_[atomID].set(orb, 1.0);
    }
}

double DfLocalize::calcG(const TlDenseGeneralMatrix_Lapack& C, const index_type startMO, const index_type endMO) {
    this->log_.info("calc G");
    double G = 0.0;

    const index_type size = endMO - startMO;
#pragma omp parallel for schedule(runtime) reduction(+ : G)
    for (index_type i = 0; i < size; ++i) {
        const index_type orb = startMO + i;
        const double QAii2 = this->calcQA_ii(C, orb);
        G += QAii2;
    }

    return G;
}

double DfLocalize::calcG_sort(const TlDenseGeneralMatrix_Lapack& C, const index_type startMO, const index_type endMO) {
    this->log_.info("calc G");
    double G = 0.0;

    const index_type size = endMO - startMO;
    this->groupMoPops_.resize(size);
#pragma omp parallel for schedule(runtime) reduction(+ : G)
    for (index_type i = 0; i < size; ++i) {
        const index_type mo = startMO + i;
        const double QAii2 = this->calcQA_ii(C, mo);
        G += QAii2;

        this->groupMoPops_[i] = MoPop(mo, QAii2);
    }

    // QAが大きい順にソート
    this->log_.info("sort QA table");
    std::sort(this->groupMoPops_.begin(), this->groupMoPops_.end(), MoPop::MoPop_sort_functor_cmp());

    return G;
}

void DfLocalize::rotateVectors(TlDenseVector_Lapack* pCpi, TlDenseVector_Lapack* pCpj,
                               const TlDenseGeneralMatrix_Lapack& rot) {
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(pCpi->getSize() == numOfAOs);
    assert(pCpj->getSize() == numOfAOs);

    TlDenseGeneralMatrix_Lapack A(2, numOfAOs);
    A.setRowVector(0, numOfAOs, pCpi->data());
    A.setRowVector(1, numOfAOs, pCpj->data());

    const TlDenseGeneralMatrix_Lapack B = rot * A;  // rotate A to form B

    *pCpi = B.getRowVector(0);
    *pCpj = B.getRowVector(1);
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

    const index_type numOfGrps = this->group_.size();

    const TlDenseVector_Lapack SCpi = this->S_ * Cpi;
    const TlDenseVector_Lapack SCpj = this->S_ * Cpj;

    double sumAij = 0.0;
    double sumBij = 0.0;

    //#pragma omp for schedule(runtime) private(sumAij, sumBij) reduction(+ : sumAij, sumBij)
    // #pragma omp for schedule(runtime)
    for (index_type grp = 0; grp < numOfGrps; ++grp) {
        TlDenseVector_Lapack tCqi = Cpi;
        TlDenseVector_Lapack tCqj = Cpj;
        tCqi.dotInPlace(this->group_[grp]);
        tCqj.dotInPlace(this->group_[grp]);

        const double QAii = tCqi * SCpi;
        const double QAjj = tCqj * SCpj;
        const double QAij = 0.5 * (tCqj * SCpi + tCqi * SCpj);

        const double QAii_QAjj = QAii - QAjj;
        sumAij += QAij * QAij - 0.25 * QAii_QAjj * QAii_QAjj;
        sumBij += QAij * QAii_QAjj;
    }

    // #pragma omp critical(calcQA_ij)
    {
        *pA_ij += sumAij;
        *pB_ij += sumBij;
    }
}

double DfLocalize::calcQA_ii(const TlDenseGeneralMatrix_Lapack& C, const index_type orb_i) {
    // const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfGrps = this->group_.size();
    double sumQAii2 = 0.0;

    const TlDenseVector_Lapack Cpi = C.getColVector(orb_i);
    const TlDenseVector_Lapack SCpi = this->S_ * Cpi;

    // #pragma omp parallel for reduction(+ : sumQAii2)
    for (index_type grp = 0; grp < numOfGrps; ++grp) {
        TlDenseVector_Lapack tCqi = Cpi;
        tCqi.dotInPlace(this->group_[grp]);

        const double QAii = tCqi * SCpi;
        sumQAii2 += QAii * QAii;
    }

    return sumQAii2;
}

std::vector<DfLocalize::TaskItem> DfLocalize::getTaskList(const index_type startMO, const index_type endMO) {
    const index_type dim = endMO - startMO + 1;
    const index_type pairs = dim * (dim - 1) / 2;
    std::vector<TaskItem> taskList(pairs);

    std::size_t taskIndex = 0;
    TaskItem item;
    for (int index1 = 1; index1 < dim; ++index1) {
        // Starting from 1 instead of 0 means that the diagonal elements are not computed.
        const int max_index2 = dim - index1;
        for (int index2 = 0; index2 < max_index2; ++index2) {
            const index_type orb_i = startMO + index2;
            const index_type orb_j = startMO + index1 + index2;

            taskList[taskIndex] = std::make_pair(orb_i, orb_j);
            ++taskIndex;
        }
    }
    assert(taskIndex == pairs);

    return taskList;
}

// for derived class
std::vector<DfLocalize::TaskItem> DfLocalize::getTaskList(const index_type startMO1, const index_type endMO1,
                                                          const index_type startMO2, const index_type endMO2) {
    assert(startMO1 < endMO1);
    assert(endMO1 < startMO2);
    assert(startMO2 < endMO2);

    const index_type pairs = (endMO1 - startMO1 + 1) * (endMO2 - startMO2 + 1);
    std::vector<TaskItem> taskList(pairs);

    std::size_t taskIndex = 0;
    for (index_type mo1 = startMO1; mo1 <= endMO1; ++mo1) {
        for (index_type mo2 = startMO2; mo2 <= endMO2; ++mo2) {
            taskList[taskIndex] = std::make_pair(mo1, mo2);
            ++taskIndex;
        }
    }
    assert(taskIndex == pairs);

    return taskList;
}

// bool DfLocalize::getTaskItem(DfLocalize::TaskItem* pTask, bool isInitialized) {
//     assert(pTask != NULL);

//     static std::vector<TaskItem>::iterator it;
//     static std::size_t counter = 0;
//     static std::size_t maxTasks = 0;
//     static std::size_t noticeNext = 0;

//     if (isInitialized == true) {
//         it = this->taskList_.begin();

//         counter = 0;
//         maxTasks = this->taskList_.size();
//         noticeNext = 1;

//         return true;
//     }

//     bool answer = false;
//     if (it != this->taskList_.end()) {
//         *pTask = *it;
//         ++it;

//         if (counter > (maxTasks / 10) * noticeNext) {
//             const double ratio = (double(counter) / double(maxTasks)) * 100.0;
//             this->log_.info(TlUtils::format("progress: %3.2f%% (%ld/%ld)", ratio, counter, maxTasks));
//             ++noticeNext;
//         }

//         ++counter;
//         answer = true;
//     } else {
//         this->log_.info("progress: 100% done");
//     }

//     return answer;
// }

// lock MO
void DfLocalize::initLockMO(const index_type numOfMOs) {
    this->lockMOs_.resize(numOfMOs);
    std::fill(this->lockMOs_.begin(), this->lockMOs_.end(), 0);
}

void DfLocalize::lockMO(const index_type mo) {
    assert(mo < static_cast<index_type>(this->lockMOs_.size()));
    this->lockMOs_[mo] = 1;
}

void DfLocalize::unlockMO(const index_type mo) {
    assert(mo < static_cast<index_type>(this->lockMOs_.size()));
    this->lockMOs_[mo] = 0;
}

bool DfLocalize::isLockedMO(const index_type mo1, const index_type mo2) const {
    assert(mo1 < static_cast<index_type>(this->lockMOs_.size()));
    assert(mo2 < static_cast<index_type>(this->lockMOs_.size()));
    bool answer = true;

    if ((this->lockMOs_[mo1] + this->lockMOs_[mo2] == 0)) {
        answer = false;
    }

    return answer;
}

// sort
double DfLocalize::localize_byPop(TlDenseGeneralMatrix_Lapack* pC) {
    this->G_ = this->calcG_sort(*pC, this->startOrb_, this->endOrb_);
    const double sumDeltaG = this->localize_core_byPop(pC, this->startOrb_, this->endOrb_);

    return sumDeltaG;
}

double DfLocalize::localize_core_byPop(TlDenseGeneralMatrix_Lapack* pC, const index_type startMO1,
                                       const index_type endMO1) {
    assert(startMO1 < endMO1);

    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(pC->getNumOfRows() == numOfAOs);

    this->log_.info("localize: start");

    const std::vector<TaskItem> taskList = this->getTaskList_byPop();

    const double rotatingThreshold = this->threshold_ * 0.01;
    double sumDeltaG = 0.0;

    this->initLockMO(endMO1 + 1);

    const std::size_t numOfTasks = taskList.size();
#pragma omp parallel for reduction(+ : sumDeltaG) schedule(runtime)
    for (std::size_t i = 0; i < numOfTasks; ++i) {
        const TaskItem& task = taskList[i];
        const index_type orb_i = task.first;
        const index_type orb_j = task.second;

        double sleep = 100.0;
        while (this->isLockedMO(orb_i, orb_j)) {
            TlTime::sleep(sleep);
            sleep = std::min(4000.0, sleep * 2.0);
        }
#pragma omp critical(lock_MO)
        {
            this->lockMO(orb_i);
            this->lockMO(orb_j);
        }

        // this->log_.info(TlUtils::format("calc: %d, %d", orb_i, orb_j));
        TlDenseVector_Lapack Cpi = pC->getColVector(orb_i);
        TlDenseVector_Lapack Cpj = pC->getColVector(orb_j);

        double A_ij = 0.0;
        double B_ij = 0.0;
        this->calcQA_ij(Cpi, Cpj, &A_ij, &B_ij);

        const double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
        const double deltaG = A_ij + normAB;

        if (std::fabs(deltaG) > rotatingThreshold) {
            sumDeltaG += deltaG;
            TlDenseGeneralMatrix_Lapack rot(2, 2);
            this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
            this->rotateVectors(&Cpi, &Cpj, rot);

            pC->setColVector(orb_i, numOfAOs, Cpi.data());
            pC->setColVector(orb_j, numOfAOs, Cpj.data());
        }

#pragma omp critical(unlock_MO)
        {
            this->unlockMO(orb_i);
            this->unlockMO(orb_j);
        }
    }

    this->log_.info("localize: end");
    return sumDeltaG;
}

std::vector<DfLocalize::TaskItem> DfLocalize::getTaskList_byPop() {
    const index_type dim = this->groupMoPops_.size();
    const index_type pairs = dim * (dim - 1) / 2;
    std::vector<TaskItem> taskList(pairs);

    std::size_t taskIndex = 0;
    TaskItem item;
    // for (int index1 = 0; index1 < dim; ++index1) {
    //     const index_type mo1 = this->groupMoPops_[index1].mo;
    //     for (int index2 = dim - 1; index2 > index1; --index2) {
    //         const index_type mo2 = this->groupMoPops_[index2].mo;

    //         taskList[taskIndex] = std::make_pair(mo1, mo2);
    //         ++taskIndex;
    //     }
    // }

    // big-big pair
    // for (int index1 = 0; index1 < dim - 1; ++index1) {
    //     const int max_index2 = dim - index1 - 1;
    //     for (int index2 = 0; index2 < max_index2; ++index2) {
    //         const index_type mo1 = this->groupMoPops_[index1 + index2 + 1].mo;
    //         const index_type mo2 = this->groupMoPops_[index2].mo;
    //         taskList[taskIndex] = std::make_pair(mo1, mo2);
    //         ++taskIndex;
    //     }
    // }

    // small-small pair
    for (int index1 = 0; index1 < dim - 1; ++index1) {
        const int max_index2 = dim - index1 - 1;
        for (int index2 = 0; index2 < max_index2; ++index2) {
            const index_type mo1 = this->groupMoPops_[dim - (index1 + index2) - 1].mo;
            const index_type mo2 = this->groupMoPops_[dim - index2 - 1].mo;
            taskList[taskIndex] = std::make_pair(mo1, mo2);
            ++taskIndex;
        }
    }

    // big-small pair
    // for (int index1 = 0; index1 < dim - 1; ++index1) {
    //     const int max_index2 = dim - index1 - 1;
    //     for (int index2 = 0; index2 < max_index2; ++index2) {
    //         const index_type mo1 = this->groupMoPops_[index1 + index2].mo;
    //         const index_type mo2 = this->groupMoPops_[dim - index2 - 1].mo;
    //         taskList[taskIndex] = std::make_pair(mo1, mo2);
    //         ++taskIndex;
    //     }
    // }

    assert(taskIndex == pairs);

    return taskList;
}
