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

#include "DfLocalize_Parallel.h"

#include <map>
#include <set>

#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"

DfLocalize_Parallel::DfLocalize_Parallel(TlSerializeData* pPdfParam) : DfLocalize(pPdfParam) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.setProcID(rComm.getRank());
}

DfLocalize_Parallel::~DfLocalize_Parallel() {
}

void DfLocalize_Parallel::initialize() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    this->G_ = 0.0;

    // output information
    this->log_.info(TlUtils::format("number of Atoms = %d", this->m_nNumOfAtoms));
    this->log_.info(TlUtils::format("number of AOs   = %d", this->m_nNumOfAOs));
    this->log_.info(TlUtils::format("number of MOs   = %d", this->m_nNumOfMOs));
    this->log_.info(TlUtils::format("max iteration = %d", this->maxIteration_));
    this->log_.info(TlUtils::format("threshold = %10.5e", this->threshold_));

    this->makeGroup();

    // load S matrix
    bool isLoadS = false;
    const std::string S_path = this->getSpqMatrixPath();
    if (rComm.isMaster()) {
        isLoadS = this->getSMatrix(&(this->S_));
    }
    rComm.broadcast(isLoadS);

    if (isLoadS == false) {
        if (rComm.isMaster()) {
            const std::string path = this->getSpqMatrixPath();
            this->log_.critical(TlUtils::format("could not load: %s", S_path.c_str()));
        }
        abort();
    } else {
        rComm.broadcast(&(this->S_));
    }
}

void DfLocalize_Parallel::exec() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->initialize();

    TlDenseGeneralMatrix_Lapack C;
    if (rComm.isMaster()) {
        this->getCMatrix(&C);
    }

    rComm.broadcast(this->lo_iteration_);
    const int maxIteration = this->maxIteration_;
    for (int i = this->lo_iteration_ + 1; i < maxIteration; ++i) {
        const double sumDeltaG = this->localize(&C);

        const std::string msg = TlUtils::format("%d th: G=%10.5e, delta_G=%10.5e", i, this->G_, sumDeltaG);
        this->log_.info(msg);
        if (rComm.isMaster()) {
            std::cout << msg << std::endl;
        }

        if (rComm.isMaster()) {
            DfObject::saveCloMatrix(RUN_RKS, i, C);
        }

        if (sumDeltaG < this->threshold_) {
            if (rComm.isMaster()) {
                (*(this->pPdfParam_))["lo/num_of_iterations"] = i;
                (*(this->pPdfParam_))["lo/satisfied"] = "yes";
            }
            break;
        }
    }
}

double DfLocalize_Parallel::localize(TlDenseGeneralMatrix_Lapack* pC) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const index_type numOfAOs = this->m_nNumOfAOs;

    if (rComm.isMaster() == true) {
        this->makeQATable(*pC);
        this->makeJobList();

        JobItem jobItem;                         // dummy
        (void)this->getJobItem(&jobItem, true);  // initialize
    }

    enum { REQUEST_JOB = 0, SEND_RESULTS = 1 };
    enum { FINISHED_JOB = 0, ASSIGNED_JOB = 1, WAIT = 2 };

    std::vector<double> CpiCpj(numOfAOs * 2);

    double sumDeltaG = 0.0;
    if (rComm.isMaster() == true) {
        int jobRequest = 0;
        int src = 0;

        bool isEscapeLoop = false;
        int numOfFinishedProc = 0;
        std::map<int, JobItem> rank_job_db;

        TlDenseVector_Lapack vec_i(numOfAOs);
        TlDenseVector_Lapack vec_j(numOfAOs);
        do {
            rComm.receiveDataFromAnySource(jobRequest, &src);
            // std::cerr << TlUtils::format("[0] recv %d from=%d", jobRequest, src) << std::endl;
            switch (jobRequest) {
                case REQUEST_JOB:
                    // slaveが仕事を要求
                    {
                        JobItem jobItem;
                        const bool hasJob = this->getJobItem(&jobItem);
                        if (hasJob) {
                            rank_job_db[src] = jobItem;
                            const index_type orb_i = jobItem.orb_i;
                            const index_type orb_j = jobItem.orb_j;

                            if ((orb_i >= 0) && (orb_j >= 0)) {
                                // assigned job
                                this->lockMO(orb_i);
                                this->lockMO(orb_j);
                                TlDenseVector_Lapack vec_i = pC->getColVector(orb_i);
                                TlDenseVector_Lapack vec_j = pC->getColVector(orb_j);

                                // std::cerr << TlUtils::format("[0] send ASSIGNED_JOB to=%d", src) << std::endl;
                                rComm.sendData(ASSIGNED_JOB, src);
                                // std::cerr << TlUtils::format("[0] send DATA to=%d", src) << std::endl;
                                std::copy(vec_i.data(), vec_i.data() + numOfAOs, &(CpiCpj[0]));
                                std::copy(vec_j.data(), vec_j.data() + numOfAOs, &(CpiCpj[numOfAOs]));

                                rComm.sendDataX(&(CpiCpj[0]), numOfAOs * 2, src);
                                // std::cerr << "[0] end JOB assign" << std::endl;
                            } else {
                                // wait
                                { rComm.sendData(WAIT, src); }
                            }
                        } else {
                            // finished
                            {
                                rComm.sendData(FINISHED_JOB, src);
                                ++numOfFinishedProc;
                                if (numOfFinishedProc == rComm.getNumOfProc() - 1) {
                                    isEscapeLoop = true;
                                }
                            }
                        }
                    }
                    break;

                case SEND_RESULTS:
                    // slaveが結果を返す
                    {
                        std::map<int, JobItem>::const_iterator it = rank_job_db.find(src);
                        if (it == rank_job_db.end()) {
                            this->log_.critical("program error.");
                        }

                        rComm.receiveDataX(&(CpiCpj[0]), numOfAOs * 2, src);
                        const index_type orb_i = it->second.orb_i;
                        const index_type orb_j = it->second.orb_j;
                        rank_job_db.erase(src);

                        // double deltaG = 0.0;
                        // rComm.receiveData(deltaG, src);
                        // sumDeltaG += deltaG;

                        // 行列の格納
                        for (index_type row = 0; row < numOfAOs; ++row) {
                            pC->set(row, orb_i, CpiCpj[row]);
                            pC->set(row, orb_j, CpiCpj[numOfAOs + row]);
                        }

                        // 行列ロックの解除
                        this->unlockMO(orb_i);
                        this->unlockMO(orb_j);
                    }
                    break;

                default:
                    // something wrong...
                    abort();
                    break;
            }
        } while (isEscapeLoop == false);
    } else {
        // for Slave
        const int root = 0;
        int request = REQUEST_JOB;
        int hasJob = 0;

        rComm.sendData(request, root);
        // std::cerr << TlUtils::format("[%d] send request", rank) << std::endl;
        rComm.receiveData(hasJob, root);
        // std::cerr << TlUtils::format("[%d] recv JOB %d", rank, hasJob) << std::endl;

        TlDenseVector_Lapack Cpi(numOfAOs);
        TlDenseVector_Lapack Cpj(numOfAOs);
        while (hasJob != FINISHED_JOB) {
            if (hasJob == ASSIGNED_JOB) {
                const index_type numOfAOs = this->m_nNumOfAOs;

                rComm.receiveDataX(&(CpiCpj[0]), numOfAOs * 2, root);
                // std::cerr << TlUtils::format("[%d] recv DATA", rank) << std::endl;
                std::copy(&(CpiCpj[0]), &(CpiCpj[0]) + numOfAOs, Cpi.data());
                std::copy(&(CpiCpj[0]) + numOfAOs, &(CpiCpj[0]) + numOfAOs * 2, Cpj.data());

                // std::cerr << Cpi << std::endl;
                // std::cerr << Cpj << std::endl;

                double A_ij = 0.0;
                double B_ij = 0.0;
                this->calcQA_ij(Cpi, Cpj, &A_ij, &B_ij);
                const double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
                const double deltaG = A_ij + normAB;

                if (std::fabs(deltaG) > 1.0E-16) {
                    sumDeltaG += deltaG;
                    // std::cerr << "sumDeltaG = " << sumDeltaG << std::endl;
                    TlDenseGeneralMatrix_Lapack rot(2, 2);
                    this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
                    this->rotateVectors(&Cpi, &Cpj, rot);
                }

                request = SEND_RESULTS;
                std::copy(Cpi.data(), Cpi.data() + numOfAOs, &(CpiCpj[0]));
                std::copy(Cpj.data(), Cpj.data() + numOfAOs, &(CpiCpj[numOfAOs]));

                rComm.sendData(request, root);
                rComm.sendDataX(&(CpiCpj[0]), numOfAOs * 2, root);
                // rComm.sendData(deltaG, root);
            } else if (hasJob == WAIT) {
                TlTime::sleep(4000);  // wait 4000 ms.
            }

            request = REQUEST_JOB;
            rComm.sendData(request, root);
            rComm.receiveData(hasJob, root);
        }
    }

    // summarize
    rComm.allReduce_SUM(sumDeltaG);
    // std::cout << TlUtils::format("[%d] sumDeltaG=%f", rank, sumDeltaG) << std::endl;
    // bool isBreak = false;
    // if (rComm.isMaster() == true) {
    // std::cout << TlUtils::format("itr: %d sum of delta_g: %10.5e\n", num_iteration + 1, sumDeltaG);
    // this->C_.save(TlUtils::format("./lo_Work/fl_Mtr_C.lo.occu.rks%d",
    // num_iteration +1));
    // DfObject::saveCloMatrix(RUN_RKS, num_iteration + 1, this->C_);

    // if (sumDeltaG < this->threshold_) {
    // std::cout << "number of iteration: " << num_iteration + 1 << std::endl;
    // (*(this->pPdfParam_))["lo/num_of_iterations"] = num_iteration + 1;
    // (*(this->pPdfParam_))["lo/satisfied"] = "yes";
    // isBreak = true;
    // }
    // }

    // rComm.broadcast(isBreak);
    // if (isBreak == true) {
    //     break;
    // }

    return sumDeltaG;
}

// void DfLocalize_Parallel::localize(const std::string& inputCMatrixPath, const int restartIteration) {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const std::size_t numOfAOs = this->m_nNumOfAOs;

//     if (rComm.isMaster() == true) {
//         std::cout << TlUtils::format("number of Atoms = %ld", this->m_nNumOfAtoms) << std::endl;
//         std::cout << TlUtils::format("number of AOs   = %ld", this->m_nNumOfAOs) << std::endl;
//         std::cout << TlUtils::format("number of MOs   = %ld", this->m_nNumOfMOs) << std::endl;
//     }
//     this->makeGroup();

//     this->S_.resize(this->m_nNumOfAOs);
//     // this->C_.resize(this->m_nNumOfAOs, this->m_nNumOfMOs);
//     if (rComm.isMaster() == true) {
//         const std::string S_path = this->getSpqMatrixPath();
//         if (this->S_.load(S_path) == false) {
//             std::cerr << "could not load: " << S_path << std::endl;
//             abort();
//         }

//         std::string CMatrixPath = inputCMatrixPath;
//         if (CMatrixPath.empty() == true) {
//             CMatrixPath = this->getCMatrixPath(DfObject::RUN_RKS, this->m_nIteration);
//         }
//         if (this->C_.load(CMatrixPath) == false) {
//             std::cerr << "could not load: " << CMatrixPath << std::endl;
//             abort();
//         }
//     }
//     rComm.broadcast(&(this->S_));

//     enum { REQUEST_JOB = 0, SEND_RESULTS = 1 };

//     enum { FINISHED_JOB = 0, ASSIGNED_JOB = 1, WAIT = 2 };

//     if (rComm.isMaster() == true) {
//         std::cerr << "begin loop." << std::endl;
//     }

//     int startIteration = std::max(0, restartIteration);
//     if (rComm.isMaster()) {
//         if (restartIteration > 0) {
//             std::cerr << TlUtils::format("restart from iteration=%d", startIteration + 1) << std::endl;
//             this->C_ = DfObject::getCloMatrix<TlDenseGeneralMatrix_Lapack>(DfObject::RUN_RKS, startIteration + 1);
//         }
//     }
//     rComm.broadcast(startIteration);
//     rComm.broadcast(&(this->C_));

//     const int maxIteration = this->maxIteration_;
//     for (int num_iteration = startIteration; num_iteration < maxIteration; ++num_iteration) {
//         if (rComm.isMaster() == true) {
//             this->makeQATable();
//             this->makeJobList();
//             JobItem jobItem;  // dummy
//             (void)this->getJobItem(&jobItem, true);
//         }

//         double sumDeltaG = 0.0;
//         if (rComm.isMaster() == true) {
//             int jobRequest = 0;
//             int src = 0;

//             JobItem jobItem;
//             bool isEscapeLoop = false;
//             int numOfFinishedProc = 0;
//             std::size_t orb_i = 0;
//             std::size_t orb_j = 0;
//             TlDenseVector_Lapack vec_i(numOfAOs);
//             TlDenseVector_Lapack vec_j(numOfAOs);
//             do {
//                 rComm.receiveDataFromAnySource(jobRequest, &src);
//                 // std::cerr << TlUtils::format("[0] recv %d from=%d",
//                 // jobRequest, src)
//                 // << std::endl;
//                 switch (jobRequest) {
//                     case REQUEST_JOB:
//                         // slaveが仕事を要求
//                         {
//                             const int hasJob = this->getJobItem(&jobItem);

//                             switch (hasJob) {
//                                 case 0:
//                                     // finished
//                                     {
//                                         rComm.sendData(FINISHED_JOB, src);
//                                         ++numOfFinishedProc;
//                                         if (numOfFinishedProc == rComm.getNumOfProc() - 1) {
//                                             isEscapeLoop = true;
//                                         }
//                                     }
//                                     break;
//                                 case 1:
//                                     // assigned job
//                                     {
//                                         rComm.sendData(ASSIGNED_JOB, src);
//                                         const std::size_t orb_i = jobItem.orb_i;
//                                         const std::size_t orb_j = jobItem.orb_j;
//                                         rComm.sendData(orb_i, src);
//                                         rComm.sendData(orb_j, src);
//                                         TlDenseVector_Lapack vec_i =
//                                         this->C_.getColVector<TlDenseVector_Lapack>(orb_i); TlDenseVector_Lapack
//                                         vec_j = this->C_.getColVector<TlDenseVector_Lapack>(orb_j);
//                                         rComm.sendData(vec_i, src);
//                                         rComm.sendData(vec_j, src);
//                                     }
//                                     break;
//                                 case 2:
//                                     // wait
//                                     { rComm.sendData(WAIT, src); }
//                                     break;

//                                 default:
//                                     // something wrong
//                                     abort();
//                                     break;
//                             }
//                         }
//                         break;

//                     case SEND_RESULTS:
//                         // slaveが結果を返す
//                         {
//                             rComm.receiveData(orb_i, src);
//                             rComm.receiveData(orb_j, src);
//                             rComm.receiveData(vec_i, src);
//                             rComm.receiveData(vec_j, src);
//                             // std::cerr << TlUtils::format("[0] recv job
//                             // from=%d", src) << std::endl;

//                             // 行列の格納
//                             for (std::size_t row = 0; row < numOfAOs; ++row) {
//                                 this->C_.set(row, orb_i, vec_i.get(row));
//                                 this->C_.set(row, orb_j, vec_j.get(row));
//                             }

//                             // 行列ロックの解除
//                             const std::size_t index_i = orb_i - this->startOrb_;
//                             const std::size_t index_j = orb_j - this->startOrb_;
//                             this->jobOccupiedOrb_[index_i] = false;
//                             this->jobOccupiedOrb_[index_j] = false;
//                         }
//                         break;

//                     default:
//                         // something wrong...
//                         abort();
//                         break;
//                 }
//             } while (isEscapeLoop == false);
//         } else {
//             // for Slave
//             const int root = 0;
//             int request = REQUEST_JOB;
//             int hasJob = 0;

//             rComm.sendData(request, root);
//             rComm.receiveData(hasJob, root);

//             std::size_t orb_i = 0;
//             std::size_t orb_j = 0;
//             TlDenseVector_Lapack vec_i(numOfAOs);
//             TlDenseVector_Lapack vec_j(numOfAOs);
//             TlDenseGeneralMatrix_Lapack rot(2, 2);
//             while (hasJob != FINISHED_JOB) {
//                 if (hasJob == ASSIGNED_JOB) {
//                     const std::size_t numOfAOs = this->m_nNumOfAOs;
//                     rComm.receiveData(orb_i, root);
//                     rComm.receiveData(orb_j, root);
//                     rComm.receiveData(vec_i, root);
//                     rComm.receiveData(vec_j, root);
//                     // std::cerr << TlUtils::format("[%d] recv job",
//                     // rComm.getRank()) << std::endl;

//                     assert(vec_i.getSize() == static_cast<TlVectorAbstract::size_type>(numOfAOs));
//                     assert(vec_j.getSize() == static_cast<TlVectorAbstract::size_type>(numOfAOs));
//                     for (std::size_t row = 0; row < numOfAOs; ++row) {
//                         this->C_.set(row, orb_i, vec_i.get(row));
//                         this->C_.set(row, orb_j, vec_j.get(row));
//                     }

//                     double A_ij = 0.0;
//                     double B_ij = 0.0;
//                     this->calcQA(orb_i, orb_j, &A_ij, &B_ij);
//                     double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
//                     double deltaG = A_ij + normAB;
//                     if (std::fabs(deltaG) > 1.0E-16) {
//                         sumDeltaG += deltaG;
//                         this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
//                         this->rotateCmatrix(orb_i, orb_j, rot);
//                     }

//                     vec_i = this->C_.getColVector<TlDenseVector_Lapack>(orb_i);
//                     vec_j = this->C_.getColVector<TlDenseVector_Lapack>(orb_j);

//                     request = SEND_RESULTS;
//                     rComm.sendData(request, root);
//                     rComm.sendData(orb_i, root);
//                     rComm.sendData(orb_j, root);
//                     rComm.sendData(vec_i, root);
//                     rComm.sendData(vec_j, root);
//                     // std::cerr << TlUtils::format("[%d] send job",
//                     // rComm.getRank()) << std::endl;
//                 } else if (hasJob == WAIT) {
//                     TlTime::sleep(1000);  // wait 4000 ms.
//                 }

//                 request = REQUEST_JOB;
//                 rComm.sendData(request, root);
//                 rComm.receiveData(hasJob, root);
//             }
//         }

//         // summarize
//         rComm.allReduce_SUM(sumDeltaG);
//         bool isBreak = false;
//         if (rComm.isMaster() == true) {
//             std::cout << TlUtils::format("itr: %d sum of delta_g: %10.5e\n", num_iteration + 1, sumDeltaG);
//             // this->C_.save(TlUtils::format("./lo_Work/fl_Mtr_C.lo.occu.rks%d",
//             // num_iteration +1));
//             DfObject::saveCloMatrix(RUN_RKS, num_iteration + 1, this->C_);

//             if (sumDeltaG < this->threshold_) {
//                 std::cout << "number of iteration: " << num_iteration + 1 << std::endl;
//                 (*(this->pPdfParam_))["lo/num_of_iterations"] = num_iteration + 1;
//                 (*(this->pPdfParam_))["lo/satisfied"] = "yes";
//                 isBreak = true;
//             }
//         }

//         rComm.broadcast(isBreak);
//         if (isBreak == true) {
//             break;
//         }
//     }
// }

bool DfLocalize_Parallel::getJobItem(DfLocalize::JobItem* pJob, bool isInitialized) {
    assert(pJob != NULL);

    static std::list<JobItem>::iterator it;
    static std::set<index_type> lockMOs;

    if (isInitialized == true) {
        it = this->jobList_.begin();

        // for parallel operation
        this->lockMOs_.clear();

        return true;
    }

    bool answer = false;
    bool locked = false;
    std::list<JobItem>::iterator itEnd = this->jobList_.end();
    for (it = this->jobList_.begin(); it != itEnd; ++it) {
        answer = true;
        if (this->isLockedMO(it->orb_i, it->orb_j)) {
            locked = true;
        } else {
            *pJob = *it;
            it = this->jobList_.erase(it);
            break;
        }
    }

    if (locked) {
        pJob->orb_i = -1;
        pJob->orb_j = -1;
    }

    return answer;
}

void DfLocalize_Parallel::lockMO(const index_type mo) {
    this->lockMOs_.insert(mo);
}

void DfLocalize_Parallel::unlockMO(const index_type mo) {
    this->lockMOs_.erase(mo);
}

bool DfLocalize_Parallel::isLockedMO(const index_type mo1, const index_type mo2) const {
    bool answer = true;

    std::set<index_type>::const_iterator itEnd = this->lockMOs_.end();
    if (this->lockMOs_.find(mo1) == itEnd) {
        if (this->lockMOs_.find(mo2) == itEnd) {
            answer = false;
        }
    }

    return answer;
}
