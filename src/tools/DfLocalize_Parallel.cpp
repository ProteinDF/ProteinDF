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
#include "TlCommunicate.h"
#include "TlUtils.h"
#include "TlTime.h"

DfLocalize_Parallel::DfLocalize_Parallel(TlSerializeData* pPdfParam)
    : DfLocalize(pPdfParam)
{
}


DfLocalize_Parallel::~DfLocalize_Parallel()
{
}


void DfLocalize_Parallel::localize(const std::string& inputCMatrixPath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const std::size_t numOfAOs = this->m_nNumOfAOs;

    if (rComm.isMaster() == true) {
        std::cout << TlUtils::format("number of Atoms = %ld", this->m_nNumOfAtoms) << std::endl;
        std::cout << TlUtils::format("number of AOs   = %ld", this->m_nNumOfAOs) << std::endl;
        std::cout << TlUtils::format("number of MOs   = %ld", this->m_nNumOfMOs) << std::endl;
    }
    this->makeGroup();

    this->S_.resize(this->m_nNumOfAOs);
    this->C_.resize(this->m_nNumOfAOs, this->m_nNumOfMOs);
    if (rComm.isMaster() == true) {
        const std::string S_path = this->getSpqMatrixPath();
        if (this->S_.load(S_path) == false) {
            std::cerr << "could not load: " << S_path << std::endl;
            abort();
        }

        std::string CMatrixPath = inputCMatrixPath;
        if (CMatrixPath.empty() == true) {
            CMatrixPath = this->getCMatrixPath(DfObject::RUN_RKS, this->m_nIteration);
        }
        if (this->C_.load(CMatrixPath) == false) {
            std::cerr << "could not load: " << CMatrixPath << std::endl;
            abort();
        }
    }
    rComm.broadcast(this->S_);

    enum {
        REQUEST_JOB = 0,
        SEND_RESULTS = 1
    };

    enum {
        FINISHED_JOB = 0,
        ASSIGNED_JOB = 1,
        WAIT = 2
    };

    if (rComm.isMaster() == true) {
        std::cerr << "begin loop." << std::endl;
    }

    const int maxIteration = this->maxIteration_;
    for (int num_iteration = 0; num_iteration < maxIteration; ++num_iteration) {
        if (rComm.isMaster() == true) {
            this->makeQATable();
            this->makeJobList();
            JobItem jobItem; // dummy
            (void)this->getJobItem(&jobItem, true);
        }

        double sumDeltaG = 0.0;
        if (rComm.isMaster() == true) {
            int jobRequest = 0;
            int src = 0;

            JobItem jobItem;
            bool isEscapeLoop = false;
            int numOfFinishedProc = 0;
            std::size_t orb_i = 0;
            std::size_t orb_j = 0;
            TlVector vec_i(numOfAOs);
            TlVector vec_j(numOfAOs);
            do {
                rComm.receiveDataFromAnySource(jobRequest, &src);
                //std::cerr << TlUtils::format("[0] recv %d from=%d", jobRequest, src) << std::endl;
                switch (jobRequest) {
                case REQUEST_JOB:
                    // slaveが仕事を要求
                    {
                        const int hasJob = this->getJobItem(&jobItem);
                        
                        switch (hasJob) {
                        case 0:
                            // finished
                            {
                                rComm.sendData(FINISHED_JOB, src);
                                ++numOfFinishedProc;
                                if (numOfFinishedProc == rComm.getNumOfProc() -1) {
                                    isEscapeLoop = true;
                                }
                            }
                            break;
                        case 1:
                            // assigned job
                            {
                                rComm.sendData(ASSIGNED_JOB, src);
                                const std::size_t orb_i = jobItem.orb_i;
                                const std::size_t orb_j = jobItem.orb_j;
                                rComm.sendData(orb_i, src);
                                rComm.sendData(orb_j, src);
                                TlVector vec_i = this->C_.getColVector(orb_i);
                                TlVector vec_j = this->C_.getColVector(orb_j);
                                rComm.sendData(vec_i, src);
                                rComm.sendData(vec_j, src);
                            }
                            break;
                        case 2:
                            // wait
                            {
                                rComm.sendData(WAIT, src);
                            }
                            break;
                            
                        default:
                            // something wrong
                            abort();
                            break;
                        }
                    }
                break;

                case SEND_RESULTS:
                    // slaveが結果を返す
                    {
                        rComm.receiveData(orb_i, src);
                        rComm.receiveData(orb_j, src);
                        rComm.receiveData(vec_i, src);
                        rComm.receiveData(vec_j, src);
                        //std::cerr << TlUtils::format("[0] recv job from=%d", src) << std::endl;
                        
                        // 行列の格納
                        for (std::size_t row = 0; row < numOfAOs; ++row) {
                            this->C_.set(row, orb_i, vec_i[row]);
                            this->C_.set(row, orb_j, vec_j[row]);
                        }
                        
                        // 行列ロックの解除
                        const std::size_t index_i = orb_i - this->startOrb_;
                        const std::size_t index_j = orb_j - this->startOrb_;
                        this->jobOccupiedOrb_[index_i] = false;
                        this->jobOccupiedOrb_[index_j] = false;
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
            rComm.receiveData(hasJob, root);

            std::size_t orb_i = 0;
            std::size_t orb_j = 0;
            TlVector vec_i(numOfAOs);
            TlVector vec_j(numOfAOs);
            TlMatrix rot(2, 2);
            while (hasJob != FINISHED_JOB) {
                if (hasJob == ASSIGNED_JOB) {
                    const std::size_t numOfAOs = this->m_nNumOfAOs;
                    rComm.receiveData(orb_i, root);
                    rComm.receiveData(orb_j, root);
                    rComm.receiveData(vec_i, root);
                    rComm.receiveData(vec_j, root);
                    //std::cerr << TlUtils::format("[%d] recv job", rComm.getRank()) << std::endl;

                    assert(vec_i.getSize() == static_cast<TlVector::size_type>(numOfAOs));
                    assert(vec_j.getSize() == static_cast<TlVector::size_type>(numOfAOs));
                    for (std::size_t row = 0; row < numOfAOs; ++row) {
                        this->C_.set(row, orb_i, vec_i[row]);
                        this->C_.set(row, orb_j, vec_j[row]);
                    }

                    double A_ij = 0.0;
                    double B_ij = 0.0;
                    this->calcQA(orb_i, orb_j, &A_ij, &B_ij);
                    double normAB = std::sqrt(A_ij * A_ij + B_ij * B_ij);
                    double deltaG = A_ij + normAB;
                    if (std::fabs(deltaG) > 1.0E-16) {
                        sumDeltaG += deltaG;
                        this->getRotatingMatrix(A_ij, B_ij, normAB, &rot);
                        this->rotateCmatrix(orb_i, orb_j, rot);
                    }

                    vec_i = this->C_.getColVector(orb_i);
                    vec_j = this->C_.getColVector(orb_j);

                    request = SEND_RESULTS;
                    rComm.sendData(request, root);
                    rComm.sendData(orb_i, root);
                    rComm.sendData(orb_j, root);
                    rComm.sendData(vec_i, root);
                    rComm.sendData(vec_j, root);
                    //std::cerr << TlUtils::format("[%d] send job", rComm.getRank()) << std::endl;
                } else if (hasJob == WAIT) {
                    TlTime::sleep(1000); // wait 4000 ms.
                }

                request = REQUEST_JOB;
                rComm.sendData(request, root);
                rComm.receiveData(hasJob, root);
            }
        }

        // summarize
        rComm.allReduce_SUM(sumDeltaG);
        bool isBreak = false;
        if (rComm.isMaster() == true) {
            std::cout << TlUtils::format("itr: %d sum of delta_g: %10.5e\n", num_iteration +1, sumDeltaG);
            // this->C_.save(TlUtils::format("./lo_Work/fl_Mtr_C.lo.occu.rks%d", num_iteration +1));
            DfObject::saveCloMatrix(RUN_RKS, num_iteration +1, this->C_);

            if (sumDeltaG < this->threshold_) {
                std::cout << "number of iteration: " << num_iteration +1 << std::endl;
                (*(this->pPdfParam_))["lo/num_of_iterations"] = num_iteration +1;
                (*(this->pPdfParam_))["lo/satisfied"] = "yes";
                isBreak = true;
            }
        }

        rComm.broadcast(isBreak);
        if (isBreak == true) {
            break;
        }
    }
}


// @ret 0 assigned no job because the task has been finished
// @ret 1 assigned job to "pJob"
// @ret 2 please wait because the job to assign is conflict.
int DfLocalize_Parallel::getJobItem(DfLocalize::JobItem* pJob, bool isInitialized)
{
    assert(pJob != NULL);

    // static std::size_t jobIndex = 0;
    const std::size_t maxJobIndex = this->jobList_.size();

    if (isInitialized == true) {
        // jobIndex = 0;
        //maxJobIndex = this->jobList_.size();

        // for parallel operation
        this->jobFinishedList_.clear();
        this->jobFinishedList_.resize(maxJobIndex, false);
        this->jobOccupiedOrb_.clear();
        this->jobOccupiedOrb_.resize(this->orb_QA_table_.size(), false);

        return true;
    }

    int answer = 0;
    for (std::size_t i = 0; i < maxJobIndex; ++i) {
        if (this->jobFinishedList_[i] == true) {
            continue;
        }

        const JobItem item = this->jobList_[i];
        const std::size_t index_i = item.orb_i - this->startOrb_;
        const std::size_t index_j = item.orb_j - this->startOrb_;

        if ((this->jobOccupiedOrb_[index_i] == false) &&
                (this->jobOccupiedOrb_[index_j] == false)) {
            this->jobFinishedList_[i] = true;
            this->jobOccupiedOrb_[index_i] = true;
            this->jobOccupiedOrb_[index_j] = true;
            *pJob = item;
            answer = 1; // assign

            break;
        }
    }

    // check
    if (answer == 0) {
        for (std::size_t i = 0; i < maxJobIndex; ++i) {
            if (this->jobFinishedList_[i] != true) {
                answer = 2; // wait
                break;
            }
        }
    }

    return answer;
}

