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

#include <cassert>
#include <map>

#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"

DfLocalize_Parallel::DfLocalize_Parallel(TlSerializeData* pPdfParam)
    : DfLocalize(pPdfParam) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.setProcID(rComm.getRank());
}

DfLocalize_Parallel::~DfLocalize_Parallel() {
}

void DfLocalize_Parallel::initialize() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

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

    // load S matrix
    this->log_.info("load S matrix");
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

    TlDenseGeneralMatrix_Lapack C;
    if (rComm.isMaster()) {
        this->getCMatrix(&C);
    }

    this->log_.info("make group");
    this->makeGroup();

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
    const int rank = rComm.getRank();
    const int numOfProcs = rComm.getNumOfProcs();

    const index_type numOfMOs = this->endOrb_ - this->startOrb_ + 1;
    const int blocks = numOfProcs * 2;
    const index_type MOsPerBlock = numOfMOs / blocks;

    this->log_.info(TlUtils::format("numOfMOs: %d", numOfMOs));
    this->log_.info(TlUtils::format("blocks: %d", blocks));
    this->log_.info(TlUtils::format("MOsPerBlock: %d", MOsPerBlock));

    // const index_type numOfAOs = this->m_nNumOfAOs;

    JobItem jobItem;
    if (rComm.isMaster() == true) {
        this->G_ = this->calcG(*pC, this->startOrb_, this->endOrb_);
        this->makeJobList(MOsPerBlock);

        (void)this->getJobItem(&jobItem, true);
    }

    enum { REQUEST_JOB = 0,
           SEND_RESULTS = 1 };
    enum { FINISHED_JOB = 0,
           ASSIGNED_JOB = 1,
           WAIT = 2 };
    const int msgLength = 3;
    std::vector<int> msg(msgLength);

    double sumDeltaG = 0.0;
    if (rComm.isMaster() == true) {
        this->initLockBlock(blocks + 1);
        int taskRequest = 0;
        int src = 0;

        int numOfFinishedProc = 0;
        std::map<int, JobItem> rank_job_db;

        bool hasTask = false;
        JobItem reservedJobItem;
        bool isEscapeLoop = false;
        while (isEscapeLoop == false) {
            rComm.receiveDataFromAnySource(taskRequest, &src);
            this->log_.info(TlUtils::format("[0] recv %d from=%d", taskRequest, src));
            switch (taskRequest) {
                case REQUEST_JOB: {
                    if (hasTask) {
                        jobItem = reservedJobItem;
                    } else {
                        hasTask = this->getJobItem(&jobItem);
                    }
                    if (hasTask) {
                        rank_job_db[src] = jobItem;
                        const int block1 = jobItem.first;
                        const int block2 = jobItem.second;

                        this->log_.info(TlUtils::format("[0] check block: %d, %d", block1, block2));
                        if (this->isLockedBlock(block1, block2) == false) {
                            this->log_.info(TlUtils::format("[0] check OK: %d, %d", block1, block2));
                            // assigned task
                            this->lockBlock(block1);
                            this->lockBlock(block2);
                            this->log_.info(TlUtils::format("[0] lock: %d, %d", block1, block2));
                            TlDenseGeneralMatrix_Lapack blockC;
                            this->getBlockCMatrix(*pC, MOsPerBlock, block1, block2, &blockC);

                            this->log_.info(TlUtils::format("[0] send ASSIGNED_JOB to=%d", src));
                            msg[0] = ASSIGNED_JOB;
                            msg[1] = block1;
                            msg[2] = block2;
                            rComm.sendDataX(&(msg[0]), msgLength, src);
                            rComm.sendData(blockC, src);
                            // std::cerr << "[0] end JOB assign" << std::endl;
                            hasTask = false;
                        } else {
                            // wait
                            this->log_.info(TlUtils::format("[0] send WAIT to %d", src));
                            msg[0] = WAIT;
                            msg[1] = 0;
                            msg[2] = 0;
                            rComm.sendDataX(&(msg[0]), msgLength, src);
                            reservedJobItem = jobItem;
                        }
                    } else {
                        // finished
                        {
                            this->log_.info(TlUtils::format("[0] send FINISHED_JOB to=%d", src));
                            msg[0] = FINISHED_JOB;
                            msg[1] = 0;
                            msg[2] = 0;
                            rComm.sendDataX(&(msg[0]), msgLength, src);
                            ++numOfFinishedProc;
                            if (numOfFinishedProc == rComm.getNumOfProc() - 1) {
                                isEscapeLoop = true;
                            }
                        }
                    }
                } break;

                case SEND_RESULTS: {
                    this->log_.info(TlUtils::format("[0] recv SEND_RESULTS from %d", src));
                    std::map<int, JobItem>::const_iterator it = rank_job_db.find(src);
                    if (it == rank_job_db.end()) {
                        this->log_.critical("program error.");
                    }

                    TlDenseGeneralMatrix_Lapack blockC;
                    rComm.receiveData(blockC, src);
                    const int block1 = it->second.first;
                    const int block2 = it->second.second;
                    rank_job_db.erase(src);

                    this->log_.debug(TlUtils::format("[0] set results from %d: C(%d x %d), block: %d, %d ", src,
                                                     blockC.getNumOfRows(), blockC.getNumOfCols(), block1, block2));
                    this->setBlockCMatrix(MOsPerBlock, block1, block2, blockC, pC);

                    this->unlockBlock(block1);
                    this->unlockBlock(block2);
                } break;

                default:
                    // something wrong...
                    abort();
                    break;
            }
        }
    } else {
        // for Slave
        const int root = 0;
        int request = REQUEST_JOB;
        int hasTask = 0;

        rComm.sendData(request, root);
        // std::cerr << TlUtils::format("[%d] send request", rank) << std::endl;
        rComm.receiveDataX(&(msg[0]), msgLength, root);
        this->log_.debug(TlUtils::format("[%d] recv JOB %d", rank, hasTask));

        while (msg[0] != FINISHED_JOB) {
            if (msg[0] == ASSIGNED_JOB) {
                const int block1 = msg[1];
                const int block2 = msg[2];
                this->log_.debug(TlUtils::format("[%d] ASSIGNED_JOB %d, %d", rank, block1, block2));

                TlDenseGeneralMatrix_Lapack blockC;
                rComm.receiveData(blockC, root);
                this->log_.debug(
                    TlUtils::format("[%d] recv matrix: %d, %d", rank, blockC.getNumOfRows(), blockC.getNumOfCols()));

                const index_type numOfCols = blockC.getNumOfCols();
                index_type startMO1 = 0;
                index_type endMO1 = numOfCols - 1;
                index_type startMO2 = 0;
                index_type endMO2 = numOfCols - 1;
                if (numOfCols > MOsPerBlock) {
                    endMO1 = MOsPerBlock - 1;
                    startMO2 = MOsPerBlock;
                }

                this->log_.debug(
                    TlUtils::format("[%d] run localize: %d, %d, %d, %d", rank, startMO1, endMO1, startMO2, endMO2));
                sumDeltaG = DfLocalize::localize_core(&blockC, startMO1, endMO1, startMO2, endMO2);

                request = SEND_RESULTS;
                rComm.sendData(request, root);
                rComm.sendData(blockC, root);
            } else if (msg[0] == WAIT) {
                TlTime::sleep(4000);  // wait 4000 ms.
            }

            request = REQUEST_JOB;
            rComm.sendData(request, root);
            rComm.receiveDataX(&(msg[0]), msgLength, root);
        }
    }

    // summarize
    rComm.allReduce_SUM(sumDeltaG);

    return sumDeltaG;
}

void DfLocalize_Parallel::getBlockCMatrix(const TlDenseGeneralMatrix_Lapack& C, const index_type MOsPerBlock,
                                          const int block1, const int block2, TlDenseGeneralMatrix_Lapack* pBlockC) {
    const index_type numOfAOs = this->m_nNumOfAOs;

    {
        const index_type startMO1 = block1 * MOsPerBlock + this->startOrb_;
        const index_type endMO1 = std::min((block1 + 1) * MOsPerBlock - 1, this->endOrb_);

        C.block(0, startMO1, numOfAOs, endMO1 - startMO1 + 1, pBlockC);
    }

    if (block1 != block2) {
        TlDenseGeneralMatrix_Lapack tmpBlockC;
        const index_type startMO2 = block2 * MOsPerBlock + this->startOrb_;
        const index_type endMO2 = std::min((block2 + 1) * MOsPerBlock - 1, this->endOrb_);
        C.block(0, startMO2, numOfAOs, endMO2 - startMO2 + 1, &tmpBlockC);

        const index_type block1endCol = pBlockC->getNumOfCols();
        pBlockC->resize(pBlockC->getNumOfRows(), pBlockC->getNumOfCols() + tmpBlockC.getNumOfCols());
        pBlockC->block(0, block1endCol, tmpBlockC);
    }
}

void DfLocalize_Parallel::setBlockCMatrix(const index_type MOsPerBlock, const int block1, const int block2,
                                          const TlDenseGeneralMatrix_Lapack& blockC, TlDenseGeneralMatrix_Lapack* pC) {
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(blockC.getNumOfRows() == numOfAOs);
    assert(pC->getNumOfRows() == numOfAOs);

    this->log_.info(TlUtils::format("blockC: %d,%d", blockC.getNumOfRows(), blockC.getNumOfCols()));
    {
        const index_type startMO1 = block1 * MOsPerBlock + this->startOrb_;
        const index_type endMO1 = std::min((block1 + 1) * MOsPerBlock - 1, this->endOrb_);
        this->log_.debug(
            TlUtils::format("DfLocalize_Parallel::setBlockCMatrix() 1: startMO=%d, endMO=%d", startMO1, endMO1));

        TlDenseGeneralMatrix_Lapack tmpC;
        blockC.block(0, 0, numOfAOs, endMO1 - startMO1 + 1, &tmpC);
        this->log_.debug(TlUtils::format("DfLocalize_Parallel::setBlockCMatrix() 1: tmpC(%d, %d)", tmpC.getNumOfRows(),
                                         tmpC.getNumOfCols()));
        pC->block(0, startMO1, tmpC);
    }

    if (block1 != block2) {
        // Processing when there is a late block
        const index_type startMO2 = block2 * MOsPerBlock + this->startOrb_;
        const index_type endMO2 = std::min((block2 + 1) * MOsPerBlock - 1, this->endOrb_);
        this->log_.debug(
            TlUtils::format("DfLocalize_Parallel::setBlockCMatrix() 2: startMO=%d, endMO=%d", startMO2, endMO2));

        TlDenseGeneralMatrix_Lapack tmpC;
        blockC.block(0, MOsPerBlock, numOfAOs, endMO2 - startMO2 + 1, &tmpC);
        assert(MOsPerBlock + (endMO2 - startMO2 + 1) == pC->getNumOfCols());

        pC->block(0, startMO2, tmpC);
    }
}

//
void DfLocalize_Parallel::makeJobList(const int dim) {
    const int pairs = dim * (dim + 1) / 2;
    this->jobList_.resize(pairs);

    int counter = 0;
    for (int index1 = 0; index1 < dim; ++index1) {
        // Starting 0 means that the diagonal elements are computed.
        const int max_index2 = dim - index1;
        for (int index2 = 0; index2 < max_index2; ++index2) {
            const int block_i = index2;
            const int block_j = index1 + index2;

            this->jobList_[counter] = std::make_pair(block_i, block_j);
            ++counter;
        }
    }
    assert(counter == pairs);
}

bool DfLocalize_Parallel::getJobItem(DfLocalize_Parallel::JobItem* pJob, bool isInitialized) {
    assert(pJob != NULL);

    static std::vector<JobItem>::iterator it;
    static std::size_t progressCounter = 0;
    static std::size_t progressMaxCounter = 0;
    static std::size_t progressNoticeNext = 0;

    if (isInitialized == true) {
        it = this->jobList_.begin();

        progressCounter = 0;
        progressMaxCounter = this->jobList_.size();
        progressNoticeNext = 1;

        return true;
    }

    bool answer = false;
    if (it != this->jobList_.end()) {
        *pJob = *it;
        ++it;

        // show progress
        if (progressCounter > (progressMaxCounter * 0.1) * progressNoticeNext) {
            const double ratio = (double(progressCounter) / double(progressMaxCounter)) * 100.0;
            this->log_.info(TlUtils::format("progress: %3.2f%% (%ld/%ld)", ratio, progressCounter, progressMaxCounter));
            ++progressNoticeNext;
        }

        ++progressCounter;
        answer = true;
    } else {
        this->log_.info("progress: 100% done");
    }

    return answer;
}

// lock Block
void DfLocalize_Parallel::initLockBlock(const int numOfBlocks) {
    this->lockBlocks_.resize(numOfBlocks);
    std::fill(this->lockBlocks_.begin(), this->lockBlocks_.end(), 0);
}

void DfLocalize_Parallel::lockBlock(const int block) {
    assert(block < static_cast<index_type>(this->lockBlocks_.size()));
    this->lockBlocks_[block] = 1;
}

void DfLocalize_Parallel::unlockBlock(const int block) {
    assert(block < static_cast<index_type>(this->lockBlocks_.size()));
    this->lockBlocks_[block] = 0;
}

bool DfLocalize_Parallel::isLockedBlock(const int block1, const int block2) const {
    assert(block1 < static_cast<index_type>(this->lockBlocks_.size()));
    assert(block2 < static_cast<index_type>(this->lockBlocks_.size()));
    bool answer = true;

    if ((this->lockBlocks_[block1] + this->lockBlocks_[block2]) == 0) {
        answer = false;
    }

    return answer;
}
