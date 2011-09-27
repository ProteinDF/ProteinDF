#include <iostream>

#include "DfHpq_Parallel.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlLogX.h"
#include "TlUtils.h"
#include "TlTime.h"
#include "TlFileSymmetricMatrix.h"
#include "TlFile.h"

DfHpq_Parallel::DfHpq_Parallel(TlSerializeData* pPdfParam) : DfHpq(pPdfParam)
{
    this->MS_blockSize_ = 100;
}


DfHpq_Parallel::~DfHpq_Parallel()
{
}


std::vector<DfObject::IJShellPair>
DfHpq_Parallel::getQueue(const int maxNumOfPairs, bool initialize)
{
    if (this->isMasterSlave_ == true) {
        return this->getQueue_MS(maxNumOfPairs, initialize);
    } else {
        return this->getQueue_DC(maxNumOfPairs, initialize);
    }
}


std::vector<DfObject::IJShellPair>
DfHpq_Parallel::getQueue_DC(const int maxNumOfPairs, bool initialize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const std::size_t nProc = rComm.getNumOfProc();
    const std::size_t nRank = rComm.getRank();

    static std::size_t count = 0;
    if (initialize == true) {
        count = 0;
    }

    std::vector<IJShellPair> ij_part;
    std::vector<IJShellPair> ij = DfHpq::getQueue(maxNumOfPairs, initialize);
    for (std::vector<IJShellPair>::const_iterator p = ij.begin(); p != ij.end(); ++p) {
        if ((count % nProc) == nRank) {
            ij_part.push_back(*p);
        }
        ++count;
    }

    return ij_part;
}


std::vector<DfObject::IJShellPair>
DfHpq_Parallel::getQueue_MS(const int maxNumOfPairs, bool initialize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    //const std::size_t nRank = rComm.getRank();

    enum {
        JOB_REQUEST,
        CHECK_OUT
    };

    std::vector<IJShellPair> ij_part;
    if (rComm.isMaster() == true) {
        // jobが完了したslaveの数
        int numOfCheckOutSlave = 0;
        if (initialize == true) {
            numOfCheckOutSlave = 0;
        }
        
        this->loggerTime(" start: ");
        int prevProgress = -1;
        while (numOfCheckOutSlave < (nProc -1)) {
            int msg = 0;
            int rank = 0;
            rComm.receiveDataFromAnySource(msg, &rank, TAG_HPQ_REQUEST_TASK);

            switch (msg) {
            case JOB_REQUEST: {
                ij_part = DfHpq::getQueue(this->MS_blockSize_, initialize);
                initialize = false;
                const int numOfJobs = ij_part.size();
                std::vector<index_type> jobList(numOfJobs * 2);
                for (int i = 0; i < numOfJobs; ++i) {
                    jobList[i * 2   ] = ij_part[i].nIShell;
                    jobList[i * 2 +1] = ij_part[i].nJShell;
                }
                rComm.sendData(jobList, rank, TAG_HPQ_SEND_TASK);

                int currentProgress = this->getProgress() / 10;
                if (currentProgress > prevProgress) {
                    prevProgress = currentProgress;
                    this->loggerTime(TlUtils::format("  %3d%% done.", currentProgress * 10));
                }
            }
            break;

            case CHECK_OUT:
                ++numOfCheckOutSlave;
                break;

            default:
                abort();
            }
        }

        ij_part.clear();
    } else {
        const int root = 0;
        int msg = JOB_REQUEST;
        rComm.sendData(msg, root, TAG_HPQ_REQUEST_TASK);
        std::vector<index_type> jobList;
        rComm.receiveData(jobList, root, TAG_HPQ_SEND_TASK);

        const std::size_t numOfJobs = jobList.size() / 2;
        ij_part.resize(numOfJobs);
        for (std::size_t i = 0; i < numOfJobs; ++i) {
            ij_part[i].nIShell = jobList[i * 2   ];
            ij_part[i].nJShell = jobList[i * 2 +1];
        }

        if (ij_part.empty() == true) {
            msg = CHECK_OUT;
            rComm.sendData(msg, root, TAG_HPQ_REQUEST_TASK);
        }
    }

    return ij_part;
}


void DfHpq_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfHpq::logger(str);
    }
}

void DfHpq_Parallel::parallelLogger(const std::string& str)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const std::string tmp = TlUtils::format("[%3d/%3d]> ", rComm.getRank(), rComm.getNumOfProc()) + str;

    TlLogX& log = TlLogX::getInstance();
    log << tmp;
}


// void DfHpq_Parallel::getHpq()
// {
// #ifdef HAVE_SCALAPACK
//     if (this->m_bUsingSCALAPACK == true) {
//         this->getHpq_ScaLAPACK();
//         return;
//     }
// #endif // HAVE_SCALAPACK

//     this->getHpq_LAPACK();
// }


void DfHpq_Parallel::getHpq(TlSymmetricMatrix* pHpq,
                            TlSymmetricMatrix* pHpq2)
{
    assert(pHpq != NULL);
    assert(pHpq != NULL);
    pHpq->resize(this->m_nNumOfAOs);
    pHpq2->resize(this->m_nNumOfAOs);

    if (this->isMasterSlave_ == true) {
        this->logger(" H_pq set up Master-Slave parallel computation.\n");
    } else {
        this->logger(" H_pq set up Divide & Conquer parallel computation.\n");
    }

    this->loggerTime(" make table");
    this->makeTable();

    this->loggerTime(" set aux");
    this->auxSet();

    this->loggerTime(" integral");
    this->resetCounter();
    this->getHpq_core(pHpq, pHpq2);
    if (this->m_nChargeExtrapolateNumber != 0) {
        (*pHpq2) /= static_cast<double>(this->m_nChargeExtrapolateNumber);
    }

    this->loggerTime(" finalize");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pHpq);
    rComm.allReduce_SUM(*pHpq2);
}


void DfHpq_Parallel::getHpq(TlDistributeSymmetricMatrix* pHpq,
                            TlDistributeSymmetricMatrix* pHpq2)
{
    assert(pHpq != NULL);
    assert(pHpq != NULL);
    pHpq->resize(this->m_nNumOfAOs);
    pHpq2->resize(this->m_nNumOfAOs);

    if (this->isMasterSlave_ == true) {
        this->logger(" H_pq set up Master-Slave parallel computation.\n");
    } else {
        this->logger(" H_pq set up Divide & Conquer parallel computation.\n");
    }

    TlSparseSymmetricMatrix Hpq_part(this->m_nNumOfAOs);
    TlSparseSymmetricMatrix Hpq2_part(this->m_nNumOfAOs);

    this->loggerTime(" make table");
    this->makeTable();

    this->loggerTime(" set aux");
    this->auxSet();

    this->loggerTime(" integral");
    this->resetCounter();
    this->getHpq_core(&Hpq_part, &Hpq2_part);
    if (this->m_nChargeExtrapolateNumber != 0) {
        (*pHpq2) /= static_cast<double>(this->m_nChargeExtrapolateNumber);
    }

    this->loggerTime(" finalize");
    pHpq->mergeSparseMatrix(Hpq_part);
    pHpq2->mergeSparseMatrix(Hpq2_part);
}


