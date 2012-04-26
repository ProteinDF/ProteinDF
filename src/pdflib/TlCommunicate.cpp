#include <mpi.h>
#include <iostream>
#include <cstring>
#include <cassert>

#include "TlCommunicate.h"
#include "TlVector.h"
#include "TlMatrixObject.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlFileMatrix.h"
#include "TlFileSymmetricMatrix.h"
#include "TlParameter.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"
#include "TlDistributeMatrix.h"

// minimum work memory size is 400 MB
#define DEFAULT_WORK_MEM_SIZE (400UL * 1024UL * 1024UL)

// =====================================================================
// static member parameters
//
TlCommunicate* TlCommunicate::m_pTlCommunicateInstance = NULL;
TlCommunicate::NonBlockingCommParamTableType TlCommunicate::nonBlockingCommParamTable_;

// =====================================================================
// static member functions
//
TlCommunicate& TlCommunicate::getInstance(int argc, char* argv[])
{
    if (TlCommunicate::m_pTlCommunicateInstance == NULL) {
        TlCommunicate::m_pTlCommunicateInstance = new TlCommunicate();

        // 初期化
        TlCommunicate::m_pTlCommunicateInstance->initialize(argc, argv);
    }

    return *(TlCommunicate::m_pTlCommunicateInstance);
}

TlCommunicate& TlCommunicate::getInstance()
{
    assert(TlCommunicate::m_pTlCommunicateInstance != NULL);

    return *(TlCommunicate::m_pTlCommunicateInstance);
}

// =====================================================================
// public
//
void TlCommunicate::setWorkMemSize(const std::size_t workMemSize)
{
    this->workMemSize_ = std::max<std::size_t>(workMemSize, DEFAULT_WORK_MEM_SIZE);
}


std::size_t TlCommunicate::getWorkMemSize() const
{
    return this->workMemSize_;
}


std::string TlCommunicate::getReport() const
{
    std::string answer = TlUtils::format(" #%d:\n",
                                         this->getRank());
    answer += TlUtils::format(" barrier:    %16.2f sec. %16ld times\n",
                              this->time_barrier_.getElapseTime(),
                              this->counter_barrier_);
    answer += TlUtils::format(" test:       %16.2f sec. %16ld times\n",
                              this->time_test_.getElapseTime(),
                              this->counter_test_);
    answer += TlUtils::format(" wait:       %16.2f sec. %16ld times\n",
                              this->time_wait_.getElapseTime(),
                              this->counter_wait_);
    answer += TlUtils::format(" all_reduce: %16.2f sec. %16ld times\n",
                              this->time_allreduce_.getElapseTime(),
                              this->counter_allreduce_);
    
    return answer;
}


int TlCommunicate::barrier(bool isDebugOut) const
{
    this->time_barrier_.start();
    ++(this->counter_barrier_);
    
    if (isDebugOut == true) {
        TlLogging& log = TlLogging::getInstance();
        log.debug(TlUtils::format("barrier called. times=%ld",
                                  this->counter_barrier_));
    }

    const int answer = MPI_Barrier(MPI_COMM_WORLD);

    this->time_barrier_.stop();
    return answer; 
}


double TlCommunicate::getTime()
{
    this->barrier();
    return MPI_Wtime();
}


bool TlCommunicate::checkNonBlockingTableCollision(uintptr_t key,
                                                   const NonBlockingCommParam& param,
                                                   int line) const
{
    bool answer = true;

#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        if (this->nonBlockingCommParamTable_.find(key) != this->nonBlockingCommParamTable_.end()) {
            const char isSendRecv = ((param.property & NonBlockingCommParam::SEND) != 0) ? 'S' : 'R';
            std::cerr << TlUtils::format("[%5d/%5d WARN] non-blocking table collision(%c) found in TlCommunicate: tag=%d, line=%d",
                                         this->getRank(), this->getNumOfProcs() -1,
                                         isSendRecv,
                                         param.tag,
                                         line)
                      << std::endl;
            answer = false;
        }
    }
    
    return answer;
}


bool TlCommunicate::checkNonBlockingCommunications() const
{
    bool answer = true;
    TlLogging& log = TlLogging::getInstance();

#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        if (this->nonBlockingCommParamTable_.empty() == false) {
            answer = false;
            NonBlockingCommParamTableType::const_iterator itEnd = this->nonBlockingCommParamTable_.end();
            for (NonBlockingCommParamTableType::const_iterator it = this->nonBlockingCommParamTable_.begin();
                 it != itEnd; ++it) {
                const char isSendRecv = ((it->second.property & NonBlockingCommParam::SEND) != 0) ? 'S' : 'R';
                log.warn(TlUtils::format("[%5d/%5d WARN] rest waiting communication(%c) in TlCommunicate: TAG=%d",
                                         this->getRank(), this->getNumOfProcs() -1,
                                         isSendRecv, it->second.tag));
            }
        }
    }

    return answer;
}


// =============================================================================
template<typename T>
int TlCommunicate::reduce(T* pData, const MPI_Datatype mpiType,
                          const std::size_t start, const std::size_t end,
                          const MPI_Op mpiOp, const int root)
{
    int answer = 0;
    const long length = end - start;
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    T* pBuf = new T[bufCount];

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Reduce(pData + startIndex, pBuf, bufCount,
                            mpiType, mpiOp, root, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + bufCount, pData + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const std::size_t remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Reduce(pData + startIndex, pBuf, remain,
                            mpiType, mpiOp, root, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + remain, pData + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    delete[] pBuf;
    pBuf = NULL;

    return answer;
}


int TlCommunicate::reduce_SUM(unsigned int* pData, std::size_t size, int root)
{
    return this->reduce(pData, MPI_UNSIGNED, 0, size, MPI_SUM, root);
}

int TlCommunicate::reduce_SUM(unsigned long* pData, std::size_t size, int root)
{
    return this->reduce(pData, MPI_UNSIGNED_LONG, 0, size, MPI_SUM, root);
}


// =============================================================================
template<typename T>
int TlCommunicate::allReduce(T* pData, const MPI_Datatype mpiType,
                             const std::size_t start, const std::size_t end,
                             const MPI_Op mpiOp)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);
    
    int answer = 0;
    const long length = end - start;
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    T* pBuf = new T[bufCount];

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Allreduce(pData + startIndex, pBuf, bufCount,
                               mpiType, mpiOp, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + bufCount, pData + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const std::size_t remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Allreduce(pData + startIndex, pBuf, remain,
                               mpiType, mpiOp, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + remain, pData + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    delete[] pBuf;
    pBuf = NULL;

    this->time_allreduce_.stop();
    return answer;
}


int TlCommunicate::allReduce_SUM(int& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    const int nSize = 1;

    int nReceive = 0;
    int nError = MPI_Allreduce(&rData, &nReceive, nSize,
                               MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    rData = nReceive;

    this->time_allreduce_.stop();
    return nError;
}


int TlCommunicate::allReduce_SUM(unsigned int& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    const int nSize = 1;

    int nReceive = 0;
    int nError = MPI_Allreduce(&rData, &nReceive, nSize,
                               MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    rData = nReceive;

    this->time_allreduce_.stop();
    return nError;
}


int TlCommunicate::allReduce_SUM(long& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    const int nSize = 1;

    long nReceive = 0;
    int nError = MPI_Allreduce(&rData, &nReceive, nSize,
                               MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    rData = nReceive;

    this->time_allreduce_.stop();
    return nError;
}


int TlCommunicate::allReduce_SUM(unsigned long& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    const int nSize = 1;

    unsigned long nReceive = 0;
    int nError = MPI_Allreduce(&rData, &nReceive, nSize,
                               MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    rData = nReceive;

    this->time_allreduce_.stop();
    return nError;
}



int TlCommunicate::allReduce_SUM(double& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    const int nSize = 1;

    double dReceive = 0.0;
    int nError = MPI_Allreduce(&rData, &dReceive, nSize,
                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    rData = dReceive;

    this->time_allreduce_.stop();
    return nError;
}


int TlCommunicate::allReduce_SUM(std::vector<int>& rData)
{
    return this->allReduce_SUM(rData, MPI_INT, 0, rData.size());
}


int TlCommunicate::allReduce_SUM(std::vector<unsigned int>& rData)
{
    return this->allReduce_SUM(rData, MPI_UNSIGNED, 0, rData.size());
}


int TlCommunicate::allReduce_SUM(std::vector<long>& rData)
{
    return this->allReduce_SUM(rData, MPI_LONG, 0, rData.size());
}


int TlCommunicate::allReduce_SUM(std::vector<unsigned long>& rData)
{
    return this->allReduce_SUM(rData, MPI_UNSIGNED_LONG, 0, rData.size());
}


int TlCommunicate::allReduce_SUM(std::vector<double>& rData)
{
    return this->allReduce_SUM(rData, MPI_DOUBLE, 0, rData.size());
}

template<typename T>
int TlCommunicate::allReduce_SUM(std::vector<T>& data, const MPI_Datatype mpiType,
                                 const std::size_t start, const std::size_t end)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    std::vector<T> buf(bufCount);
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Allreduce(&(data[startIndex]), &(buf[0]), bufCount,
                               mpiType, MPI_SUM, MPI_COMM_WORLD);
        std::copy(buf.begin(), buf.end(), data.begin() + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Allreduce(&(data[startIndex]), &(buf[0]), remain,
                               mpiType, MPI_SUM, MPI_COMM_WORLD);
        std::copy(buf.begin(), buf.begin() + remain, data.begin() + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    this->time_allreduce_.stop();
    return answer;
}

int TlCommunicate::allReduce_SUM(std::valarray<double>& rData)
{
    return this->allReduce_SUM(rData, 0, rData.size());
}

int TlCommunicate::allReduce_SUM(std::valarray<double>& data,
                                 const std::size_t start, const std::size_t end)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    double* pBuf = new double[bufCount];

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Allreduce((void*)&(data[startIndex]), pBuf, bufCount,
                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + bufCount, &(data[0]) + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Allreduce((void*)&(data[startIndex]), pBuf, remain,
                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::copy(pBuf, pBuf + remain, &(data[0]) + startIndex);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    delete[] pBuf;
    pBuf = NULL;

    this->time_allreduce_.stop();
    return answer;
}


int TlCommunicate::allReduce_SUM(int* pData, std::size_t length)
{
    return this->allReduce(pData, MPI_INT, 0, length, MPI_SUM);
}

int TlCommunicate::allReduce_SUM(double* pData, std::size_t length)
{
    return this->allReduce(pData, MPI_DOUBLE, 0, length, MPI_SUM);
}


// int TlCommunicate::allReduce_SUM(double* pData,
//                                  const std::size_t start, const std::size_t end)
// {
//     int answer = 0;
//     const long length = end - start;
//     const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
//     const ldiv_t tmp = std::ldiv(length, bufCount);
//     double* pBuf = new double[bufCount];

//     // 作業用メモリ分のループ
//     for (long i = 0; i < tmp.quot; ++i) {
//         std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
//         answer = MPI_Allreduce(pData + startIndex, pBuf, bufCount,
//                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         std::copy(pBuf, pBuf + bufCount, pData + startIndex);
//         if (answer != 0) {
//             std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
//                       << " anwer=" << answer
//                       << std::endl;
//         }
//     }

//     // 残り分のループ
//     if (tmp.rem != 0) {
//         const std::size_t remain = tmp.rem;
//         std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
//         answer = MPI_Allreduce(pData + startIndex, pBuf, remain,
//                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         std::copy(pBuf, pBuf + remain, pData + startIndex);
//         if (answer != 0) {
//             std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
//                       << " anwer=" << answer
//                       << std::endl;
//         }
//     }

//     delete[] pBuf;
//     pBuf = NULL;

//     return answer;
// }


int TlCommunicate::allReduce_SUM(TlVector& rVector)
{
    return this->allReduce_SUM(rVector.data_, rVector.getSize()); // this class is friend class of TlVector.
}


int TlCommunicate::allReduce_SUM(TlMatrix& rMatrix)
{
    return this->allReduce_SUM(rMatrix.data_, rMatrix.getNumOfElements());
}


int TlCommunicate::allReduce_SUM(TlSymmetricMatrix& rMatrix)
{
    return this->allReduce_SUM(rMatrix.data_, rMatrix.getNumOfElements());
}


// int TlCommunicate::allReduce_SUM(TlMmapMatrix& rMatrix)
// {
//     const std::size_t numOfElements = rMatrix.getNumOfRows() * rMatrix.getNumOfCols();

//     return this->allReduce_SUM(rMatrix.data_, numOfElements);
// }


// int TlCommunicate::allReduce_SUM(TlMmapSymmetricMatrix& rMatrix)
// {
//     // rMatrix.getNumOfRows() == rMatrix.getNumOfCols()
//     const std::size_t dim = rMatrix.getNumOfRows();
//     const std::size_t numOfElements = dim * (dim +1) / 2;

//     return this->allReduce_SUM(rMatrix.data_, numOfElements);
// }


int TlCommunicate::allReduce_AND(bool& rData)
{
    int nData = (rData == true) ? 1 : 0;
    int nError = this->allReduce_SUM(nData);

    rData = (nData != 0) ?  true : false;

    return nError;
}

int TlCommunicate::allReduce_MAX(int& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    int receive = 0;
    const int count = 1;
    int nError = MPI_Allreduce(&rData, &receive, count,
                               MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    rData = receive;

    this->time_allreduce_.stop();
    return nError;
}

int TlCommunicate::allReduce_MIN(int& rData)
{
    this->time_allreduce_.start();
    ++(this->counter_allreduce_);

    int receive = 0;
    const int count = 1;
    int nError = MPI_Allreduce(&rData, &receive, count,
                               MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    rData = receive;

    this->time_allreduce_.stop();
    return nError;
}

int TlCommunicate::gatherToMaster(TlSparseMatrix& rMatrix)
{
    const int nProc = this->getNumOfProcs();
    const int nRank = this->getRank();

    int count = 1;
    int nDiv = 1;
    while (nDiv < nProc) {
        nDiv *= 2;

        const div_t d = std::div(nRank, nDiv);
        const int nReceiver = d.quot * nDiv;
        const int nSender = nReceiver +(nDiv /2);

        if (nSender < nProc) {
            const int nTag = nSender;
            if (nRank == nSender) {
                // my rank is sender
                std::vector<unsigned long> indexData;
                std::vector<double> valData;
                TlSparseMatrix::const_iterator pEnd = rMatrix.end();
                for (TlSparseMatrix::const_iterator p = rMatrix.begin(); p != pEnd; ++p) {
                    //int nRow = p->first.row;
                    //int nCol = p->first.col;
                    unsigned long index = p->first;
                    double dValue = p->second;
                    indexData.push_back(index);
                    valData.push_back(dValue);
                }
                this->sendData(indexData, nReceiver, nTag);
                this->sendData(valData, nReceiver, nTag);

            } else if (nRank == nReceiver) {
                // my rank is receiver
                std::vector<unsigned long> indexData;
                std::vector<double> valData;
                this->receiveData(indexData, nSender, nTag);
                this->receiveData(valData, nSender, nTag);
                assert(indexData.size() == valData.size());
                long nSize = indexData.size();
                for (long i = 0; i < nSize; ++i) {
                    rMatrix.m_aMatrix[indexData[i]] += valData[i];
                    //rMatrix(rowData[i], colData[i]) += valData[i];
                }
            }
        }

        count++;
        this->barrier();
    }

    return 0;
}


// 1:1
int TlCommunicate::sendData(bool data, int destination, int tag)
{
    const int yn = (data == false) ? 0 : 1;
    return this->sendData(yn, destination, tag);
}


int TlCommunicate::sendData(int data, int nDestination, int nTag)
{
    int nErr = MPI_Send((void*)&data, 1, MPI_INT, nDestination, nTag, MPI_COMM_WORLD);
    return nErr;
}


int TlCommunicate::sendData(unsigned int data, int nDestination, int nTag)
{
//     if (isDebugOut == true) {
//         std::cerr << TlUtils::format("[%d] sendData(unsigned int) to [%d] data=%u, using tag=[%d]",
//                                      this->getRank(), nDestination, data, nTag) << std::endl;
//     }
    int nErr = MPI_Send((void*)&data, 1, MPI_UNSIGNED, nDestination, nTag, MPI_COMM_WORLD);
    return nErr;
}

int TlCommunicate::sendData(long data, int nDestination, int nTag)
{
    int nErr = MPI_Send((void*)&data, 1, MPI_LONG, nDestination, nTag, MPI_COMM_WORLD);
    return nErr;
}


int TlCommunicate::sendData(unsigned long data, int nDestination, int nTag)
{
    int nErr = MPI_Send((void*)&data, 1, MPI_UNSIGNED_LONG, nDestination, nTag, MPI_COMM_WORLD);
    return nErr;
}

int TlCommunicate::sendData(double data, int nDestination, int nTag)
{
    int nErr = MPI_Send((void*)&data, 1, MPI_DOUBLE, nDestination, nTag, MPI_COMM_WORLD);
    return nErr;
}


template<typename T>
int TlCommunicate::sendData(const std::vector<T>& data, const MPI_Datatype mpiType,
                            const std::size_t start, const std::size_t end,
                            const int destination, const int tag)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Send((void*)&(data[startIndex]), bufCount,
                          mpiType, destination, tag, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Send((void*)&(data[startIndex]), remain,
                          mpiType, destination, tag, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}


int TlCommunicate::sendData(const std::vector<int>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, MPI_INT, 0, size, destination, tag);
    }
    return answer;
}

int TlCommunicate::sendData(const std::vector<unsigned int>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, MPI_UNSIGNED, 0, size, destination, tag);
    }
    return answer;
}

int TlCommunicate::sendData(const std::vector<long>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, MPI_LONG, 0, size, destination, tag);
    }
    return answer;
}

int TlCommunicate::sendData(const std::vector<unsigned long>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, MPI_UNSIGNED_LONG, 0, size, destination, tag);
    }
    return answer;
}

int TlCommunicate::sendData(const std::vector<double>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, MPI_DOUBLE, 0, size, destination, tag);
    }
    return answer;
}

int TlCommunicate::sendData(const std::valarray<double>& data, int destination, int tag)
{
    const size_t size = data.size();
    int answer = this->sendData(size, destination, tag);

    if ((answer == 0) && (size > 0)) {
        answer = this->sendData(data, 0, size, destination, tag);
    }

    return answer;
}

int TlCommunicate::sendData(const std::valarray<double>& data,
                            const std::size_t start, const std::size_t end,
                            int destination, int tag)
{
    // this const_cast for PGI compiler compile error:
    // "error: expression must be an lvalue or a function designator"
    std::valarray<double>& data_tmp = const_cast<std::valarray<double>& >(data);

    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Send((void*)&(data_tmp[startIndex]), bufCount,
                          MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Send((void*)&(data_tmp[startIndex]), remain,
                          MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}


int TlCommunicate::sendData(const std::string& data, int nDestination, int nTag)
{
    const int nSize = static_cast<int>(data.length());
    char* pBuff = new char[nSize];
    std::strncpy(pBuff, data.c_str(), nSize);

    this->sendData(nSize, nDestination, nTag);
    int nErr = MPI_Send((void*)pBuff, nSize, MPI_CHAR, nDestination, nTag, MPI_COMM_WORLD);

    delete[] pBuff;
    pBuff = NULL;

    return nErr;
}


int TlCommunicate::sendData(const TlVector& data, int destination, int nTag)
{
    const std::size_t dim = data.getSize();
    int nErr = this->sendData(dim, destination, nTag);
    if (nErr == 0) {
        nErr = this->sendDataX(data.data_, MPI_DOUBLE, 0, dim, destination, nTag);
    }

    return nErr;
}

int TlCommunicate::sendData(const TlMatrix& data, const int destination, const int tag)
{
    const int headerSize = 2;
    TlMatrixObject::index_type* pHeader = new TlMatrixObject::index_type[headerSize];
    pHeader[0] = data.getNumOfRows();
    pHeader[1] = data.getNumOfCols();
    int err = this->sendDataX(pHeader, headerSize, destination, tag);
    this->log_.debug(TlUtils::format("TlCommunicate::sendData(TlMatrix): dest=%d, tag=%d, row=%d, col=%d, err=%d.",
                                     destination, tag, data.getNumOfRows(), data.getNumOfCols(), err));
    delete[] pHeader;
    pHeader = NULL;
    if (err != 0) {
        return err;
    }

    err = this->sendDataX(data.data_, data.getNumOfElements(), destination, tag);
    this->log_.debug(TlUtils::format("TlCommunicate::sendData(TlMatrix): dest=%d, tag=%d, size=%ld, err=%d.",
                                     destination, tag, data.getNumOfElements(), err));
    return err;
}

int TlCommunicate::sendData(const TlSymmetricMatrix& data, int nDestination, int nTag)
{
    const std::size_t dim = data.getNumOfRows();
    int nErr = this->sendData(dim, nDestination, nTag);
    if (nErr == 0) {
        nErr = this->sendDataX(data.data_, MPI_DOUBLE, 0, data.getNumOfElements(), nDestination, nTag);
    }

    return nErr;
}

int TlCommunicate::sendData(const TlSparseSymmetricMatrix& data, int nDestination, int nTag)
{
    std::vector<std::size_t> buf(2);
    buf[0] = data.getNumOfRows(); // == getNumOfCols()
    buf[1] = data.getSize();
    int nErr = this->sendData(buf, nDestination, nTag);

    if (nErr == 0) {
        const std::size_t size = data.getSize();
        std::vector<unsigned long> index(size);
        std::vector<double> value(size);
        int count = 0;
        for (TlSparseSymmetricMatrix::const_iterator p = data.begin(); p != data.end(); ++p) {
            index[count] = p->first;
            value[count] = p->second;
            ++count;
        }
        assert(std::size_t(count) == size);

        nErr = this->sendData(index, nDestination, nTag);
        nErr = this->sendData(value, nDestination, nTag);
    }

    return nErr;
}

int TlCommunicate::sendData(const TlPartialSymmetricMatrix& data, int nDestination, int nTag)
{
    std::vector<std::size_t> buf(6);
    buf[0] = data.getNumOfRows();
    buf[1] = data.getNumOfCols();
    buf[2] = data.getStartRow();
    buf[3] = data.getStartCol();
    buf[4] = data.getRowRange();
    buf[5] = data.getColRange();
    int nErr = this->sendData(buf, nDestination, nTag);

    if (nErr == 0) {
        const std::size_t bufCount = data.getRowRange() * data.getColRange();
        nErr = this->sendDataX(data.pData_, MPI_DOUBLE, 0, bufCount, nDestination, nTag);
    }

    return nErr;
}

// sendDataX ===================================================================
int TlCommunicate::sendDataX(const int* pData, const std::size_t size,
                             const int dest, const int tag)
{
    return this->sendDataX(pData, MPI_INT, 0, size, dest, tag);
}

int TlCommunicate::sendDataX(const unsigned int* pData, const std::size_t size,
                             const int dest, const int tag)
{
    return this->sendDataX(pData, MPI_UNSIGNED, 0, size, dest, tag);
}

int TlCommunicate::sendDataX(const double* pData, const std::size_t size,
                             const int dest, const int tag)
{
    return this->sendDataX(pData, MPI_DOUBLE, 0, size, dest, tag);
}

template<typename T>
int TlCommunicate::sendDataX(const T* pData, const MPI_Datatype mpiType,
                             const std::size_t start, const std::size_t end,
                             const int dest, const int tag)
{
    const int answer = MPI_Send((void*)(pData + start), (end - start),
                                mpiType, dest, tag, MPI_COMM_WORLD);
    if (answer != 0) {
        std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                  << " anwer=" << answer
                  << std::endl;
    }

    return answer;
}

// receiveData =================================================================
int TlCommunicate::receiveData(bool& data, int src, int tag)
{
    int yn = 0;
    int err = this->receiveData(yn, src, tag);
    data = (yn != 0);

    return err;
}


int TlCommunicate::receiveData(int& rData, int nSrc, int nTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&rData, 1, MPI_INT, nSrc, nTag, MPI_COMM_WORLD, &status);

    return nErr;
}


int TlCommunicate::receiveData(unsigned int& rData, int nSrc, int nTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&rData, 1, MPI_UNSIGNED, nSrc, nTag, MPI_COMM_WORLD, &status);

    return nErr;
}

int TlCommunicate::receiveData(long& rData, int nSrc, int nTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&rData, 1, MPI_LONG, nSrc, nTag, MPI_COMM_WORLD, &status);

    return nErr;
}

int TlCommunicate::receiveData(unsigned long& rData, int nSrc, int nTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&rData, 1, MPI_UNSIGNED_LONG, nSrc, nTag, MPI_COMM_WORLD, &status);

    return nErr;
}

int TlCommunicate::receiveData(double& rData, int nSrc, int nTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&rData, 1, MPI_DOUBLE, nSrc, nTag, MPI_COMM_WORLD, &status);

    return nErr;
}

template<typename T>
int TlCommunicate::receiveData(std::vector<T>& data, const MPI_Datatype mpiType,
                               const std::size_t start, const std::size_t end,
                               const int src, const int tag)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        MPI_Status status;
        answer = MPI_Recv(&(data[startIndex]), bufCount, mpiType,
                          src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        MPI_Status status;
        answer = MPI_Recv(&(data[startIndex]), remain, mpiType,
                          src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer = 0;
}


// int TlCommunicate::receiveData(double* pData,
//                                const std::size_t start, const std::size_t end,
//                                const int src, const int tag)
// {
//     int answer = 0;
//     const long length = static_cast<long>(end - start);
//     const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
//     const ldiv_t tmp = std::ldiv(length, bufCount);
//     MPI_Status status;

//     // 作業用メモリ分のループ
//     for (long i = 0; i < tmp.quot; ++i) {
//         const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
//         answer = MPI_Recv((void*)(pData + startIndex), bufCount, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
//         if (answer != 0) {
//             std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
//                       << " anwer=" << answer
//                       << std::endl;
//         }
//     }

//     // 残り分のループ
//     if (tmp.rem != 0) {
//         const int remain = tmp.rem;
//         const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
//         answer = MPI_Recv((void*)(pData + startIndex), remain, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
//         if (answer != 0) {
//             std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
//                       << " anwer=" << answer
//                       << std::endl;
//         }
//     }

//     return answer;
// }


int TlCommunicate::receiveData(std::vector<int>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, MPI_INT, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::vector<unsigned int>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, MPI_UNSIGNED, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::vector<long>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, MPI_LONG, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::vector<unsigned long>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, MPI_UNSIGNED_LONG, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::vector<double>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, MPI_DOUBLE, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::valarray<double>& data,
                               const std::size_t start, const std::size_t end,
                               const int src, const int tag)
{
    int answer = 0;
    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    MPI_Status status;

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Recv((void*)&(data[startIndex]), bufCount,
                          MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Recv((void*)&(data[startIndex]), remain,
                          MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}

int TlCommunicate::receiveData(std::valarray<double>& data, int src, int tag)
{
    std::size_t size = 0;
    this->receiveData(size, src, tag);
    data.resize(size);

    return this->receiveData(data, 0, size, src, tag);
}

int TlCommunicate::receiveData(std::string& rData, int nSrc, int nTag)
{
    int nErr = 0;

    int nSize = 0;
    this->receiveData(nSize, nSrc, nTag);

    char* pBuff = new char[nSize];
    MPI_Status status;
    nErr = MPI_Recv((void*)pBuff, nSize, MPI_CHAR, nSrc, nTag, MPI_COMM_WORLD, &status);
    rData = std::string(pBuff, nSize);

    delete[] pBuff;
    pBuff = NULL;

    return nErr;
}


int TlCommunicate::receiveData(TlVector& rData, int src, int tag)
{
    std::size_t dim = 0;
    int nErr = this->receiveData(dim, src, tag);
    if (nErr == 0) {
        rData.resize(dim);
        nErr = this->receiveDataX(rData.data_, MPI_DOUBLE, 0, dim, src, tag);
    }

    return nErr;
}

int TlCommunicate::receiveData(TlMatrix& data, const int src, const int tag)
{
    const int headerSize = 2;
    TlMatrixObject::index_type* pHeader = new TlMatrixObject::index_type[headerSize];
    int err = this->receiveDataX(pHeader, headerSize, src, tag);
    const TlMatrixObject::index_type numOfRows = pHeader[0];
    const TlMatrixObject::index_type numOfCols = pHeader[1];
    this->log_.debug(TlUtils::format("TlCommunicate::recvData(TlMatrix): src=%d, tag=%d, row=%d, col=%d, err=%d.",
                                     src, tag, numOfRows, numOfCols, err));
    
    delete[] pHeader;
    pHeader = NULL;
    if (err != 0) {
        return err;
    }
    
    data.resize(numOfRows, numOfCols);
    err = this->receiveDataX(data.data_, MPI_DOUBLE, 0, data.getNumOfElements(), src, tag);
    this->log_.debug(TlUtils::format("TlCommunicate::recvData(TlMatrix): src=%d, tag=%d, size=%ld, err=%d.",
                                     src, tag, data.getNumOfElements(), err));
    
    return err;
}

int TlCommunicate::receiveData(TlSymmetricMatrix& rData, int nSrc, int nTag)
{
    std::size_t dim = 0;
    int nErr = this->receiveData(dim, nSrc, nTag);
    if (nErr == 0) {
        rData.resize(dim);
        nErr = this->receiveDataX(rData.data_, MPI_DOUBLE, 0, rData.getNumOfElements(), nSrc, nTag);
    }

    return nErr;
}


int TlCommunicate::receiveData(TlSparseSymmetricMatrix& rData, int nSrc, int nTag)
{
    std::vector<std::size_t> header;
    int nErr = this->receiveData(header, nSrc, nTag);

    if (nErr == 0) {
        assert(header.size() == 2);

        const std::size_t numOfSize = header[0];
        const std::size_t size = header[1];

        std::vector<unsigned long> index;
        std::vector<double> value;
        nErr = this->receiveData(index, nSrc, nTag);
        assert(index.size() == size);
        if (nErr == 0) {
            nErr = this->receiveData(value, nSrc, nTag);
            assert(value.size() == size);
            rData.clear();
            rData.resize(numOfSize);
            for (std::size_t i = 0; i < size; ++i) {
                rData.set(std::pair<unsigned long, double>(index[i], value[i]));
            }
        }
    }

    return nErr;
}

int TlCommunicate::receiveData(TlPartialSymmetricMatrix& rData, int nSrc, int nTag)
{
    std::vector<std::size_t> header;
    int nErr = this->receiveData(header, nSrc, nTag);

    if (nErr == 0) {
        assert(header.size() == 6);
        const std::size_t globalNumOfRows = header[0];
        //const std::size_t globalNumOfCols = header[1];
        const std::size_t startRow = header[2];
        const std::size_t startCol = header[3];
        const std::size_t rowRange = header[4];
        const std::size_t colRange = header[5];

        TlPartialSymmetricMatrix tmp(globalNumOfRows, startRow, startCol, rowRange, colRange);
        const std::size_t bufCount = tmp.getRowRange() * tmp.getColRange();
        nErr = this->receiveDataX(tmp.pData_, MPI_DOUBLE, 0, bufCount, nSrc, nTag);
        rData = tmp;
    }

    return nErr;
}

// =============================================================================
int TlCommunicate::receiveDataX(int* pData, const std::size_t size,
                                const int src, const int tag)
{
    return this->receiveDataX(pData, MPI_INT,
                              0, size, src, tag);
}

int TlCommunicate::receiveDataX(unsigned int* pData, const std::size_t size,
                                const int src, const int tag)
{
    return this->receiveDataX(pData, MPI_UNSIGNED,
                              0, size, src, tag);
}

int TlCommunicate::receiveDataX(double* pData, const std::size_t size,
                                const int src, const int tag)
{
    return this->receiveDataX(pData, MPI_DOUBLE,
                              0, size, src, tag);
}

template<typename T>
int TlCommunicate::receiveDataX(T* pData, const MPI_Datatype mpiType,
                               const std::size_t start, const std::size_t end,
                               const int src, const int tag)
{
    int answer = 0;
    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    MPI_Status status;

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Recv((void*)(pData + startIndex), bufCount,
                          mpiType, src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Recv((void*)(pData + startIndex), remain,
                          mpiType, src, tag, MPI_COMM_WORLD, &status);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}


// =============================================================================
template<typename T>
int TlCommunicate::receiveDataFromAnySource(T& data, const MPI_Datatype mpiType,
                                            int* pSrc, const int tag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&data, 1, mpiType, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

    assert(pSrc != NULL);
    *pSrc = status.MPI_SOURCE;

    return nErr;
}

template<typename T>
int TlCommunicate::receiveDataFromAnySource(T& data, const MPI_Datatype mpiType,
                                            int* pSrc, int* pTag)
{
    MPI_Status status;
    int nErr = MPI_Recv((void*)&data, 1, mpiType, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    assert(pSrc != NULL);
    *pSrc = status.MPI_SOURCE;
    if (pTag != NULL) {
        *pTag = status.MPI_TAG;
    }

    return nErr;
}


int TlCommunicate::receiveDataFromAnySource(int& data, int* pSrc, const int tag)
{
    return this->receiveDataFromAnySource(data, MPI_INT, pSrc, tag);
}

int TlCommunicate::receiveDataFromAnySource(int& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_INT, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(unsigned int& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_UNSIGNED, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(long& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_LONG, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(unsigned long& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_UNSIGNED_LONG, pSrc, pTag);
}

template<typename T>
int TlCommunicate::receiveDataFromAnySource(std::vector<T>& data, const MPI_Datatype mpiType,
                                            int* pSrc, int*pTag)
{
    std::size_t size = 0;
    int src = 0;
    int tag = 0;
    this->receiveDataFromAnySource(size, &src, &tag);
    data.resize(size);

    const int answer = this->receiveData(data, mpiType, 0, size, src, tag);

    assert(pSrc != NULL);
    *pSrc = src;
    if (pTag != NULL) {
        *pTag = tag;
    }

    return answer;
}

int TlCommunicate::receiveDataFromAnySource(std::vector<int>& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_INT, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(std::vector<unsigned int>& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_UNSIGNED, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(std::vector<long>& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_LONG, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(std::vector<unsigned long>& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_UNSIGNED_LONG, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(std::vector<double>& data, int* pSrc, int* pTag)
{
    return this->receiveDataFromAnySource(data, MPI_DOUBLE, pSrc, pTag);
}

int TlCommunicate::receiveDataFromAnySource(TlMatrix& data, int* pSrc, int tag)
{
    assert(pSrc != NULL);

    const int headerSize = 2;
    TlMatrixObject::index_type* pHeader = new TlMatrixObject::index_type[headerSize];
    int src = 0;
    int err = this->receiveDataFromAnySourceX(pHeader, headerSize, &src, tag);
    if (err != 0) {
        return err;
    }
    
    *pSrc = src;

    const TlMatrixObject::index_type numOfRows = pHeader[0];
    const TlMatrixObject::index_type numOfCols = pHeader[1];
    delete[] pHeader;
    pHeader = NULL;
    
    data.resize(numOfRows, numOfCols);
    return this->receiveDataX(data.data_, MPI_DOUBLE, 0, data.getNumOfElements(), src, tag);
}

int TlCommunicate::receiveDataFromAnySource(TlSparseSymmetricMatrix& rData, int* pSrc, int* pTag)
{
    std::vector<std::size_t> header;
    int nSrc = 0;
    int nTag = 0;
    int nErr = this->receiveDataFromAnySource(header, &nSrc, &nTag);

    assert(pSrc != NULL);
    *pSrc = nSrc;
    if (pTag != NULL) {
        *pTag = nTag;
    }

    if (nErr == 0) {
        assert(header.size() == 2);

        const std::size_t numOfSize = header[0];
        const std::size_t size = header[1];

        std::vector<unsigned long> index;
        std::vector<double> value;
        nErr = this->receiveData(index, nSrc, nTag);
        assert(index.size() == size);
        if (nErr == 0) {
            nErr = this->receiveData(value, nSrc, nTag);
            assert(value.size() == size);
            rData.clear();
            rData.resize(numOfSize);
            for (std::size_t i = 0; i < size; ++i) {
                rData.set(std::pair<unsigned long, double>(index[i], value[i]));
            }
        }
    }

    return nErr;
}

int TlCommunicate::receiveDataFromAnySource(TlPartialSymmetricMatrix& rData, int* pSrc, int* pTag)
{
    std::vector<std::size_t> header;
    int nSrc = 0;
    int nTag = 0;
    int nErr = this->receiveDataFromAnySource(header, &nSrc, &nTag);

    if (nErr == 0) {
        assert(header.size() == 6);
        const std::size_t globalNumOfRows = header[0];
        //const std::size_t globalNumOfCols = header[1];
        const std::size_t startRow = header[2];
        const std::size_t startCol = header[3];
        const std::size_t rowRange = header[4];
        const std::size_t colRange = header[5];

        TlPartialSymmetricMatrix tmp(globalNumOfRows, startRow, startCol, rowRange, colRange);
        const std::size_t bufCount = tmp.getRowRange() * tmp.getColRange();
        nErr = this->receiveDataX(tmp.pData_, MPI_DOUBLE, 0, bufCount, nSrc, nTag);
        rData = tmp;

        assert(pSrc != NULL);
        *pSrc = nSrc;
        if (pTag != NULL) {
            *pTag = nTag;
        }
    }

    return nErr;
}

// =============================================================================

template<typename T>
int TlCommunicate::receiveDataFromAnySourceX(T* pData, const MPI_Datatype mpiType,
                                             const std::size_t start, const std::size_t end,
                                             int* pSrc, const int tag)
{
    int answer = 0;
    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T));
    const ldiv_t tmp = std::ldiv(length, bufCount);
    MPI_Status status;
    int src = 0;
    bool isSrcDefined = false;

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        if (isSrcDefined != true) {
            answer = MPI_Recv((void*)(pData + startIndex), bufCount,
                              mpiType, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            src = status.MPI_SOURCE;
            isSrcDefined = true;
        } else {
            answer = MPI_Recv((void*)(pData + startIndex), bufCount,
                              mpiType, src, tag, MPI_COMM_WORLD, &status);
        }
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        if (isSrcDefined != true) {
            answer = MPI_Recv((void*)(pData + startIndex), remain,
                              mpiType, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            src = status.MPI_SOURCE;
            isSrcDefined = true;
        } else {
            answer = MPI_Recv((void*)(pData + startIndex), remain,
                              mpiType, src, tag, MPI_COMM_WORLD, &status);
        }
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    if (pSrc != NULL) {
        *pSrc = src;
    }
    
    return answer;
}


int TlCommunicate::receiveDataFromAnySourceX(int* pData, const std::size_t size,
                                             int* pSrc, const int tag)
{
    return this->receiveDataFromAnySourceX(pData, MPI_INT, 0, size, pSrc, tag);
}

// =============================================================================
template<typename T>
int TlCommunicate::iSendData(const T& data, const MPI_Datatype mpiType,
                             const int destination, const int tag)
{
    MPI_Request* pRequest = new MPI_Request;
    const int answer = MPI_Isend((void*)&data, 1, mpiType, destination, tag, MPI_COMM_WORLD, pRequest);
    if (answer == 0) {
        std::vector<uintptr_t> requests;
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        const uintptr_t key = reinterpret_cast<uintptr_t>((void*)&data);
        const NonBlockingCommParam param(requests, tag,
                                         NonBlockingCommParam::SEND);
        this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
        {
            this->nonBlockingCommParamTable_[key] = param;
        }
    }

    // don't call "delete pRequest"
    // because the object which is pointed will be deleted on calling wait().
    
    return answer;
}


int TlCommunicate::iSendData(const int& data, const int destination, const int tag)
{
    return this->iSendData(data, MPI_INT, destination, tag);
}

int TlCommunicate::iSendData(const unsigned int& data, const int destination, const int tag)
{
    return this->iSendData(data, MPI_UNSIGNED, destination, tag);
}

int TlCommunicate::iSendData(const long& data, const int destination, const int tag)
{
    return this->iSendData(data, MPI_LONG, destination, tag);
}

int TlCommunicate::iSendData(const unsigned long& data, const int destination, const int tag)
{
    return this->iSendData(data, MPI_UNSIGNED_LONG, destination, tag);
}

template<typename T>
int TlCommunicate::iSendData(const std::vector<T>& data, const MPI_Datatype mpiType,
                             const std::size_t start, const std::size_t end,
                             const int destination, const int tag)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    std::vector<T> buf(bufCount);
    const ldiv_t tmp = std::ldiv(length, bufCount);

    std::vector<uintptr_t> requests;

    // 作業用メモリ分のループ
    if (tmp.quot != 0) {
        for (long i = 0; i < tmp.quot; ++i) {
            std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
            MPI_Request* pRequest = new MPI_Request;
            answer = MPI_Isend((void*)&(data[startIndex]), bufCount,
                               mpiType, destination, tag, MPI_COMM_WORLD, pRequest);
            requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
            if (answer != 0) {
                std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                          << " anwer=" << answer
                          << std::endl;
            }
            // don't delete pRequest!
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        MPI_Request* pRequest = new MPI_Request;
        answer = MPI_Isend((void*)&(data[startIndex]), remain,
                           mpiType, destination, tag, MPI_COMM_WORLD, pRequest);
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
        // don't delete pRequest!
    }

    const uintptr_t key = reinterpret_cast<uintptr_t>((void*)&data[start]);
    const NonBlockingCommParam param(requests, tag,
                                     NonBlockingCommParam::SEND);
    this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        this->nonBlockingCommParamTable_[key] = param;
    }

    return answer;
}

int TlCommunicate::iSendData(const std::vector<int>& data, const int destination, const int tag)
{
    const std::size_t size = data.size();
    int answer = this->iSendData(size, destination, tag);
    this->wait(size);

    answer |= this->iSendData(data, MPI_INT, 0, size, destination, tag);
    return answer;
}

int TlCommunicate::iSendData(const std::vector<double>& data, int destination, int tag)
{
    const std::size_t size = data.size();
    int answer = this->iSendData(size, destination, tag);
    this->wait(size);

    answer |= this->iSendData(data, MPI_DOUBLE, 0, size, destination, tag);
    return answer;
}


// int TlCommunicate::iSendData(const TlVector& data, int dest, int tag)
// {
//     std::vector<uintptr_t> requests;

//     MPI_Request* pRequest1 = new MPI_Request;
//     int answer = MPI_Isend(&(data.size_), 1,
//                            MPI_UNSIGNED_LONG, dest, tag,
//                            MPI_COMM_WORLD, pRequest1);
//     requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest1));

//     if (answer == 0) {
//         MPI_Request* pRequest2 = new MPI_Request;
//         answer = MPI_Isend(&(data.data_), data.size_,
//                            MPI_DOUBLE, dest, tag,
//                            MPI_COMM_WORLD, pRequest2);
//         requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest2));
        
//         uintptr_t key = reinterpret_cast<uintptr_t>((void*)&data);
//         this->nonBlockingCommParamTable_[key] = NonBlockingCommParam(NonBlockingCommParam::SEND,
//                                                                      requests);
//     }

//     return answer;
// }


// =============================================================================
template<typename T>
int TlCommunicate::iSendDataX(const T* pData, const MPI_Datatype mpiType,
                              const std::size_t start, const std::size_t end,
                              const int dest, const int tag)
{
    MPI_Request* pRequest = new MPI_Request;
    const int answer = MPI_Isend((void*)(pData + start), (end - start),
                                 mpiType, dest, tag,
                                 MPI_COMM_WORLD, pRequest);
    
    if (answer == 0) {
        std::vector<uintptr_t> requests;
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        const uintptr_t key = reinterpret_cast<uintptr_t>((void*)pData);
        const NonBlockingCommParam param(requests, tag,
                                         NonBlockingCommParam::SEND);
        this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
        {
            this->nonBlockingCommParamTable_[key] = param;
        }
    }

    return answer;
}


int TlCommunicate::iSendDataX(const int* pData, const std::size_t size,
                              const int dest, const int tag)
{
    return this->iSendDataX(pData, MPI_INT, 0, size, dest, tag);
}


int TlCommunicate::iSendDataX(const double* pData, const std::size_t size,
                              const int dest, const int tag)
{
    return this->iSendDataX(pData, MPI_DOUBLE, 0, size, dest, tag);
}


// =============================================================================
template<typename T>
int TlCommunicate::iReceiveData(T& data, const MPI_Datatype mpiType,
                                const int src, const int tag)
{
    MPI_Request* pRequest = new MPI_Request;
    int err = MPI_Irecv((void*)&data, 1, mpiType, src, tag, MPI_COMM_WORLD, pRequest);
    if (err == 0) {
        std::vector<uintptr_t> requests;
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        const uintptr_t key = reinterpret_cast<uintptr_t>((void*)&data);
        const NonBlockingCommParam param(requests, tag,
                                         NonBlockingCommParam::RECV);
        this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
        {
            this->nonBlockingCommParamTable_[key] = param;
        }
    }
    // don't delete pRequest

    return err;
}


int TlCommunicate::iReceiveData(int& data, const int src, const int tag)
{
    return this->iReceiveData(data, MPI_INT, src, tag);
}

int TlCommunicate::iReceiveData(unsigned int& data, const int src, const int tag)
{
    return this->iReceiveData(data, MPI_UNSIGNED, src, tag);
}

int TlCommunicate::iReceiveData(long& data, const int src, const int tag)
{
    return this->iReceiveData(data, MPI_LONG, src, tag);
}

int TlCommunicate::iReceiveData(unsigned long& data, const int src, const int tag)
{
    return this->iReceiveData(data, MPI_UNSIGNED_LONG, src, tag);
}

// iReceiveData (array) ========================================================
template<typename T>
int TlCommunicate::iReceiveData(std::vector<T>& data, const MPI_Datatype mpiType,
                                const std::size_t start, const std::size_t end,
                                const int src, const int tag)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    std::vector<T> buf(bufCount);
    const ldiv_t tmp = std::ldiv(length, bufCount);

    std::vector<uintptr_t> requests;

    // 作業用メモリ分のループ
    if (tmp.quot != 0) {
        for (long i = 0; i < tmp.quot; ++i) {
            std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
            MPI_Request* pRequest = new MPI_Request;
            answer = MPI_Irecv(&(data[startIndex]), bufCount, mpiType,
                               src, tag, MPI_COMM_WORLD, pRequest);
            requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
            if (answer != 0) {
                std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                          << " anwer=" << answer
                          << std::endl;
            }
            // don't delete pRequest
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        MPI_Request* pRequest = new MPI_Request;
        answer = MPI_Irecv(&(data[startIndex]), remain, mpiType,
                           src, tag, MPI_COMM_WORLD, pRequest);
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
        // don't delete pRequest
    }

    const uintptr_t key = reinterpret_cast<uintptr_t>((void*)&(data[start]));
    const NonBlockingCommParam param(requests, tag,
                                     NonBlockingCommParam::RECV);
    this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        this->nonBlockingCommParamTable_[key] = param;
    }

    return answer;
}


int TlCommunicate::iReceiveData(std::vector<int>& data, const int src, const int tag)
{
    std::size_t size = 0;
    this->iReceiveData(size, src, tag);
    this->wait(size);
    data.resize(size);

    return this->iReceiveData(data, MPI_INT, 0, size, src, tag);
}


int TlCommunicate::iReceiveData(std::vector<double>& data, int src, int tag)
{
    std::size_t size = 0;
    this->iReceiveData(size, src, tag);
    this->wait(size);
    data.resize(size);

    return this->iReceiveData(data, MPI_DOUBLE, 0, size, src, tag);
}


// int TlCommunicate::iReceiveData(TlVector& data, int src, int tag)
// {
//     std::vector<uintptr_t> requests;

//     MPI_Request* pRequest1 = new MPI_Request;
//     int answer = MPI_Isend(&(data.size_), 1,
//                            MPI_UNSIGNED_LONG, dest, tag,
//                            MPI_COMM_WORLD, pRequest1);
//     requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest1));

//     if (answer == 0) {
//         MPI_Request* pRequest2 = new MPI_Request;
//         answer = MPI_Isend(&(data.data_), data.size_,
//                            MPI_DOUBLE, dest, tag,
//                            MPI_COMM_WORLD, pRequest2);
//         requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest2));
        
//         uintptr_t key = reinterpret_cast<uintptr_t>((void*)&data);
//         this->nonBlockingCommParamTable_[key] = NonBlockingCommParam(NonBlockingCommParam::SEND,
//                                                                      requests);
//     }

//     return answer;
// }


// iReceiveDataX (array) =======================================================
template<typename T>
int TlCommunicate::iReceiveDataX(T* pData, const MPI_Datatype mpiType,
                                 const std::size_t start, const std::size_t end,
                                 const int src, const int tag)
{
    MPI_Request* pRequest = new MPI_Request;
    const int answer = MPI_Irecv((void*)(pData + start), (end - start),
                                 mpiType, src, tag,
                                 MPI_COMM_WORLD, pRequest);
    
    if (answer == 0) {
        std::vector<uintptr_t> requests;
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        const uintptr_t key = reinterpret_cast<uintptr_t>((void*)pData);
        const NonBlockingCommParam param(requests, tag,
                                         NonBlockingCommParam::RECV);
        this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
        {
            this->nonBlockingCommParamTable_[key] = param;
        }
    }

    return answer;
}


int TlCommunicate::iReceiveDataX(int* pData, const std::size_t size,
                                 const int src, const int tag)
{
    return this->iReceiveDataX(pData, MPI_INT, 0, size, src, tag);
}


int TlCommunicate::iReceiveDataX(double* pData, const std::size_t size,
                                 const int src, const int tag)
{
    return this->iReceiveDataX(pData, MPI_DOUBLE, 0, size, src, tag);
}


// iReceiveDataFromAnySource ===================================================

int TlCommunicate::iReceiveDataFromAnySource(int& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_INT, tag);
}


int TlCommunicate::iReceiveDataFromAnySource(unsigned int& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_UNSIGNED, tag);
}


int TlCommunicate::iReceiveDataFromAnySource(long& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_LONG, tag);
}


int TlCommunicate::iReceiveDataFromAnySource(unsigned long& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_UNSIGNED_LONG, tag);
}


template<typename T>
int TlCommunicate::iReceiveDataFromAnySource(T& data,
                                             const MPI_Datatype mpiType,
                                             const int tag)
{
    return this->iReceiveData(data, mpiType, MPI_ANY_SOURCE, tag);
}


int TlCommunicate::iReceiveDataFromAnySource(std::vector<int>& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_INT, tag);
}


int TlCommunicate::iReceiveDataFromAnySource(std::vector<double>& data, const int tag)
{
    return this->iReceiveDataFromAnySource(data, MPI_DOUBLE, tag);
}


template<typename T>
int TlCommunicate::iReceiveDataFromAnySource(std::vector<T>& data,
                                             const MPI_Datatype mpiType,
                                             const int tag)
{
    std::size_t size = 0;
    this->iReceiveDataFromAnySource(size, tag);
    int src = 0;
    this->wait(size, &src);
    data.resize(size);
    
    return this->iReceiveData(data, mpiType, 0, size, src, tag);
}


// =============================================================================
template<typename T>
int TlCommunicate::iReceiveDataFromAnySourceX(T* pData, const MPI_Datatype mpiType,
                                              const std::size_t start, const std::size_t end,
                                              const int tag)
{
    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T));
    const ldiv_t tmp = std::ldiv(length, bufCount);

    MPI_Request* pRequest = new MPI_Request;
    const int answer = MPI_Irecv((void*)(pData + start), (end - start), mpiType,
                                 MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, pRequest);

    if (answer == 0) {
        std::vector<uintptr_t> requests;
        requests.push_back(reinterpret_cast<uintptr_t>((void*)pRequest));
        const uintptr_t key = reinterpret_cast<uintptr_t>((void*)pData);
        const NonBlockingCommParam param(requests, tag,
                                         NonBlockingCommParam::RECV);
        this->checkNonBlockingTableCollision(key, param, __LINE__);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
        {
            this->nonBlockingCommParamTable_[key] = param;
        }
    }
    
    return answer;
}


int TlCommunicate::iReceiveDataFromAnySourceX(int* pData, const std::size_t size,
                                              const int tag)
{
    return this->iReceiveDataFromAnySourceX(pData, MPI_INT, 0, size, tag);
}


// =============================================================================
bool TlCommunicate::cancel(void* pData)
{
    bool answer = false;
    
    uintptr_t key = reinterpret_cast<uintptr_t>(pData);
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        NonBlockingCommParamTableType::iterator it = this->nonBlockingCommParamTable_.find(key);
        if (it != this->nonBlockingCommParamTable_.end()) {
            bool isComplete = true;
            std::vector<uintptr_t>::iterator reqEnd = it->second.requests.end();
            for (std::vector<uintptr_t>::iterator req = it->second.requests.begin(); req != reqEnd; ++req) {
                MPI_Request* pRequest = reinterpret_cast<MPI_Request*>(*req);
                const int retval = MPI_Cancel(pRequest);
                if (retval != MPI_SUCCESS) {
                    isComplete = false;
                    break;
                }
            }
            answer = isComplete;
        } else {
            std::cerr << "PROGRAM ERROR(TlCommunicate::cancel()): unknown key. file:" << __FILE__ << ", l." << __LINE__ << std::endl;
            std::abort();
        }
    }
    if (answer == true) {
        this->wait(pData);
    }
    
    return answer;
}


bool TlCommunicate::test(void* pData, int* pSrc)
{
    this->time_test_.start();
    ++(this->counter_test_);

    bool answer = false;
    
    uintptr_t key = reinterpret_cast<uintptr_t>(pData);

    NonBlockingCommParamTableType::iterator it;
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        it = this->nonBlockingCommParamTable_.find(key);
        if (it != this->nonBlockingCommParamTable_.end()) {
            const NonBlockingCommParam& param = it->second;
            if ((param.property & NonBlockingCommParam::COMPLETE) != 0) {
                answer = true;
                if (pSrc != NULL) {
                    *pSrc = param.source;
                }
            } else {
                bool isComplete = true;
                int source = -1;

                const std::size_t numOfRequests = param.requests.size();
                assert(numOfRequests == param.requestStates.size());
                for (std::size_t i = 0; i < numOfRequests; ++i) {
                    const unsigned int state = param.requestStates[i];
                    if (state == 0) {
                        MPI_Request* pRequest = reinterpret_cast<MPI_Request*>(param.requests[i]);
                        int flag = 0;
                        MPI_Status status;
                        MPI_Test(pRequest, &flag, &status);
                        if (flag == 0) {
                            isComplete = false;
                        } else {
                            it->second.requestStates[i] = 1;
                            source = status.MPI_SOURCE;
                        }
                    }
                }
                
                if (isComplete == true) {
                    it->second.property |= NonBlockingCommParam::COMPLETE;
                    it->second.source = source;
                    if (pSrc != NULL) {
                        *pSrc = source;
                    }
                }
                answer = isComplete;
            }
        } else {
            std::cerr << "PROGRAM ERROR(TlCommunicate::test()): unknown key. file:" << __FILE__ << ", l." << __LINE__ << std::endl;
            std::abort();
        }
    }

    this->time_test_.stop();
    return answer;
}


int TlCommunicate::wait(void* pData, int* pSrc)
{
    this->time_wait_.start();
    ++(this->counter_wait_);

    int answer = 0;

    uintptr_t key = reinterpret_cast<uintptr_t>(pData);
    MPI_Status status;

    NonBlockingCommParamTableType::iterator it;
#pragma omp critical (TlCommunicate_nonBlockingCommParamTable_update)
    {
        it = this->nonBlockingCommParamTable_.find(key);

        if (it != this->nonBlockingCommParamTable_.end()) {
            std::vector<uintptr_t>::iterator reqEnd = it->second.requests.end();
            for (std::vector<uintptr_t>::iterator req = it->second.requests.begin(); req != reqEnd; ++req) {
                MPI_Request* pRequest = reinterpret_cast<MPI_Request*>(*req);
                
                int err = MPI_Wait(pRequest, &status);
                if (err != 0) {
                    std::cerr << " MPI wait error. " << __FILE__ << ":" << __LINE__
                              << " err=" << err
                              << std::endl;
                }
                answer |= err;
                
                delete pRequest;
                pRequest = NULL;
            }
            
            this->nonBlockingCommParamTable_.erase(it);
            
            if (pSrc != NULL) {
                *pSrc = status.MPI_SOURCE;
            }
        }
    }
    
    this->time_wait_.stop();
    return answer;
}


// =====================================================================
// private
TlCommunicate::TlCommunicate()
    : workMemSize_(DEFAULT_WORK_MEM_SIZE), log_(TlLogging::getInstance())
{
}

TlCommunicate::TlCommunicate(const TlCommunicate& rhs)
    : log_(TlLogging::getInstance())
{
}

TlCommunicate::~TlCommunicate()
{
    // 後始末
    (void)this->finalize();
}

int TlCommunicate::initialize(int argc, char* argv[])
{
    int nAnswer = MPI_Init(&argc, &argv);

    // nproc & rank
    MPI_Comm_size(MPI_COMM_WORLD, &(this->m_nProc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(this->m_nRank));

    this->counter_barrier_ = 0;
    this->counter_test_ = 0;
    this->counter_wait_ = 0;
    this->counter_allreduce_ = 0;
    this->time_barrier_.stop();
    this->time_test_.stop();
    this->time_wait_.stop();
    this->time_allreduce_.stop();
    this->time_barrier_.reset();
    this->time_test_.reset();
    this->time_wait_.reset();
    this->time_allreduce_.reset();
    
    return nAnswer;
}

int TlCommunicate::finalize()
{
    //this->barrier();
    
    TlScalapackContext::finalize();

    return MPI_Finalize();
}


int TlCommunicate::abort(const int errorCode)
{
    int answer = MPI_Abort(MPI_COMM_WORLD, errorCode);
    this->finalize();
    return answer;
}


int TlCommunicate::getNumOfProcs() const
{
    return this->m_nProc;
}

int TlCommunicate::getRank() const
{
    return this->m_nRank;
}

/**
 *  rank が 0 だったらtrueを返す
 */
bool TlCommunicate::isMaster() const
{
    return (this->m_nRank == 0);
}

bool TlCommunicate::isSlave() const
{
    return !(this->isMaster());
}

// =====================================================================
// BROADCAST
// =====================================================================

int TlCommunicate::broadcast(bool& rData)
{
    int nData = 0;

    if (this->isMaster() == true) {
        nData = (rData == true) ? 1 : 0;
    }

    int nError = this->broadcast(nData);

    rData = ((nData == 1) ? true : false);

    return nError;
}

template<typename T>
int TlCommunicate::broadcast(T& data, const MPI_Datatype mpiType, const int root)
{
    this->log_.debug(TlUtils::format("TlCommunicate::broadcast(): type=%s, root=%d",
                                     TlUtils::xtos(mpiType).c_str(), root));
    return MPI_Bcast(&data, 1, mpiType, root, MPI_COMM_WORLD);
}

int TlCommunicate::broadcast(int& data, int root)
{
    return this->broadcast(data, MPI_INT, root);
}

int TlCommunicate::broadcast(unsigned int& data, int root)
{
    return this->broadcast(data, MPI_UNSIGNED, root);
}

int TlCommunicate::broadcast(long& data, int root)
{
    return this->broadcast(data, MPI_LONG, root);
}

int TlCommunicate::broadcast(unsigned long& data, int root)
{
    return this->broadcast(data, MPI_UNSIGNED_LONG, root);
}

int TlCommunicate::broadcast(double& data, int root)
{
    return this->broadcast(data, MPI_DOUBLE, root);
}

int TlCommunicate::broadcast(std::string& rData)
{
    // バッファサイズの送受信
    int nCount = 0;
    if (this->isMaster() == true) {
        nCount = rData.size();
    }
    this->broadcast(nCount);

    char* pBuffer = new char[nCount];
    if (this->isMaster() == true) {
        rData.copy(pBuffer, nCount);
    }
    const int nRoot = 0;

    int nError = MPI_Bcast(pBuffer, nCount, MPI_CHAR, nRoot, MPI_COMM_WORLD);

    rData = std::string(pBuffer, nCount);

    delete[] pBuffer;
    pBuffer = NULL;

    return nError;
}

template<typename T>
int TlCommunicate::broadcast(std::vector<T>& data, const MPI_Datatype mpiType,
                             const std::size_t start, const std::size_t end,
                             int root)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(T)); // T型配列の配列数
    std::vector<T> buf(bufCount);
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Bcast(&(data[startIndex]), bufCount, mpiType,
                           root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Bcast(&(data[startIndex]), remain, mpiType,
                           root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}

int TlCommunicate::broadcast(std::vector<int>& rData, int nRoot)
{
    std::size_t size = 0;
    if (this->getRank() == nRoot) {
        size = rData.size();
    }
    this->broadcast(size, nRoot);
    rData.resize(size);

    return this->broadcast(rData, MPI_INT, 0, size, nRoot);
}

int TlCommunicate::broadcast(std::vector<long>& data, int root)
{
    std::size_t size = 0;
    if (this->isMaster() == true) {
        size = data.size();
    }
    this->broadcast(size);
    data.resize(size);

    return this->broadcast(data, MPI_LONG, 0, size, root);
}

int TlCommunicate::broadcast(std::vector<unsigned long>& rData, const int nRoot)
{
    std::size_t size = 0;
    if (this->getRank() == nRoot) {
        size = rData.size();
    }
    this->broadcast(size, nRoot);
    rData.resize(size);

    return this->broadcast(rData, MPI_UNSIGNED_LONG, 0, size, nRoot);
}

int TlCommunicate::broadcast(std::vector<double>& rData, const int nRoot)
{
    std::size_t size = 0;
    if (this->getRank() == nRoot) {
        size = rData.size();
    }
    this->broadcast(size, nRoot);
    rData.resize(size);

    return this->broadcast(rData, MPI_DOUBLE, 0, size, nRoot);
}


int TlCommunicate::broadcast(std::valarray<double>& data,
                             const std::size_t start, const std::size_t end,
                             int root)
{
    int answer = 0;
    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Bcast((void*)(&data[startIndex]), bufCount, MPI_DOUBLE, root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Bcast((void*)(&data[startIndex]), remain, MPI_DOUBLE, root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}

int TlCommunicate::broadcast(std::valarray<double>& rData, const int nRoot)
{
    std::size_t size = 0;
    if (this->getRank() == nRoot) {
        size = rData.size();
    }
    int answer = this->broadcast(size, nRoot);

    /// この排他処理はvalarrayのために必要
    if (this->getRank() != nRoot) {
        rData.resize(size);
    }

    if (answer == 0) {
        answer = this->broadcast(rData, 0, size, nRoot);
    }

    return answer;
}

int TlCommunicate::broadcast(std::vector<std::string>& rData)
{
    int nArraySize = 0;
    if (this->isMaster() == true) {
        nArraySize = rData.size();
    }
    this->broadcast(nArraySize);
    rData.resize(nArraySize);

    for (int i = 0; i < nArraySize; ++i) {
        std::string str = "";
        if (this->isMaster() == true) {
            str = rData[i];
        }
        this->broadcast(str);

        rData[i] = str;
    }

    return 0;
}

int TlCommunicate::broadcast(TlVector& data, const int root)
{
    std::size_t size = 0;
    if (this->getRank() == root) {
        size = data.getSize();
    }
    int answer = this->broadcast(size, root);
    data.resize(size);

    if (answer == 0) {
        answer = this->broadcast(data.data_, MPI_DOUBLE, 0, data.getSize(), root);
    }

    return answer;
}

int TlCommunicate::broadcast(TlMatrix& data, const int root)
{
    TlMatrixObject::index_type row = 0;
    TlMatrixObject::index_type col = 0;

    if (this->getRank() == root) {
        row = data.getNumOfRows();
        col = data.getNumOfCols();
    }
    int answer = 0;
    answer = this->broadcast(row, root);
    answer = this->broadcast(col, root);
    data.resize(row, col);

    if (answer == 0) {
        this->broadcast(data.data_, MPI_DOUBLE, 0, data.getNumOfElements(), root);
    }

    return answer;
}

int TlCommunicate::broadcast(TlSymmetricMatrix& data, int root)
{
    TlMatrixObject::index_type dim = 0;

    if (this->getRank() == root) {
        dim = data.getNumOfRows();
    }
    int answer = this->broadcast(dim, root);
    data.resize(dim);

    if (answer == 0) {
        answer = this->broadcast(data.data_, MPI_DOUBLE, 0, data.getNumOfElements(), root);
    }

    return answer;
}

int TlCommunicate::broadcast(TlSparseMatrix& rData, const int nRoot)
{
    std::size_t nRow = 0;
    std::size_t nCol = 0;

    if (this->getRank() == nRoot) {
        nRow = rData.getNumOfRows();
        nCol = rData.getNumOfCols();
    }
    this->broadcast(nRow, nRoot);
    this->broadcast(nCol, nRoot);

    std::vector<unsigned long> indexData;
    std::vector<double> valData;
    if (this->getRank() == nRoot) {
        TlSparseMatrix::const_iterator pEnd = rData.end();
        for (TlSparseMatrix::const_iterator p = rData.begin(); p != pEnd; ++p) {
            indexData.push_back(p->first);
            valData.push_back(p->second);
        }
    }
    this->broadcast(indexData, nRoot);
    this->broadcast(valData, nRoot);

    //if (this->getRank() == nRoot){
    TlSparseMatrix m(nRow, nCol);
    const std::size_t nSize= indexData.size();
    for (std::size_t i = 0; i < nSize; ++i) {
        const unsigned long index = indexData[i];
        const double val = valData[i];
        //m(row, col) = val;
        m.m_aMatrix[index] = val;
    }
    rData = m;
    //}

    return 0;
}

// int TlCommunicate::broadcast(TlSparseHashMatrix& rData, const int nRoot) const{
//   int nRow = 0;
//   int nCol = 0;

//   if (this->getRank() == nRoot){
//     nRow = rData.getNumOfRows();
//     nCol = rData.getNumOfCols();
//   }
//   this->broadcast(nRow, nRoot);
//   this->broadcast(nCol, nRoot);

//   std::vector<unsigned long> indexData;
//   std::vector<double> valData;
//   if (this->getRank() == nRoot){
//     TlSparseHashMatrix::const_iterator pEnd = rData.end();
//     for (TlSparseHashMatrix::const_iterator p = rData.begin(); p != pEnd; ++p){
//       indexData.push_back(p->first);
//       valData.push_back(p->second);
//     }
//   }
//   this->broadcast(indexData, nRoot);
//   this->broadcast(valData, nRoot);

//   TlSparseHashMatrix m(nRow, nCol);
//   const int nSize= static_cast<int>(indexData.size());
//   for (int i = 0; i < nSize; ++i){
//     const unsigned long index = indexData[i];
//     const double val = valData[i];
//     m.container[index] = val;
//   }
//   rData = m;

//   return 0;
// }


// int TlCommunicate::broadcast(TlMmapMatrix& data, int root)
// {
//     std::size_t row = 0;
//     std::size_t col = 0;
//     if (this->getRank() == root) {
//         row = data.getNumOfRows();
//         col = data.getNumOfCols();
//     }
//     int answer = 0;
//     answer = this->broadcast(row, root);
//     answer = this->broadcast(col, root);
//     if (this->getRank() != root) {
//         data.resize(row, col);
//     }

//     if (answer == 0) {
//         this->broadcast(data.data_, 0, data.getNumOfElements(), root);
//     }

//     return answer;
// }


// int TlCommunicate::broadcast(TlMmapSymmetricMatrix& data, int root)
// {
//     std::size_t dim = 0;
//     if (this->getRank() == root) {
//         dim = data.getNumOfRows();
//     }
//     int answer = 0;
//     answer = this->broadcast(dim, root);
//     if (this->getRank() != root) {
//         data.resize(dim);
//     }

//     if (answer == 0) {
//         this->broadcast(data.data_, 0, data.getNumOfElements(), root);
//     }

//     return answer;
// }

template<typename T>
int TlCommunicate::broadcast(T* pData, const MPI_Datatype mpiType,
                             const std::size_t start, const std::size_t end,
                             int root)
{
    int answer = 0;

    const long length = static_cast<long>(end - start);
    const int bufCount = static_cast<int>(this->workMemSize_ / sizeof(double));
    const ldiv_t tmp = std::ldiv(length, bufCount);

    // 作業用メモリ分のループ
    for (long i = 0; i < tmp.quot; ++i) {
        const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * i);
        answer = MPI_Bcast((void*)(pData + startIndex), bufCount, mpiType, root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    // 残り分のループ
    if (tmp.rem != 0) {
        const int remain = tmp.rem;
        const std::size_t startIndex = start + static_cast<std::size_t>(bufCount * tmp.quot);
        answer = MPI_Bcast((void*)(pData + startIndex), remain, mpiType, root, MPI_COMM_WORLD);
        if (answer != 0) {
            std::cerr << " MPI error. " << __FILE__ << ":" << __LINE__
                      << " anwer=" << answer
                      << std::endl;
        }
    }

    return answer;
}


int TlCommunicate::broadcast(double* pData, const std::size_t size, int root)
{
    return this->broadcast(pData, MPI_DOUBLE, 0, size, root);
}


int TlCommunicate::broadcast(TlParameter& rParam)
{
    std::vector<std::string> aGroups;
    if (this->isMaster() == true) {
        aGroups = rParam.getGroups();
    }
    this->broadcast(aGroups);

    const int nNumOfGroups = aGroups.size();
    for (int nGroup = 0; nGroup < nNumOfGroups; ++nGroup) {
        const std::string sGroup = aGroups[nGroup];

        std::vector<std::string> aKeywords;
        if (this->isMaster() == true) {
            aKeywords = rParam[sGroup].getKeywords();
        }
        this->broadcast(aKeywords);

        const int nNumOfKeywords = aKeywords.size();
        for (int nKeyword = 0; nKeyword < nNumOfKeywords; ++nKeyword) {
            const std::string sKeyword = aKeywords[nKeyword];

            std::string sValue= "";
            if (this->isMaster() == true) {
                sValue = rParam[sGroup][sKeyword];
            }
            this->broadcast(sValue);

            rParam[sGroup][sKeyword] = sValue;
        }
    }

    return 0;
}


int TlCommunicate::broadcast(TlSerializeData& data)
{
    std::string buf = "";
    if (this->isMaster() == true) {
        TlMsgPack mpac(data);
        buf = mpac.dump();
    }

    int ans = this->broadcast(buf);

    if ((ans == 0) && (this->isMaster() != true)) {
        TlMsgPack mpac;
        mpac.pack(buf);
        data = mpac.getSerializeData();
    }

    return ans;
}


int TlCommunicate::allReduce_SUM(const TlFileMatrix& fromLocalMatrix,
                                 const std::string& toMatrixFilePath)
{
    const std::size_t numOfRows = fromLocalMatrix.getNumOfRows();
    const std::size_t numOfCols = fromLocalMatrix.getNumOfCols();

    const std::size_t dataSize = numOfRows * numOfCols;
    std::size_t maxBufferIndex = 10 * 1024 * 1024 / sizeof(double); // 10 MB分
    std::vector<double> buf(maxBufferIndex);

    const std::size_t startFromPos = fromLocalMatrix.startPos_;
    fromLocalMatrix.fs_.seekg(static_cast<std::fstream::pos_type>(startFromPos), std::ios_base::beg);

    size_t currentPos = 0;

    if (this->isMaster() == true) {
        // master用ルーチン
        // 書き込みルーチンが含まれている点がslaveと異なる
        TlFileMatrix toMatrix(toMatrixFilePath, numOfRows, numOfCols);
        const std::size_t startToPos = toMatrix.startPos_;
        toMatrix.fs_.seekp(static_cast<std::fstream::pos_type>(startToPos), std::ios_base::beg);

        while (currentPos < dataSize) {
            const size_t readSize = std::min(maxBufferIndex, dataSize - currentPos);
            const size_t bufferSize = sizeof(double) * readSize;
            fromLocalMatrix.fs_.read(reinterpret_cast<char*>(&(buf[0])), bufferSize);
            this->allReduce_SUM(buf);
            currentPos += readSize;

            toMatrix.fs_.write(reinterpret_cast<const char*>(&(buf[0])), bufferSize);
        }
    } else {
        while (currentPos < dataSize) {
            const size_t readSize = std::min(maxBufferIndex, dataSize - currentPos);
            const size_t bufferSize = sizeof(double) * readSize;
            fromLocalMatrix.fs_.read(reinterpret_cast<char*>(&(buf[0])), bufferSize);
            this->allReduce_SUM(buf);
            currentPos += readSize;
        }
    }

    return 0;
}


int TlCommunicate::allReduce_SUM(const TlFileSymmetricMatrix& fromLocalMatrix,
                                 const std::string& toMatrixFilePath)
{
    const std::size_t numOfDims = fromLocalMatrix.getNumOfRows();
    const std::size_t dataSize = numOfDims * (numOfDims +1) / 2;
    std::size_t maxBufferIndex = 10 * 1024 * 1024 / sizeof(double); // 10 MB分
    std::vector<double> buf(maxBufferIndex);

    const std::size_t startFromPos = fromLocalMatrix.startPos_;
    fromLocalMatrix.fs_.seekg(static_cast<std::fstream::pos_type>(startFromPos), std::ios_base::beg);

    size_t currentPos = 0;

    if (this->isMaster() == true) {
        // master用ルーチン
        // 書き込みルーチンが含まれている点がslaveと異なる
        TlFileSymmetricMatrix toMatrix(toMatrixFilePath, numOfDims);
        const std::size_t startToPos = toMatrix.startPos_;
        toMatrix.fs_.seekp(static_cast<std::fstream::pos_type>(startToPos), std::ios_base::beg);

        while (currentPos < dataSize) {
            const size_t readSize = std::min(maxBufferIndex, dataSize - currentPos);
            const size_t bufferSize = sizeof(double) * readSize;
            fromLocalMatrix.fs_.read(reinterpret_cast<char*>(&(buf[0])), bufferSize);
            this->allReduce_SUM(buf);
            currentPos += readSize;

            toMatrix.fs_.write(reinterpret_cast<const char*>(&(buf[0])), bufferSize);
        }
    } else {
        while (currentPos < dataSize) {
            const size_t readSize = std::min(maxBufferIndex, dataSize - currentPos);
            const size_t bufferSize = sizeof(double) * readSize;
            fromLocalMatrix.fs_.read(reinterpret_cast<char*>(&(buf[0])), bufferSize);
            this->allReduce_SUM(buf);
            currentPos += readSize;
        }
    }

    return 0;
}

