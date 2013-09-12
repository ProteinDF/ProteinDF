#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfCD_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "DfOverlapEngine.h"
#include "DfEriEngine.h"
#include "TlCommunicate.h"
#include "TlTime.h"

#define TRANS_MEM_SIZE (1 * 1024 * 1024 * 1024) // 1GB
// #define CD_DEBUG

DfCD_Parallel::DfCD_Parallel(TlSerializeData* pPdfParam) 
    : DfCD(pPdfParam) {

    this->isDebugSaveL_ = false;
    if (! (*this->pPdfParam_)["debug/saveL"].getStr().empty()) {
        this->isDebugSaveL_ = (*this->pPdfParam_)["debug/saveL"].getBoolean();
    }
}


DfCD_Parallel::~DfCD_Parallel()
{
}

void DfCD_Parallel::calcCholeskyVectorsForJK()
{
    // for J & K
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets"]);
    // productive code
    const TlRowVectorMatrix2 Ljk 
        = DfCD::calcCholeskyVectorsOnTheFly<DfEriEngine>(orbInfo,
                                                         this->getI2pqVtrPath());
    this->saveL(Ljk, DfObject::getLjkMatrixPath());
}

void DfCD_Parallel::calcCholeskyVectorsForGridFree()
{
    // for XC(gridfree)
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);

    if (this->isDedicatedBasisForGridFree_) {
        const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                      (*this->pPdfParam_)["basis_sets_GF"]);
        // productive code
        const TlRowVectorMatrix2 Lxc 
            = DfCD::calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                 orbInfo_q,
                                                                 this->getI2pqVtrXCPath());
        this->saveL(Lxc, DfObject::getLxcMatrixPath());
    } else {
        // productive code
        this->log_.info("build Lxc matrix by on-the-fly method.");
        const TlRowVectorMatrix2 Lxc 
            = DfCD::calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                 this->getI2pqVtrXCPath());
        this->saveL(Lxc, DfObject::getLxcMatrixPath());
    }
}

DfTaskCtrl* DfCD_Parallel::getDfTaskCtrlObject() const
{
    DfTaskCtrl *pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}

void DfCD_Parallel::finalize(TlSymmetricMatrix *pMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMat);
}

void DfCD_Parallel::finalize(TlSparseSymmetricMatrix *pMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.gatherToMaster(*pMat);
    rComm.broadcast(*pMat);
}


void DfCD_Parallel::finalize_I2PQ(PQ_PairArray* pI2PQ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    
    // gather to master
    if (rComm.isMaster() == true) {
        const int slaves = numOfProcs - 1;
        int proc = 0;
        std::vector<index_type> shellArray;
        for (int i = 0; i < slaves; ++i) {
            rComm.receiveDataFromAnySource(shellArray, &proc);

            const std::size_t I2PQ_size = shellArray.size() / 2;
            PQ_PairArray i2pq_tmp(I2PQ_size);
            for (std::size_t i = 0; i < I2PQ_size; ++i) {
                i2pq_tmp[i] = Index2(shellArray[i*2   ],
                                     shellArray[i*2 +1]);
            }
            
            pI2PQ->insert(pI2PQ->end(),
                          i2pq_tmp.begin(), i2pq_tmp.end());
        }
        //std::sort(pI2PQ->begin(), pI2PQ->end(), PQ_Pair_less());
        std::sort(pI2PQ->begin(), pI2PQ->end());
    } else {
        const std::size_t I2PQ_size = pI2PQ->size();
        std::vector<index_type> shellArray(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i*2   ] = (*pI2PQ)[i].index1();
            shellArray[i*2 +1] = (*pI2PQ)[i].index2();
        }
        rComm.sendData(shellArray);
    }
    
    // broadcast
    {
        std::vector<index_type> tmp;
        if (rComm.isMaster() == true) {
            const std::size_t size = pI2PQ->size();
            tmp.resize(size * 2);
            for (std::size_t i = 0; i < size; ++i) {
                tmp[i*2   ] = (*pI2PQ)[i].index1();
                tmp[i*2 +1] = (*pI2PQ)[i].index2();
            }
        }
        rComm.broadcast(tmp);
        if (rComm.isMaster() != true) {
            const std::size_t size = tmp.size() / 2;
            pI2PQ->resize(size);
            for (std::size_t i = 0; i < size; ++i) {
                (*pI2PQ)[i] = Index2(tmp[i*2   ],
                                     tmp[i*2 +1]);
            }
        }
    }
}


void DfCD_Parallel::saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCD::saveI2PQ(I2PQ, filepath);
    }
}


DfCD::PQ_PairArray DfCD_Parallel::getI2PQ(const std::string& filepath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("distribute I2PQ table.");

    PQ_PairArray I2PQ;
    std::vector<index_type> shellArray;
    if (rComm.isMaster() == true) {
        I2PQ = DfCD::getI2PQ(filepath);

        const std::size_t I2PQ_size = I2PQ.size(); 
        shellArray.resize(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i*2   ] = I2PQ[i].index1();
            shellArray[i*2 +1] = I2PQ[i].index2();
        }
    }
    
    rComm.broadcast(shellArray);
    if (rComm.isMaster() != true) {
        const std::size_t I2PQ_size = shellArray.size() / 2;
        I2PQ.resize(I2PQ_size);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            I2PQ[i] = Index2(shellArray[i*2   ],
                             shellArray[i*2 +1]);
        }
    }

    return I2PQ;
}

// void DfCD_Parallel::saveLjk(const TlMatrix& L)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfCD::saveLjk(L);
//     }
// }


// TlMatrix DfCD_Parallel::getLjk()
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     TlMatrix L;
//     if (rComm.isMaster() == true) {
//         L = DfCD::getLjk();
//     }
//     rComm.broadcast(L);

//     return L;
// }


TlSymmetricMatrix DfCD_Parallel::getPMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlMatrix P;
    if (rComm.isMaster() == true) {
        P = DfCD::getPMatrix();
    }
    rComm.broadcast(P);
    return P;
}


void DfCD_Parallel::divideCholeskyBasis(const index_type numOfCBs,
                                        index_type *pStart, index_type *pEnd)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();
    
    const index_type range = (numOfCBs + numOfProcs -1) / numOfProcs;
    *pStart = range * rank;
    *pEnd = std::min(range * (rank +1), numOfCBs);
}


TlDistributeSymmetricMatrix 
DfCD_Parallel::getCholeskyVector_distribute(const TlVector& L_col,
                                            const PQ_PairArray& I2PQ)
{
    // TlCommunicate& rComm = TlCommunicate::getInstance();
    const index_type numOfItilde = L_col.getSize();
    TlDistributeSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        // std::cerr << TlUtils::format("[%d] L_col[%d] => LI(%d, %d)",
        //                              rComm.getRank(),
        //                              i, I2PQ[i].shellIndex1, I2PQ[i].shellIndex2)
        //           << std::endl;
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

TlRowVectorMatrix2 DfCD_Parallel::calcCholeskyVectorsOnTheFlyS(const TlOrbitalInfoObject& orbInfo,
                                                               const std::string& I2PQ_path)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("call on-the-fly Cholesky Decomposition routine (MPI parallel)");
    assert(this->pEngines_ != NULL);
    this->initializeCutoffStats(orbInfo.getMaxShellType());

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs +1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector global_diagonals; // 対角成分

    assert(this->pEngines_ != NULL);
    this->calcDiagonals(orbInfo, &I2PQ, &schwartzTable, &global_diagonals);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();

    // prepare variables (parallel)
    TlRowVectorMatrix2 L(N, 1,
                         rComm.getNumOfProcs(),
                         rComm.getRank(),
                         this->isEnableMmap_); // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
    const index_type local_N = L.getNumOfLocalRows();
    std::vector<double> L_pm(N);
    std::vector<int> global_pivot(N); // ERI計算リストを作成するために必要
    std::vector<int> reverse_pivot(N); // global_pivotの逆引き
    std::vector<int> local_pivot(local_N);
    std::vector<double> local_diagonals(local_N);
    const int myRank = rComm.getRank();

    double error = 0.0;
    index_type error_global_loc = 0;
    index_type error_local_loc = 0;
    {
        int local_i = 0;
        for (int global_i = 0; global_i < N; ++global_i) {
            global_pivot[global_i] = global_i;
            reverse_pivot[global_i] = global_i;
            if (L.getPEinChargeByRow(global_i) == myRank) {
                local_pivot[local_i] = global_i;
                local_diagonals[local_i] = global_diagonals[global_i];
                if (error < local_diagonals[local_i]) {
                    error_local_loc = local_i;
                    error_global_loc = local_pivot[local_i];
                    error = local_diagonals[local_i];
                }
                ++local_i;
            }
        }
        assert(local_i == local_N);
    }
    rComm.allReduce_MAXLOC(&error, &error_global_loc);

    index_type m = 0;
    index_type local_m = 0;
    if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
        std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
        ++local_m;
    }

    int progress = 0;
    index_type division = std::max<index_type>(N * 0.01, 100);
    while (error > threshold) {
#ifdef CD_DEBUG
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif // CD_DEBUG

        // progress 
        //CD_resizeL_time.start();
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e",
                                            m, error));
            ++progress;

            // メモリの確保
            L.reserve_cols(progress * division);
        }
        L.resize(N, m+1);
        //CD_resizeL_time.stop();

        // pivot
        std::swap(global_pivot[m], global_pivot[error_global_loc]);
        reverse_pivot[global_pivot[m]] = m;
        reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

        const double l_m_pm = std::sqrt(error);
        const index_type pivot_m = global_pivot[m];
        L.set(pivot_m, m, l_m_pm); // 通信発生せず。関係無いPEは値を捨てる。
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        std::vector<double> G_pm;
        // const index_type numOf_G_cols = N -(m+1);
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = global_pivot[m+1 +i]; // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(orbInfo,
                                                pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        {
            // 全PEに分配
            const int PEinCharge = L.getPEinChargeByRow(pivot_m);
#ifndef NDEBUG
            if (PEinCharge == rComm.getRank()) {
                const index_type copySize = L.getRowVector(pivot_m, &(L_pm[0]), m +1);
                assert(copySize == m +1);
            }
#endif // NDEBUG
            rComm.broadcast(&(L_pm[0]), m +1, PEinCharge);
        }

        error = 0.0;
#pragma omp parallel
        {
            std::vector<double> L_pi(m +1);
            double my_error = 0.0;
            int my_error_global_loc = 0;
            int my_error_local_loc = 0;

#pragma omp for schedule(runtime)
            for (int i = local_m; i < local_N; ++i) {
                const int pivot_i = local_pivot[i];
                const index_type copySize = L.getRowVector(pivot_i, &(L_pi[0]), m +1);
                assert(copySize == m +1);
                double sum_ll = 0.0;
                for (index_type j = 0; j < m; ++j) {
                    sum_ll += L_pm[j] * L_pi[j];
                }

                const int G_pm_index = reverse_pivot[pivot_i] - (m+1);
                const double l_m_pi = (G_pm[G_pm_index] - sum_ll) * inv_l_m_pm;
#pragma omp critical(DfCD_Parallel__calcCholeskyVectors_onTheFly)
                {
                    L.set(pivot_i, m, l_m_pi);
                }

                const double ll = l_m_pi * l_m_pi;
#pragma omp atomic
                global_diagonals[pivot_i] -= ll;

                if (global_diagonals[pivot_i] > my_error) {
                    my_error = global_diagonals[pivot_i];
                    my_error_global_loc = reverse_pivot[pivot_i]; // == m +1 + i
                    my_error_local_loc = i;
                }
            }

#ifdef _OPENMP
            const int numOfThreads = omp_get_num_threads();
            const int myThreadID = omp_get_thread_num();
            for (int thread = 0; thread < numOfThreads; ++thread) {
                if (thread == myThreadID) {
                    if (error < my_error) {
                        error = my_error;
                        error_global_loc = my_error_global_loc;
                        error_local_loc = my_error_local_loc;
                    }
                }
#pragma omp flush(error, error_global_loc)
            }
#else
            error = my_error;
            error_global_loc = my_error_global_loc;
            error_local_loc = my_error_local_loc;
#endif // _OPENMP
        }

        rComm.allReduce_MAXLOC(&error, &error_global_loc);
        global_diagonals[global_pivot[error_global_loc]] = error;

        ++m;
        if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
            std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
            ++local_m;
        }
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(orbInfo.getMaxShellType());

    return L;
}


TlRowVectorMatrix2 DfCD_Parallel::calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
                                                               const TlOrbitalInfoObject& orbInfo_q,
                                                               const std::string& I2PQ_path)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("call on-the-fly Cholesky Decomposition routine (MPI parallel)");
    assert(this->pEngines_ != NULL);
    this->initializeCutoffStats(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfOrbs_p * numOfOrbs_q;
    this->log_.info(TlUtils::format("number of orbitals1: %d", numOfOrbs_p));
    this->log_.info(TlUtils::format("number of orbitals2: %d", numOfOrbs_q));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));
    PQ_PairArray I2PQ;
    TlSparseMatrix schwartzTable(numOfOrbs_p, numOfOrbs_q);
    TlVector global_diagonals; // 対角成分

    assert(this->pEngines_ != NULL);
    this->calcDiagonalsA(orbInfo_p, orbInfo_q,
                         &I2PQ, &schwartzTable, &global_diagonals);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();

    // prepare variables (parallel)
    TlRowVectorMatrix2 L(N, 1,
                         rComm.getNumOfProcs(),
                         rComm.getRank(),
                         this->isEnableMmap_); // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
    const index_type local_N = L.getNumOfLocalRows();
    std::vector<double> L_pm(N);
    std::vector<int> global_pivot(N); // ERI計算リストを作成するために必要
    std::vector<int> reverse_pivot(N); // global_pivotの逆引き
    std::vector<int> local_pivot(local_N);
    std::vector<double> local_diagonals(local_N);
    const int myRank = rComm.getRank();

    double error = 0.0;
    index_type error_global_loc = 0;
    index_type error_local_loc = 0;
    {
        int local_i = 0;
        for (int global_i = 0; global_i < N; ++global_i) {
            global_pivot[global_i] = global_i;
            reverse_pivot[global_i] = global_i;
            if (L.getPEinChargeByRow(global_i) == myRank) {
                local_pivot[local_i] = global_i;
                local_diagonals[local_i] = global_diagonals[global_i];
                if (error < local_diagonals[local_i]) {
                    error_local_loc = local_i;
                    error_global_loc = local_pivot[local_i];
                    error = local_diagonals[local_i];
                }
                ++local_i;
            }
        }
        assert(local_i == local_N);
    }
    rComm.allReduce_MAXLOC(&error, &error_global_loc);

    index_type m = 0;
    index_type local_m = 0;
    if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
        std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
        ++local_m;
    }

    int progress = 0;
    index_type division = std::max<index_type>(N * 0.01, 100);
    while (error > threshold) {
#ifdef CD_DEBUG
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif // CD_DEBUG

        // progress 
        //CD_resizeL_time.start();
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e",
                                            m, error));
            ++progress;

            // メモリの確保
            L.reserve_cols(progress * division);
        }
        L.resize(N, m+1);

        // pivot
        std::swap(global_pivot[m], global_pivot[error_global_loc]);
        reverse_pivot[global_pivot[m]] = m;
        reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

        const double l_m_pm = std::sqrt(error);
        const index_type pivot_m = global_pivot[m];
        L.set(pivot_m, m, l_m_pm); // 通信発生せず。関係無いPEは値を捨てる。
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // calc
        std::vector<double> G_pm;
        // const index_type numOf_G_cols = N -(m+1);
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = global_pivot[m+1 +i]; // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElementsA(orbInfo_p, orbInfo_q,
                                                 pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        {
            // 全PEに分配
            const int PEinCharge = L.getPEinChargeByRow(pivot_m);
            if (PEinCharge == rComm.getRank()) {
                const index_type copySize = L.getRowVector(pivot_m, &(L_pm[0]), m +1);
                assert(copySize == m +1);
            }
            rComm.broadcast(&(L_pm[0]), m +1, PEinCharge);
        }

        error = 0.0;
#pragma omp parallel
        {
            std::vector<double> L_pi(m +1);
            double my_error = 0.0;
            int my_error_global_loc = 0;
            int my_error_local_loc = 0;

#pragma omp for schedule(runtime)
            for (int i = local_m; i < local_N; ++i) {
                const int pivot_i = local_pivot[i];
                const index_type copySize = L.getRowVector(pivot_i, &(L_pi[0]), m +1);
                assert(copySize == m +1);
                double sum_ll = 0.0;
                for (index_type j = 0; j < m; ++j) {
                    sum_ll += L_pm[j] * L_pi[j];
                }

                const int G_pm_index = reverse_pivot[pivot_i] - (m+1);
                const double l_m_pi = (G_pm[G_pm_index] - sum_ll) * inv_l_m_pm;
#pragma omp critical(DfCD_Parallel__calcCholeskyVectors_onTheFly)
                {
                    L.set(pivot_i, m, l_m_pi);
                }

                const double ll = l_m_pi * l_m_pi;
#pragma omp atomic
                global_diagonals[pivot_i] -= ll;

                if (global_diagonals[pivot_i] > my_error) {
                    my_error = global_diagonals[pivot_i];
                    my_error_global_loc = reverse_pivot[pivot_i]; // == m +1 + i
                    my_error_local_loc = i;
                }
            }

#ifdef _OPENMP
            const int numOfThreads = omp_get_num_threads();
            const int myThreadID = omp_get_thread_num();
            for (int thread = 0; thread < numOfThreads; ++thread) {
                if (thread == myThreadID) {
                    if (error < my_error) {
                        error = my_error;
                        error_global_loc = my_error_global_loc;
                        error_local_loc = my_error_local_loc;
                    }
                }
#pragma omp flush(error, error_global_loc)
            }
#else
            error = my_error;
            error_global_loc = my_error_global_loc;
            error_local_loc = my_error_local_loc;
#endif // _OPENMP
        }

        rComm.allReduce_MAXLOC(&error, &error_global_loc);
        global_diagonals[global_pivot[error_global_loc]] = error;

        ++m;
        if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
            std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
            ++local_m;
        }
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    return L;
}


std::vector<double>
DfCD_Parallel::getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                                      const index_type G_row,
                                      const std::vector<index_type>& G_col_list,
                                      const PQ_PairArray& I2PQ,
                                      const TlSparseSymmetricMatrix& schwartzTable)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const std::size_t numOfG_cols = G_col_list.size();
    const std::size_t range = (numOfG_cols + numOfProcs -1) / numOfProcs;
    const std::size_t start = std::min(range * rank, numOfG_cols);
    const std::size_t end = std::min(range * (rank +1), numOfG_cols);

    std::vector<double> elements;
    if (end - start > 0) {
        std::vector<index_type> G_col_list_local(end - start);
        std::copy(G_col_list.begin() + start,
                  G_col_list.begin() + end,
                  G_col_list_local.begin());
        elements = DfCD::getSuperMatrixElements(orbInfo, 
                                                G_row,
                                                G_col_list_local,
                                                I2PQ,
                                                schwartzTable);
    }
    assert(elements.size() == (end - start));

    // gather to master
    const int TAG_SME_HEADER = 1001;
    const int TAG_SME_DATA   = 1002;
    std::vector<double> answer(numOfG_cols);
    std::vector<std::size_t> header(2);
    if (rank == 0) {
        // copy from myself
        std::copy(elements.begin(), elements.end(),
                  answer.begin());

        // receive data from slaves
        std::vector<bool> recvCheck(numOfProcs, false);
        std::vector<double> buf(range);
        for (int pe = 1; pe < numOfProcs; ++pe) {
            int src = 0;
            rComm.receiveDataFromAnySourceX(&(header[0]), 2, &src, TAG_SME_HEADER);

            const std::size_t local_start = header[0];
            const std::size_t local_end   = header[1];
            const std::size_t size = local_end - local_start;
            rComm.receiveDataX(&(buf[0]), size, src, TAG_SME_DATA);

            std::copy(buf.begin(), buf.begin() + size,
                      answer.begin() + local_start);
        }

    } else {
        header[0] = start;
        header[1] = end;
        rComm.sendDataX(&(header[0]), 2, 0, TAG_SME_HEADER);
        rComm.sendDataX(&(elements[0]), elements.size(), 0, TAG_SME_DATA);
    }

    //rComm.broadcast(answer);
    rComm.broadcast(&(answer[0]), numOfG_cols, 0);

    return answer;
}


void DfCD_Parallel::saveL(const TlRowVectorMatrix2& L,
                          const std::string& path)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix2 colVecL(numOfRows, numOfCols, numOfProcs, rank,
                               isUsingMemManager);

    const div_t turns = std::div(numOfRows, numOfProcs);
    const index_type localRows = turns.quot + 1;

    const std::size_t colMemSize = numOfCols * sizeof(double);
    const std::size_t transMemSize = TRANS_MEM_SIZE;
    const int transRowsPerCycle = std::max<int>(transMemSize / colMemSize, 1);
    const int transCycle = localRows / transRowsPerCycle + 1;
    std::vector<double> buf(numOfCols * transRowsPerCycle);
    for (int proc = 0; proc < numOfProcs; ++proc) {
        for (int cycle = 0; cycle < transCycle; ++cycle) {
            if (proc == rank) {
                std::vector<double> v(numOfCols);
                for (index_type r = 0; r < transRowsPerCycle; ++r) {
                    const index_type row = (cycle * transRowsPerCycle + r) * numOfProcs + rank;
                    if (row < numOfRows) {
                        L.getRowVector(row, &(v[0]), numOfCols);
                        std::copy(v.begin(),
                                  v.begin() + numOfCols,
                                  buf.begin() + numOfCols * r);
                    }
                }
            }
            rComm.broadcast(&(buf[0]), transRowsPerCycle * numOfCols, proc);
            
            // set
            for (index_type r = 0; r < transRowsPerCycle; ++r) {
                index_type row = (cycle * transRowsPerCycle + r) * numOfProcs + proc;
                if (row < numOfRows) {
                    for (index_type col = 0; col < numOfCols; ++col) {
                        colVecL.set(row, col, buf[numOfCols * r + col]);
                    }
                }
            }
        }
    }
    
    colVecL.save(path);

    // if (this->isDebugSaveL_ == true) {
    //     {
    //         TlMatrix tmpL = L.getTlMatrix();
    //         rComm.allReduce_SUM(tmpL);
    //         if (rComm.isMaster()) {
    //             tmpL.save("L.mat");
    //         }
    //     }

    //     {
    //         TlMatrix tmpL = colVecL.getTlMatrix();
    //         rComm.allReduce_SUM(tmpL);
    //         if (rComm.isMaster()) {
    //             tmpL.save("L2.mat");
    //         }
    //     }
    // }
}


TlMatrix DfCD_Parallel::mergeL(const TlRowVectorMatrix2& L)
{
    this->log_.info("merge L(row): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfRows, numOfProcs);
    const index_type localRows = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type r = 0; r < localRows; ++r) {
        const index_type row = r * numOfProcs + rank;
        TlVector rowVec = L.getRowVector(row);
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, rowVec[col]);
        }
    }
    rComm.allReduce_SUM(answer);
    
    this->log_.info("merge L: end");
    return answer;
}


TlMatrix DfCD_Parallel::mergeL(const TlColVectorMatrix2& L)
{
    this->log_.info("merge L(col): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfCols, numOfProcs);
    const index_type localCols = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type c = 0; c < localCols; ++c) {
        const index_type col = c * numOfProcs + rank;
        TlVector colVec = L.getColVector(col);
        for (index_type row = 0; row < numOfRows; ++row) {
            answer.set(row, col, colVec[row]);
        }
    }
    rComm.allReduce_SUM(answer);
    
    this->log_.info("merge L: end");
    return answer;
}


void DfCD_Parallel::getJ(TlSymmetricMatrix* pJ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc J by CD method (parallel).");
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath());
    assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();

    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);

            TlSymmetricMatrix LI = 
                this->getCholeskyVector(cv, I2PQ);
            assert(LI.getNumOfRows() == this->m_nNumOfAOs);
            assert(LI.getNumOfCols() == this->m_nNumOfAOs);

            TlSymmetricMatrix QI = LI;
            QI.dot(P);
            const double qi = QI.sum();
            
            *pJ += qi*LI;
        }
    }

    rComm.allReduce_SUM(*pJ);
}


void DfCD_Parallel::getK(const RUN_TYPE runType,
                         TlSymmetricMatrix* pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K by CD method (parallel).");

    // cholesky vector
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath());
    assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P = 0.5 * this->getPMatrix(); // RKS
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);

            TlSymmetricMatrix l = 
                this->getCholeskyVector(cv, I2PQ);
            
            TlMatrix X = l * C;
            TlMatrix Xt = X;
            Xt.transpose();
            
            TlSymmetricMatrix XX = X * Xt;
            *pK += XX;
        }
    }
    rComm.allReduce_SUM(*pK);
    
    *pK *= -1.0;
}


void DfCD_Parallel::getJ_D(TlDistributeSymmetricMatrix* pJ)
{
    this->log_.info("calc J by CD method (parallel; distributed).");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlDistributeSymmetricMatrix P = 
        DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

    // cholesky vector
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank(), 
                         isUsingMemManager);
    L.load(DfObject::getLjkMatrixPath());
    assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    assert(L.getRank() == rComm.getRank());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();

    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);

        TlDistributeSymmetricMatrix LI = 
            this->getCholeskyVector_distribute(cv, I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlDistributeSymmetricMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();
        
        *pJ += qi*LI;
    }
}


void DfCD_Parallel::getK_D(const RUN_TYPE runType,
                           TlDistributeSymmetricMatrix* pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K by CD method (parallel; distributed).");


    // cholesky vector
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                         isUsingMemManager);
    L.load(DfObject::getLjkMatrixPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();

    this->log_.info("calc CD of density matrix");
    this->log_.info(TlUtils::format("epsilon = %8.3e", this->epsilon_));
    TlDistributeSymmetricMatrix P = 
        0.5 * DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1); // RKS
    //const TlDistributeMatrix C = P.choleskyFactorization(this->epsilon_);
    const TlDistributeMatrix C = P.choleskyFactorization_mod2(this->epsilon_);

    TlTime time_all;
    TlTime time_bcast;
    TlTime time_translate;
    TlTime time_multimat1;
    TlTime time_multimat2;
    TlTime time_sym2gen;
    TlTime time_transpose;
    TlTime time_add;
    
    time_all.start();
    this->log_.info("start loop");
    int progress = 0;
    const int division = numOfCVs * 0.1;
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        // progress
        if (I >= progress * division) {
            const double rate = double(I) / double(numOfCVs) * 100.0;
            this->log_.info(TlUtils::format("K loop progress: %5.2f%%", rate));
            ++progress;
        }

        time_bcast.start();
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);
        time_bcast.stop();

        time_translate.start();
        TlDistributeSymmetricMatrix l = 
            this->getCholeskyVector_distribute(cv, I2PQ);
        time_translate.stop();

        time_multimat1.start();
        TlDistributeMatrix X = l * C;
        time_multimat1.stop();

        time_sym2gen.start();
        TlDistributeMatrix Xt = X;
        time_sym2gen.stop();

        time_transpose.start();
        Xt.transpose();
        time_transpose.stop();
        
        time_multimat2.start();
        TlDistributeSymmetricMatrix XX = X * Xt;
        time_multimat2.stop();

        time_add.start();
        *pK += XX;
        time_add.stop();
    }
    
    *pK *= -1.0;
    time_all.stop();

    // timing data
    this->log_.info(TlUtils::format("K all:       %10.1f sec.", time_all.getElapseTime()));
    this->log_.info(TlUtils::format("K bcast:     %10.1f sec.", time_bcast.getElapseTime()));
    this->log_.info(TlUtils::format("K translate: %10.1f sec.", time_translate.getElapseTime()));
    this->log_.info(TlUtils::format("K multi1:    %10.1f sec.", time_multimat1.getElapseTime()));
    this->log_.info(TlUtils::format("K multi2:    %10.1f sec.", time_multimat2.getElapseTime()));
    this->log_.info(TlUtils::format("K sym2gen:   %10.1f sec.", time_sym2gen.getElapseTime()));
    this->log_.info(TlUtils::format("K transpose: %10.1f sec.", time_transpose.getElapseTime()));
    this->log_.info(TlUtils::format("K add:       %10.1f sec.", time_add.getElapseTime()));
}

void DfCD_Parallel::getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    DfCD::getM(P, pM);
}

void DfCD_Parallel::getM_S(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (symmetric routine; parallel)");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets"]);
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pM->resize(numOfAOs);

    // cholesky vector
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath());
    assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();

    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);

            TlSymmetricMatrix LI = 
                this->getCholeskyVector(cv, I2PQ);
            assert(LI.getNumOfRows() == this->m_nNumOfAOs);
            assert(LI.getNumOfCols() == this->m_nNumOfAOs);

            TlSymmetricMatrix QI = LI;
            QI.dot(P);
            const double qi = QI.sum();
            
            *pM += qi*LI;
        }
    }

    rComm.allReduce_SUM(*pM);
}

void DfCD_Parallel::getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (asymmetric routine; parallel)");

    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets_GF"]);
    // const index_type numOfAOs = orbInfo_p.getNumOfOrbitals();
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath());
    assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();

    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);

            TlSymmetricMatrix l = 
                this->getCholeskyVector(cv, I2PQ);
            
            TlMatrix X = l * C;
            TlMatrix Xt = X;
            Xt.transpose();
            
            TlSymmetricMatrix XX = X * Xt;
            *pM += XX;
        }
    }
    rComm.allReduce_SUM(*pM);
}

void DfCD_Parallel::getM(const TlDistributeSymmetricMatrix& P,
                         TlDistributeSymmetricMatrix* pM)
{
    this->log_.critical("sorry, NO IMPLEMENTED!");
    this->log_.critical("DfCD_Parallel::getM(const TlDistributeSymmetricMatrix&, TlDistributeSymmetricMatrix*)");
    abort();
    // if (this->isDedicatedBasisForGridFree_) {
    //     this->getM_A(P, pM);
    // } else {
    //     this->getM_S(P, pM);
    // }
}

void DfCD_Parallel::getM_S(const TlDistributeSymmetricMatrix& P,
                           TlDistributeSymmetricMatrix* pM)
{
    this->log_.critical("sorry, NO IMPLEMENTED!");
    this->log_.critical("DfCD_Parallel::getM_S(const TlDistributeSymmetricMatrix&, TlDistributeSymmetricMatrix*)");
    abort();
}

void DfCD_Parallel::getM_A(const TlDistributeSymmetricMatrix& P,
                           TlDistributeSymmetricMatrix* pM)
{
    this->log_.critical("sorry, NO IMPLEMENTED!");
    this->log_.critical("DfCD_Parallel::getM_A(const TlDistributeSymmetricMatrix&, TlDistributeSymmetricMatrix*)");
    abort();
}

