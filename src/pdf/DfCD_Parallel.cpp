#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfCD_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
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


void DfCD_Parallel::finalize_I2PQ(I2PQ_Type* pI2PQ)
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
            I2PQ_Type i2pq_tmp(I2PQ_size);
            for (std::size_t i = 0; i < I2PQ_size; ++i) {
                i2pq_tmp[i] = IndexPair2(shellArray[i*2   ],
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
                (*pI2PQ)[i] = IndexPair2(tmp[i*2   ],
                                         tmp[i*2 +1]);
            }
        }
    }
}


void DfCD_Parallel::saveI2PQ(const I2PQ_Type& I2PQ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCD::saveI2PQ(I2PQ);
    }
}


DfCD::I2PQ_Type DfCD_Parallel::getI2PQ()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("distribute I2PQ table.");

    I2PQ_Type I2PQ;
    std::vector<index_type> shellArray;
    if (rComm.isMaster() == true) {
        I2PQ = DfCD::getI2PQ();

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
            I2PQ[i] = IndexPair2(shellArray[i*2   ],
                                 shellArray[i*2 +1]);
        }
    }

    return I2PQ;
}

void DfCD_Parallel::saveL(const TlMatrix& L)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCD::saveL(L);
    }
}


TlMatrix DfCD_Parallel::getL()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlMatrix L;
    if (rComm.isMaster() == true) {
        L = DfCD::getL();
    }
    rComm.broadcast(L);

    return L;
}


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


void DfCD_Parallel::getJ_distributed(TlDistributeSymmetricMatrix* pJ)
{
    this->log_.info("calc J by CD method on distributed matrix.");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlDistributeSymmetricMatrix P = 
        DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

    // cholesky vector
    TlDistributeMatrix L = DfObject::getLMatrix<TlDistributeMatrix>();
    const index_type numOfCBs = L.getNumOfCols();

    const I2PQ_Type I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = numOfCBs;
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlDistributeSymmetricMatrix LI = 
            this->getCholeskyVector_distribute(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);

        TlDistributeSymmetricMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();
        
        *pJ += qi*LI;
    }

    //this->finalize(pJ);
}


void DfCD_Parallel::getK_distributed(const RUN_TYPE runType,
                                     TlDistributeSymmetricMatrix* pK)
{
    this->log_.info("calc K by CD method on distributed matrix.");

    TlDistributeMatrix L = DfObject::getLMatrix<TlDistributeMatrix>();
    const index_type numOfCBs = L.getNumOfCols();
    
    // for RKS
    TlDistributeSymmetricMatrix P = 0.5 * DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(runType,
                                                                                              this->m_nIteration -1);
    const TlDistributeMatrix C = P.choleskyFactorization(this->epsilon_);
    
    const I2PQ_Type I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = numOfCBs;
    //this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlDistributeSymmetricMatrix l = this->getCholeskyVector_distribute(L.getColVector(I), I2PQ);

        TlDistributeMatrix X = l * C;
        TlDistributeMatrix Xt = X;
        Xt.transpose();
        
        TlDistributeSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
    //this->finalize(pK);
}


TlDistributeSymmetricMatrix 
DfCD_Parallel::getCholeskyVector_distribute(const TlVector& L_col,
                                            const I2PQ_Type& I2PQ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
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


void DfCD_Parallel::makeSuperMatrix_distribute()
{
    const index_type numOfPQs = this->numOfPQs_;

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    // calc (pq|pq)
    TlSparseSymmetricMatrix schwarzTable;
    PQ_PairArray I2PQ; // I~ to (pq) index table; size of (I2PQ) is the number of I~.
    this->calcPQPQ(orbitalInfo, &schwarzTable, &I2PQ);
    this->saveI2PQ(I2PQ);
    const index_type numOfItilde = I2PQ.size();
    this->log_.info(TlUtils::format(" # of PQ dimension: %d", int(numOfPQs)));
    this->log_.info(TlUtils::format(" # of I~ dimension: %d", int(numOfItilde)));

    // make PQ2I from I2PQ
    PQ2I_Type PQ2I(numOfPQs, -1);
    for (size_type i = 0; i < numOfItilde; ++i) {
        const size_type PQ2I_index = I2PQ[i].index();
        assert(PQ2I_index < numOfPQs);
        PQ2I[PQ2I_index] = i;
    }

    // 
    TlDistributeSymmetricMatrix G = this->getGMatrix_distribute(orbitalInfo, schwarzTable, numOfItilde, PQ2I);

    this->makeL(G);
}


TlDistributeSymmetricMatrix 
DfCD_Parallel::getGMatrix_distribute(const TlOrbitalInfoObject& orbitalInfo, 
                                     const TlSparseSymmetricMatrix& schwarzTable,
                                     const index_type numOfItilde,
                                     const PQ2I_Type& PQ2I)
{
    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    TlDistributeSymmetricMatrix G(numOfItilde);
    TlSparseSymmetricMatrix tmpG(numOfItilde);
    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->makeSuperMatrix_kernel2(orbitalInfo,
                                      taskList,
                                      PQ2I,
                                      &tmpG);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_, &taskList);
    }

    // finalize
    G.mergeSparseMatrix(tmpG);
    //G.save("G.mat");
    //std::cerr << TlUtils::format("G(%d, %d)", G.getNumOfRows(), G.getNumOfCols()) << std::endl;
    //pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

    return G;
}


void DfCD_Parallel::makeL(const TlDistributeSymmetricMatrix& G)
{
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    TlDistributeMatrix L = G.choleskyFactorization_mod(this->epsilon_);
    //std::cerr << TlUtils::format("L(%d, %d)", L.getNumOfRows(), L.getNumOfCols()) << std::endl;
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", L.getNumOfCols()));

    DfObject::saveLMatrix(L);
}

// On the Fly method -----------------------------------------------------------
void DfCD_Parallel::calcCholeskyVectors_onTheFly()
{
    // timing data
    TlTime CD_all_time;
    TlTime CD_diagonals_time;
    TlTime CD_resizeL_time;
    TlTime CD_ERI_time;
    TlTime CD_Lpm_time;
    TlTime CD_calc_time;
    TlTime CD_d_time;
    TlTime CD_save_time;

    CD_all_time.start();
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->createEngines();
    this->initializeCutoffStats();

    CD_diagonals_time.start();
    this->log_.info(TlUtils::format("# of PQ dimension: %d", int(this->numOfPQs_)));
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector global_diagonals; // 対角成分
    this->calcDiagonals(&schwartzTable, &I2PQ, &global_diagonals);
    this->log_.info(TlUtils::format("# of I~ dimension: %d", int(I2PQ.size())));
    this->saveI2PQ(I2PQ);
    CD_diagonals_time.stop();

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();
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
    index_type division = index_type(N * 0.01);
    while (error > threshold) {
#ifdef CD_DEBUG
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif // CD_DEBUG

        // progress 
        CD_resizeL_time.start();
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e, ERI cache=%ld MB",
                                            m, N, error,
                                            this->eriCache_.size() * (sizeof(IndexPair4) 
                                                                      + sizeof(double)) / (1024*1024)));
            ++progress;

            // メモリの確保
            L.reserve_cols(progress * division);
        }
        L.resize(N, m+1);
        CD_resizeL_time.stop();

        // pivot
        std::swap(global_pivot[m], global_pivot[error_global_loc]);
        reverse_pivot[global_pivot[m]] = m;
        reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

        const double l_m_pm = std::sqrt(error);
        const index_type pivot_m = global_pivot[m];
        L.set(pivot_m, m, l_m_pm); // 通信発生せず。関係無いPEは値を捨てる。
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        CD_ERI_time.start();
        std::vector<double> G_pm;
        // const index_type numOf_G_cols = N -(m+1);
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = global_pivot[m+1 +i]; // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(G_pm.size() == numOf_G_cols);
        CD_ERI_time.stop();

        // CD calc
        CD_Lpm_time.start();
        {
            // 全PEに分配
            const int PEinCharge = L.getPEinChargeByRow(pivot_m);
            if (PEinCharge == rComm.getRank()) {
                const index_type copySize = L.getRowVector(pivot_m, &(L_pm[0]), m +1);
                assert(copySize == m +1);
            }
            rComm.broadcast(&(L_pm[0]), m +1, PEinCharge);
        }
        CD_Lpm_time.stop();

        CD_calc_time.start();
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
        CD_calc_time.stop();

        CD_d_time.start();
        rComm.allReduce_MAXLOC(&error, &error_global_loc);
        global_diagonals[global_pivot[error_global_loc]] = error;
        CD_d_time.stop();

        ++m;
        if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
            std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
            ++local_m;
        }
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->destroyEngines();
    this->schwartzCutoffReport();

    CD_save_time.start();
    this->saveL(L);
    CD_save_time.stop();

    CD_all_time.stop();

    // timing data
    this->log_.info(TlUtils::format("CD all:       %f sec.", CD_all_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD diagonals: %f sec.", CD_diagonals_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD resize L:  %f sec.", CD_resizeL_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD ERI:       %f sec.", CD_ERI_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD L(m):      %f sec.", CD_Lpm_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD calc:      %f sec.", CD_calc_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD d:         %f sec.", CD_d_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD save:      %f sec.", CD_save_time.getElapseTime()));
}


std::vector<double>
DfCD_Parallel::getSuperMatrixElements(const index_type G_row,
                                      const std::vector<index_type>& G_col_list,
                                      const I2PQ_Type& I2PQ,
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
        elements = DfCD::getSuperMatrixElements(G_row,
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


void DfCD_Parallel::saveL(const TlRowVectorMatrix2& L)
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
    
    colVecL.save(DfObject::getLMatrixPath());

    if (this->isDebugSaveL_ == true) {
        {
            TlMatrix tmpL = L.getTlMatrix();
            rComm.allReduce_SUM(tmpL);
            if (rComm.isMaster()) {
                tmpL.save("L.mat");
            }
        }

        {
            TlMatrix tmpL = colVecL.getTlMatrix();
            rComm.allReduce_SUM(tmpL);
            if (rComm.isMaster()) {
                tmpL.save("L2.mat");
            }
        }
    }
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
    this->log_.info("calc J by CD method (parallel).");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    const I2PQ_Type I2PQ = this->getI2PQ();
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLMatrixPath());
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
    const I2PQ_Type I2PQ = this->getI2PQ();
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLMatrixPath());

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
    const I2PQ_Type I2PQ = this->getI2PQ();
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank(), 
                         isUsingMemManager);
    L.load(DfObject::getLMatrixPath());
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
    const I2PQ_Type I2PQ = this->getI2PQ();
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix2 L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                         isUsingMemManager);
    L.load(DfObject::getLMatrixPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCBs = L.getNumOfCols();
    
    TlDistributeSymmetricMatrix P = 
        0.5 * DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1); // RKS
    const TlDistributeMatrix C = P.choleskyFactorization(this->epsilon_);
    
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCBs; ++I) {
        const int PEinCharge = L.getPEinChargeByCol(I);
        if (PEinCharge == rComm.getRank()) {
            const index_type copySize = L.getColVector(I, &(cv[0]), cvSize);
            assert(copySize == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);
            
        TlDistributeSymmetricMatrix l = 
            this->getCholeskyVector_distribute(cv, I2PQ);
        
        TlDistributeMatrix X = l * C;
        TlDistributeMatrix Xt = X;
        Xt.transpose();
        
        TlDistributeSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
}


