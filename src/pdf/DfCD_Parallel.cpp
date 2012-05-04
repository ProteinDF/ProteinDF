#include <algorithm>
#include "DfCD_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"

DfCD_Parallel::DfCD_Parallel(TlSerializeData* pPdfParam) 
    : DfCD(pPdfParam) {
    this->CD_all_time_.stop();
    this->CD_ERI_time_.stop();
    this->CD_calc_time_.stop();
    this->CD_calc_d1_time_.stop();
    this->CD_calc_d2_time_.stop();
}


DfCD_Parallel::~DfCD_Parallel()
{
    this->log_.info(TlUtils::format("CD all:    %f sec.", this->CD_all_time_.getElapseTime()));
    this->log_.info(TlUtils::format("CD ERI:    %f sec.", this->CD_ERI_time_.getElapseTime()));
    this->log_.info(TlUtils::format("CD calc:   %f sec.", this->CD_calc_time_.getElapseTime()));
    this->log_.info(TlUtils::format("CD d1:     %f sec.", this->CD_calc_d1_time_.getElapseTime()));
    this->log_.info(TlUtils::format("CD d2:     %f sec.", this->CD_calc_d2_time_.getElapseTime()));
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
    
    TlDistributeSymmetricMatrix P = DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(runType,
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
    pDfTaskCtrl->cutoffReport();

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
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->createEngines();
    this->initializeCutoffStats();

    this->CD_all_time_.start();

    this->log_.info(TlUtils::format("# of PQ dimension: %d", int(this->numOfPQs_)));
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector d; // 対角成分

    this->calcDiagonals(&schwartzTable, &I2PQ, &d);

    this->log_.info(TlUtils::format("# of I~ dimension: %d", int(I2PQ.size())));
    this->saveI2PQ(I2PQ);

    const index_type N = I2PQ.size();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
#pragma omp parallel for 
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    // prepare variables
    RowVectorMatrix L; // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));

    index_type m = 0;
    while (error > threshold) {
        // this->log_.info(TlUtils::format("m=%d: err=%f", m, error));
        L.resize(N, m+1);

        // pivot
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }
        
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm); // 通信発生せず。関係無いPEは値を捨てる。
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        this->CD_ERI_time_.start();
        const index_type pivot_m = pivot[m];
        std::vector<double> G_pm;
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        this->CD_ERI_time_.stop();
        assert(G_pm.size() == numOf_G_cols);

        // CD calc
        this->CD_calc_time_.start();
        TlVector L_pm;
        // L_pm = L.getRowVector(pivot_m);
        {
            // 全PEに分配
            const int PEinCharge = L.getPEinChargeByRow(pivot_m);
            if (PEinCharge == rComm.getRank()) {
                L_pm = L.getRowVector(pivot_m);
            }
            rComm.broadcast(L_pm, PEinCharge);
        }
        assert(L_pm.getSize() == (m+1));
        this->CD_calc_time_.stop();

        this->CD_calc_d1_time_.start();
        std::vector<double> tmp_d(numOf_G_cols);
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N

            if (L.getPEinChargeByRow(pivot_i) == rComm.getRank()) { // 自分がL(pivot_i, *)を持っていたら
                const TlVector L_pi = L.getRowVector(pivot_i);
                double sum_ll = 0.0;
                for (index_type j = 0; j < m; ++j) {
                    sum_ll += L_pm[j] * L_pi[j];
                }

                const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;
                L.set(pivot_i, m, l_m_pi);
            
                tmp_d[i] -= l_m_pi * l_m_pi;
                // d[pivot_i] -= l_m_pi * l_m_pi;
            }
        }
        this->CD_calc_d1_time_.start();

        this->CD_calc_d2_time_.start();
        rComm.allReduce_SUM(&(tmp_d[0]), numOf_G_cols);
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            d[pivot_i] += tmp_d[i];
        }
        this->CD_calc_d2_time_.stop();

        // calc error
        error = 0.0;
#pragma omp parallel for reduction(+: error)
        for (index_type i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }

        ++m;
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->destroyEngines();
    this->schwartzCutoffReport();

    this->CD_all_time_.stop();

    this->saveL(L);
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


void DfCD_Parallel::saveL(const RowVectorMatrix& L)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (this->m_bUsingSCALAPACK == true) {
        const TlDistributeMatrix tmpL = L.getTlDistributeMatrix();
        DfObject::saveLMatrix(tmpL);
    } else {
        const TlMatrix tmpL = L.getTlMatrix();
        if (rComm.isMaster() == true) {
            DfCD::saveL(tmpL);
        }
    }
}


// =============================================================================
DfCD_Parallel::RowVectorMatrix::RowVectorMatrix(const index_type row,
                                                const index_type col)
{
    this->globalRows_ = 0;
    this->globalCols_ = 0;
    this->resize(row, col);
}


DfCD_Parallel::RowVectorMatrix::~RowVectorMatrix()
{
}


void DfCD_Parallel::RowVectorMatrix::resize(const index_type new_globalRows,
                                            const index_type new_globalCols)
{
    const index_type prev_globalRows = this->globalRows_;

    if (prev_globalRows != new_globalRows) {
        // remake PE table
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int numOfProcs = rComm.getNumOfProcs();
        const int myRank = rComm.getRank();

        this->row_PE_table_.resize(new_globalRows);
        int pe = 0;
        for (index_type r = 0; r < new_globalRows; ++r) {
            this->row_PE_table_[r] = pe;
            
            ++pe;
            if (pe >= numOfProcs) {
                pe = 0;
            }
        }

        std::vector<index_type> localRowTable;
        localRowTable.reserve(new_globalRows / numOfProcs +1);
        for (index_type r = 0; r < new_globalRows; ++r) {
            if (this->row_PE_table_[r] == myRank) {
                localRowTable.push_back(r);
            }
        }

        const std::size_t prev_localRowTableSize = this->data_.size();
        const std::size_t new_localRowTableSize = localRowTable.size();
        if (prev_localRowTableSize < new_localRowTableSize) {
            this->data_.resize(new_localRowTableSize);
            for (std::size_t i = prev_localRowTableSize; i < new_localRowTableSize; ++i) {
                this->data_[i].row = localRowTable[i];
            }
        } else if (prev_localRowTableSize < new_localRowTableSize) {
            this->data_.resize(new_localRowTableSize);
        }
    }

    std::vector<RowVector>::iterator itEnd = this->data_.end();
    for (std::vector<RowVector>::iterator it = this->data_.begin(); it != itEnd; ++it) {
        it->cols.resize(new_globalCols);
    }

    this->globalRows_ = new_globalRows;
    this->globalCols_ = new_globalCols;
}


void DfCD_Parallel::RowVectorMatrix::set(index_type row, index_type col, double value)
{
    std::vector<RowVector>::iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), RowVector(row));

    if ((it != this->data_.end()) && (it->row == row)) {
        it->cols.set(col, value);
    }
}


TlVector DfCD_Parallel::RowVectorMatrix::getRowVector(index_type row) const
{
    TlVector answer;
    std::vector<RowVector>::const_iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), RowVector(row));
    if (it != this->data_.end()) {
        assert(it->row == row);
        answer = it->cols;
    }

    return answer;
}


int DfCD_Parallel::RowVectorMatrix::getPEinChargeByRow(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    return this->row_PE_table_[row];
}


TlMatrix DfCD_Parallel::RowVectorMatrix::getTlMatrix() const 
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlMatrix answer(row, col);
    
    std::vector<RowVector>::const_iterator itEnd = this->data_.end();
    for (std::vector<RowVector>::const_iterator it = this->data_.begin(); it != itEnd; ++it) {
        const index_type r = it->row;
        for (index_type c = 0; c < col; ++c) {
            answer.set(r, c, it->cols.get(c));
        }
    }

    rComm.allReduce_SUM(answer);
    return answer;
}


TlDistributeMatrix DfCD_Parallel::RowVectorMatrix::getTlDistributeMatrix() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlDistributeMatrix answer(row, col);

    int blocks = 0;
    std::vector<index_type> rows;
    std::vector<double> vtr;
    for (int pe = 0; pe < numOfProcs; ++pe) {
        if (pe == myRank) {
            blocks = this->data_.size();
            vtr.resize(col * blocks);
            for (int i = 0; i < blocks; ++i) {
                rows[i] = this->data_[i].row;
                for (index_type c = 0; c < col; ++c) {
                    vtr[col * i + c] = this->data_[i].cols.get(c);
                }
            }
        }

        rComm.broadcast(blocks, pe);
        rComm.broadcast(rows, pe);
        rComm.broadcast(vtr, pe);

        // set
        for (int i = 0; i < blocks; ++i) {
            const index_type r = rows[i];
            for (index_type c = 0; c < col; ++c) {
                answer.set(r, c, vtr[col * i + c]);
            }
        }
    }

    return answer;
}


