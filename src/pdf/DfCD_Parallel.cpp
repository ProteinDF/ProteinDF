#include <algorithm>
#include "DfCD_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"

DfCD_Parallel::DfCD_Parallel(TlSerializeData* pPdfParam) 
    : DfCD(pPdfParam) {
}


DfCD_Parallel::~DfCD_Parallel()
{
}


void DfCD_Parallel::makeSuperMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->makeSuperMatrix_distribute();
    } else {
        DfCD::makeSuperMatrix_screening();
    }
#else
    {
        DfCD::makeSuperMatrix_screening();
    }
#endif // HAVE_SCALAPACK
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
                i2pq_tmp[i] = PQ_Pair(shellArray[i*2   ],
                                      shellArray[i*2 +1]);
            }
            
            pI2PQ->insert(pI2PQ->end(),
                          i2pq_tmp.begin(), i2pq_tmp.end());
        }
        std::sort(pI2PQ->begin(), pI2PQ->end(), PQ_Pair_less());
    } else {
        const std::size_t I2PQ_size = pI2PQ->size();
        std::vector<index_type> shellArray(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i*2   ] = (*pI2PQ)[i].shellIndex1;
            shellArray[i*2 +1] = (*pI2PQ)[i].shellIndex2;
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
                tmp[i*2   ] = (*pI2PQ)[i].shellIndex1;
                tmp[i*2 +1] = (*pI2PQ)[i].shellIndex2;
            }
        }
        rComm.broadcast(tmp);
        if (rComm.isMaster() != true) {
            const std::size_t size = tmp.size() / 2;
            pI2PQ->resize(size);
            for (std::size_t i = 0; i < size; ++i) {
                (*pI2PQ)[i] = PQ_Pair(tmp[i*2   ],
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
            shellArray[i*2   ] = I2PQ[i].shellIndex1;
            shellArray[i*2 +1] = I2PQ[i].shellIndex2;
        }
    }
    
    rComm.broadcast(shellArray);
    if (rComm.isMaster() != true) {
        const std::size_t I2PQ_size = shellArray.size() / 2;
        I2PQ.resize(I2PQ_size);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            I2PQ[i] = PQ_Pair(shellArray[i*2   ],
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
        answer.set(I2PQ[i].shellIndex1,
                   I2PQ[i].shellIndex2,
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
        const size_type PQ2I_index = this->pqPairIndex(I2PQ[i]);
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
                                          this->grainSize_,
                                          &taskList,
                                          true);
    while (hasTask == true) {
        this->makeSuperMatrix_kernel2(orbitalInfo,
                                      taskList,
                                      PQ2I,
                                      &tmpG);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_,
                                         &taskList);
    }

    // finalize
    G.mergeSparseMatrix(tmpG);
    //G.save("G.mat");
    //std::cerr << TlUtils::format("G(%d, %d)", G.getNumOfRows(), G.getNumOfCols()) << std::endl;

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
