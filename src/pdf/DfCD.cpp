#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <set>
#include "DfCD.h"
#include "DfEriEngine.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"


DfCD::DfCD(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEriEngines_(NULL),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_sets"])
{
    this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    

    this->CDAM_tau_ = 1.0E-5;
    if ((*pPdfParam)["CDAM_tau"].getStr().empty() != true) {
        this->CDAM_tau_ = (*pPdfParam)["CDAM_tau"].getDouble();
    }    

    this->epsilon_ = std::sqrt(this->cutoffThreshold_);
    if ((*pPdfParam)["CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["CD_epsilon"].getDouble();
    }    

    this->cutoffEpsilon3_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon3"].getStr().empty() != true) {
        this->cutoffEpsilon3_ = (*pPdfParam)["cutoff_epsilon3"].getDouble();
    }    

}

DfCD::~DfCD()
{
}

void DfCD::createEngines()
{
    assert(this->pEriEngines_ == NULL);
    
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->pEriEngines_ = new DfEriEngine[numOfThreads];
    }
#else
    this->pEriEngines_ = new DfEriEngine[1];
#endif // _OPENMP
}


void DfCD::destroyEngines()
{
    if (this->pEriEngines_ != NULL) {
        delete[] this->pEriEngines_;
    }
    this->pEriEngines_ = NULL;
}

void DfCD::calcCholeskyVectors()
{
    this->calcCholeskyVectors_onTheFly();
    //this->makeSuperMatrix_screening();
    //this->makeSuperMatrix_noScreening();
}

void DfCD::makeSuperMatrix_screening()
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
    TlSymmetricMatrix G = this->getGMatrix(orbitalInfo, schwarzTable, numOfItilde, PQ2I);

    this->makeL(G);
}


void DfCD::calcPQPQ(const TlOrbitalInfoObject& orbitalInfo,
                    TlSparseSymmetricMatrix *pSchwarzTable,
                    PQ_PairArray *pI2PQ)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(numOfAOs == orbitalInfo.getNumOfOrbitals());

    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format(" CDAM tau: %e", tau));

    this->createEngines();
    this->log_.info(TlUtils::format(" pGTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(this->numOfPQs_);
    pSchwarzTable->clear();
    pSchwarzTable->resize(numOfAOs);

    // task
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcPQPQ_kernel(orbitalInfo,
                              taskList,
                              &(*pSchwarzTable), &(*pI2PQ));
        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    
    // finalize
    pDfTaskCtrl->cutoffReport();
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

    this->finalize(pSchwarzTable);
    this->finalize_I2PQ(pI2PQ);
}


void DfCD::calcPQPQ_kernel(const TlOrbitalInfoObject& orbitalInfo,
                           const std::vector<DfTaskCtrl::Task2>& taskList,
                           TlSparseSymmetricMatrix *pSchwarzTable,
                           PQ_PairArray *pI2PQ)
{
    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseSymmetricMatrix local_schwarzTable(pSchwarzTable->getNumOfRows());
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = 
                this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                          shellIndexP,
                                                          shellIndexQ,
                                                          pairwisePGTO_cutoffThreshold);
            this->pEriEngines_[threadID].calc(queryPQ, queryPQ, PQ, PQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;

                    if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                        const int pq_index = p * maxStepsQ + q;
                        const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                        
                        const double value = std::fabs(this->pEriEngines_[threadID].WORK[pqpq_index]);
                        
                        // for schwartz
                        maxValue = std::max(maxValue, value);
                        
                        // for I~ to pq table
                        if (value > tau) {
                            local_I2PQ.push_back(PQ_Pair(indexP, indexQ));
                        }
                    }
                }
            }
            local_schwarzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
            const int numOfThreads = omp_get_num_threads();
            for (int i = 0; i < numOfThreads; ++i) {
                if (threadID == i) {
                    pI2PQ->insert(pI2PQ->end(),
                                  local_I2PQ.begin(), local_I2PQ.end());
                    pSchwarzTable->merge(local_schwarzTable);
                }
#pragma omp barrier                
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pSchwarzTable = localSchwarzTable;
        }
#endif // _OPENMP
    }
}


TlSymmetricMatrix DfCD::getGMatrix(const TlOrbitalInfoObject& orbitalInfo, 
                                   const TlSparseSymmetricMatrix& schwarzTable,
                                   const index_type numOfItilde,
                                   const PQ2I_Type& PQ2I)
{
    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    TlSymmetricMatrix G(numOfItilde);
#ifdef CHECK_LOOP
    this->check.resize(numOfItilde);
#endif // CHECK_LOOP
    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->makeSuperMatrix_kernel2(orbitalInfo,
                                      taskList,
                                      PQ2I,
                                      &G);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_, &taskList);
    }

    this->finalize(&G);
    //G.save("G.mat");
    //std::cerr << TlUtils::format("G(%d, %d)", G.getNumOfRows(), G.getNumOfCols()) << std::endl;
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

#ifdef CHECK_LOOP
    {
        this->check.save("check.mat");
        const index_type dim = this->check.getNumOfRows();
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                if (std::fabs(this->check.get(i, j) - 1.0) > 1.0E-5) {
                    std::cerr << TlUtils::format("count err: (%d, %d)=%f", i, j, this->check.get(i, j))
                              << std::endl;
                }
            }
        }
    }
#endif // CHECK_LOOP

    return G;
}


void DfCD::makeSuperMatrix_kernel2(const TlOrbitalInfoObject& orbitalInfo,
                                   const std::vector<DfTaskCtrl::Task4>& taskList,
                                   const PQ2I_Type& PQ2I,
                                   TlMatrixObject* pG)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const index_type shellIndexR = taskList[i].shellIndex3;
            const index_type shellIndexS = taskList[i].shellIndex4;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
            const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;

            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexR,
                                                                                        shellIndexS,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
            
            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
            this->storeG2(shellIndexP, maxStepsP,
                          shellIndexQ, maxStepsQ,
                          shellIndexR, maxStepsR,
                          shellIndexS, maxStepsS,
                          PQ2I,
                          this->pEriEngines_[threadID], pG);
        }
    }
}

void DfCD::storeG2(const index_type shellIndexP, const int maxStepsP,
                   const index_type shellIndexQ, const int maxStepsQ,
                   const index_type shellIndexR, const int maxStepsR,
                   const index_type shellIndexS, const int maxStepsS,
                   const PQ2I_Type& PQ2I,
                   const DfEriEngine& engine,
                   TlMatrixObject* pG)
{
    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;

        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            PQ_Pair pq(indexP, indexQ);

            const size_type indexPQ = PQ2I[this->pqPairIndex(pq)];
            // if (indexPQ == -1) {
            //     index += maxStepsR * maxStepsS;
            //     continue;
            // }

            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;
                const index_type maxIndexS = (indexP == indexR) ? indexQ : indexR;

                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    PQ_Pair rs(indexR, indexS);

                    const index_type indexRS = PQ2I[this->pqPairIndex(rs)];
                    // if (indexRS == -1) {
                    //     ++index;
                    //     continue;
                    // }
                    
                    const double value = engine.WORK[index];

                    if ((indexPQ != -1) && (indexRS != -1)) {
                        if (indexQ <= indexP) {
                            if ((shellIndexQ != shellIndexS) || (indexR <= indexP)) {
                                if (indexS <= maxIndexS) {
                                    pG->set(indexPQ, indexRS, value);
#ifdef CHECK_LOOP
                                    this->check.add(indexPQ, indexRS, 1.0);
#endif // CHECK_LOOP
                                }
                            }
                        }
                    }
                    ++index;
                 }
            }
        }
    }
}

void DfCD::saveI2PQ(const PQ_PairArray& I2PQ) 
{
    std::string filepath = "I2PQ.vtr";
    std::ofstream ofs;
    ofs.open(filepath.c_str(), std::ofstream::out | std::ofstream::binary);

    const std::size_t size = I2PQ.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
    for (std::size_t i = 0; i < size; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(I2PQ[i].shellIndex1)), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&(I2PQ[i].shellIndex2)), sizeof(index_type));
    }

    ofs.close();
}

DfCD::PQ_PairArray DfCD::getI2PQ()
{
    std::string filepath = "I2PQ.vtr";
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = PQ_Pair(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}


void DfCD::makeL(const TlSymmetricMatrix& G)
{
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    TlMatrix L = G.choleskyFactorization2omp(this->epsilon_);
    //std::cerr << TlUtils::format("L(%d, %d)", L.getNumOfRows(), L.getNumOfCols()) << std::endl;
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", L.getNumOfCols()));
    this->saveL(L);
}


void DfCD::saveL(const TlMatrix& L)
{
    DfObject::saveLMatrix(L);
}


TlMatrix DfCD::getL()
{
    TlMatrix L = DfObject::getLMatrix<TlMatrix>();
    return L;
}


std::size_t DfCD::index(index_type p, index_type q) const
{
    if (p < q) {
        std::swap(p, q);
    }

    // This class treats 'L' type matrix.
    // Follows means:
    //  index = row + (2 * this->m_nRows - (col +1)) * col / 2;
    std::size_t s = this->m_nNumOfAOs;
    s = s << 1; // means 's *= 2'

    std::size_t t = (s - (q +1)) * q;
    t = t >> 1; // means 't /= 2'

    return (p + t);
}


DfTaskCtrl* DfCD::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    // pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    // pDfTaskCtrl->setCutoffEpsilon_density(0.0);
    // pDfTaskCtrl->setCutoffEpsilon_distribution(this->CDAM_tau_);

    return pDfTaskCtrl;
}

void DfCD::finalize(TlSymmetricMatrix* pMat)
{
    // do nothing
}

void DfCD::finalize(TlSparseSymmetricMatrix *pMat) 
{
    // do nothing
}

void DfCD::finalize_I2PQ(PQ_PairArray *pI2PQ)
{
    std::sort(pI2PQ->begin(), pI2PQ->end());
}


TlSparseSymmetricMatrix DfCD::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo)
{
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    TlSparseSymmetricMatrix schwarz(maxShellIndex);

    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);
    
    for (index_type shellIndexP = 0; shellIndexP < maxShellIndex; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;

        for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                  shellIndexP,
                                                                  shellIndexQ,
                                                                  0.0);
            engine.calc(queryPQ, queryPQ, PQ, PQ);

            double maxValue = 0.0;
            const int maxIndex = maxStepsP * maxStepsQ;
            for (int index = 0; index < maxIndex; ++index) {
                maxValue = std::max(maxValue, std::fabs(engine.WORK[index]));
            }
            schwarz.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
            
            shellIndexQ += maxStepsQ;
        }
        shellIndexP += maxStepsP;
    }

    return schwarz;
}

TlSymmetricMatrix DfCD::getCholeskyVector(const TlVector& L_col,
                                          const I2PQ_Type& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].shellIndex1,
                   I2PQ[i].shellIndex2,
                   L_col[i]);
    }

    return answer;
}

void DfCD::getJ(TlSymmetricMatrix* pJ)
{
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    TlMatrix L = this->getL();
    const index_type numOfCBs = L.getNumOfCols();

    const I2PQ_Type I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();

        *pJ += qi*LI;
    }

    this->finalize(pJ);
}


void DfCD::divideCholeskyBasis(const index_type numOfCBs,
                               index_type *pStart, index_type *pEnd)
{
    *pStart = 0;
    *pEnd = numOfCBs;
}


void DfCD::getK(const RUN_TYPE runType,
                TlSymmetricMatrix *pK)
{
    TlMatrix L = this->getL();
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P = this->getPMatrix(); // RKS
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    const I2PQ_Type I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix l = this->getCholeskyVector(L.getColVector(I), I2PQ);
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
    this->finalize(pK);
}


TlSymmetricMatrix DfCD::getPMatrix()
{
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    return P;
}


/// @param numOfCDAMs [in] I~の総数
void DfCD::calcCholeskyVectors_onTheFly()
{
    this->createEngines();

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector d; // 対角成分
    this->calcDiagonals(&schwartzTable, &I2PQ, &d);
    this->log_.info(TlUtils::format(" # of PQ dimension: %d", int(this->numOfPQs_)));
    this->log_.info(TlUtils::format(" # of I~ dimension: %d", int(I2PQ.size())));
    this->saveI2PQ(I2PQ);

    const index_type N = I2PQ.size();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    // clear cutoff stats
    {
        const int maxShellType = this->orbitalInfo_.getMaxShellType();
        const int numOfShellPairType = maxShellType* maxShellType;
        const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;
        this->cutoffAll_schwartz_.clear();
        this->cutoffAlive_schwartz_.clear();
        this->cutoffAll_schwartz_.resize(numOfShellQuartetType, 0);
        this->cutoffAlive_schwartz_.resize(numOfShellQuartetType, 0);
    }

    // prepare variables
    TlMatrix L;
    index_type m = 0;
    double sum_ll = 0.0;
    TlSparseSymmetricMatrix request(N);

    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    while (error > threshold) {
        L.resize(N, m+1);
        std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                       pivot.end());
        const index_type i = it - pivot.begin();
        std::swap(pivot[m], pivot[i]);
        
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // request
        TlSparseSymmetricMatrix G(N);
        const index_type pivot_m = pivot[m];
        for (index_type i = m +1; i < N; ++i) {
            const index_type pivot_i = pivot[i];
            G.set(pivot_m, pivot_i, 0.0);
        }
        this->calcERIs(schwartzTable, I2PQ, &G);

        // calc
        for (index_type i = m +1; i < N; ++i) {
            const index_type pivot_i = pivot[i];
            double sum_ll = 0.0;
            for (index_type j = 0; j < m; ++j) {
                sum_ll += L.get(pivot_m, j) * L.get(pivot_i, j);
            }

            const double l_m_pi = (G.get(pivot_m, pivot_i) - sum_ll) * inv_l_m_pm;
            L.set(pivot_i, m, l_m_pi);
            
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
            
        error = 0.0;
        for (index_type i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }

        ++m;
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->destroyEngines();
    this->schwartzCutoffReport();

    this->saveL(L);
}


void DfCD::calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                         PQ_PairArray *pI2PQ,
                         TlVector *pDiagonals)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(numOfAOs == this->orbitalInfo_.getNumOfOrbitals());

    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format(" CDAM tau: %e", tau));

    //this->createEngines();
    this->log_.info(TlUtils::format(" pGTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    TlSparseSymmetricMatrix diagonalMat(numOfAOs);
    pI2PQ->clear();
    pI2PQ->reserve(this->numOfPQs_);
    pSchwartzTable->clear();
    pSchwartzTable->resize(numOfAOs);

    // task
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_kernel(taskList,
                                   pSchwartzTable,
                                   &diagonalMat, pI2PQ);
        hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    
    // finalize
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    //this->destroyEngines();
    this->finalize_I2PQ(pI2PQ);

    // set diagonals
    {
        const index_type numOfI = pI2PQ->size();
        pDiagonals->resize(numOfI);
        for (index_type i = 0; i < numOfI; ++i) {
            const index_type row = (*pI2PQ)[i].shellIndex1;
            const index_type col = (*pI2PQ)[i].shellIndex2;
            const double value = diagonalMat.get(row, col);
            (*pDiagonals)[i] = value;
        }
    }
}


void DfCD::calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                                TlSparseSymmetricMatrix *pSchwartzTable,
                                TlSparseSymmetricMatrix *pDiagonalMat,
                                PQ_PairArray *pI2PQ)
{
    const index_type numOfAOs = this->orbitalInfo_.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseSymmetricMatrix local_diagMat(numOfAOs);
        TlSparseSymmetricMatrix local_schwartzTable(numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = 
                this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                          shellIndexP,
                                                          shellIndexQ,
                                                          pairwisePGTO_cutoffThreshold);
            this->pEriEngines_[threadID].calc(queryPQ, queryPQ, PQ, PQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                        const int pq_index = p * maxStepsQ + q;
                        const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                        
                        const double value = this->pEriEngines_[threadID].WORK[pqpq_index];

                        // for schwartz
                        maxValue = std::max(maxValue, std::fabs(value));
                        
                        // for I~ to pq table
                        if (std::fabs(value) > tau) {
                            local_diagMat.set(indexP, indexQ, value);
                            local_I2PQ.push_back(PQ_Pair(indexP, indexQ));
                        }
                    }
                }
            }
            local_schwartzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
            const int numOfThreads = omp_get_num_threads();
            for (int i = 0; i < numOfThreads; ++i) {
                if (threadID == i) {
                    pI2PQ->insert(pI2PQ->end(),
                                  local_I2PQ.begin(), local_I2PQ.end());
                    pDiagonalMat->merge(local_diagMat);
                    pSchwartzTable->merge(local_schwartzTable);
                }
#pragma omp barrier                
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
            *pSchwartzTable = local_schwartzTable;
        }
#endif // _OPENMP
    }
}


void DfCD::calcERIs(const TlSparseSymmetricMatrix& schwartzTable,
                    const I2PQ_Type& I2PQ,
                    TlSparseSymmetricMatrix* pG)
{
    static const int basisTypeBase[] = {0, 1, 4}; // s, px, dxy
    const int maxShellType = this->orbitalInfo_.getMaxShellType();
    const double threshold = this->CDAM_tau_;

    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;
    int threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

    TlSparseSymmetricMatrix::const_iterator itEnd = pG->end();
    for (TlSparseSymmetricMatrix::const_iterator it = pG->begin(); it != itEnd; ++it) {
        index_type G_row = 0;
        index_type G_col = 0;
        pG->index(it->first, &G_row, &G_col);

        const index_type indexP = I2PQ[G_row].shellIndex1;
        const index_type indexQ = I2PQ[G_row].shellIndex2;
        const index_type indexR = I2PQ[G_col].shellIndex1;
        const index_type indexS = I2PQ[G_col].shellIndex2;

        const int basisTypeP = this->orbitalInfo_.getBasisType(indexP) - basisTypeBase[this->orbitalInfo_.getShellType(indexP)];
        const int basisTypeQ = this->orbitalInfo_.getBasisType(indexQ) - basisTypeBase[this->orbitalInfo_.getShellType(indexQ)];
        const int basisTypeR = this->orbitalInfo_.getBasisType(indexR) - basisTypeBase[this->orbitalInfo_.getShellType(indexR)];
        const int basisTypeS = this->orbitalInfo_.getBasisType(indexS) - basisTypeBase[this->orbitalInfo_.getShellType(indexS)];
        const index_type shellIndexP = indexP - basisTypeP;
        const index_type shellIndexQ = indexQ - basisTypeQ;
        const index_type shellIndexR = indexR - basisTypeR;
        const index_type shellIndexS = indexS - basisTypeS;
        assert((this->orbitalInfo_.getBasisType(shellIndexP) == 0) ||
               (this->orbitalInfo_.getBasisType(shellIndexP) == 1) ||
               (this->orbitalInfo_.getBasisType(shellIndexP) == 4));
        assert((this->orbitalInfo_.getBasisType(shellIndexQ) == 0) ||
               (this->orbitalInfo_.getBasisType(shellIndexQ) == 1) ||
               (this->orbitalInfo_.getBasisType(shellIndexQ) == 4));
        assert((this->orbitalInfo_.getBasisType(shellIndexR) == 0) ||
               (this->orbitalInfo_.getBasisType(shellIndexR) == 1) ||
               (this->orbitalInfo_.getBasisType(shellIndexR) == 4));
        assert((this->orbitalInfo_.getBasisType(shellIndexS) == 0) ||
               (this->orbitalInfo_.getBasisType(shellIndexS) == 1) ||
               (this->orbitalInfo_.getBasisType(shellIndexS) == 4));

        const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
        const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
        const int shellTypeR = this->orbitalInfo_.getShellType(shellIndexR);
        const int shellTypeS = this->orbitalInfo_.getShellType(shellIndexS);

        const int shellQuartetType =
            ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
        const bool isAlive = this->isAliveBySchwartzCutoff(shellIndexP, shellIndexQ,
                                                           shellIndexR, shellIndexS,
                                                           shellQuartetType,
                                                           schwartzTable,
                                                           threshold);
        if (isAlive == true) {
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                                                        shellIndexR,
                                                                                        shellIndexS,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
            
            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
            
            const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
            const double value = this->pEriEngines_[threadID].WORK[index];
            pG->set(G_row, G_col, value);
        }
    }
}


bool DfCD::isAliveBySchwartzCutoff(const index_type shellIndexP,
                                   const index_type shellIndexQ,
                                   const index_type shellIndexR,
                                   const index_type shellIndexS,
                                   const int shellQuartetType,
                                   const TlSparseSymmetricMatrix& schwarzTable,
                                   const double threshold)
{
    bool answer = false;

    const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
    const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

    if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
        answer = true;

#pragma omp atomic
        ++(this->cutoffAlive_schwartz_[shellQuartetType]);
    }

#pragma omp atomic
    ++(this->cutoffAll_schwartz_[shellQuartetType]);

    return answer;
}


void DfCD::schwartzCutoffReport()
{
    static const char typeStr4[][5] = {
        "SSSS", "SSSP", "SSSD", "SSPS", "SSPP", "SSPD", "SSDS", "SSDP", "SSDD",
        "SPSS", "SPSP", "SPSD", "SPPS", "SPPP", "SPPD", "SPDS", "SPDP", "SPDD",
        "SDSS", "SDSP", "SDSD", "SDPS", "SDPP", "SDPD", "SDDS", "SDDP", "SDDD",
        "PSSS", "PSSP", "PSSD", "PSPS", "PSPP", "PSPD", "PSDS", "PSDP", "PSDD",
        "PPSS", "PPSP", "PPSD", "PPPS", "PPPP", "PPPD", "PPDS", "PPDP", "PPDD",
        "PDSS", "PDSP", "PDSD", "PDPS", "PDPP", "PDPD", "PDDS", "PDDP", "PDDD",
        "DSSS", "DSSP", "DSSD", "DSPS", "DSPP", "DSPD", "DSDS", "DSDP", "DSDD",
        "DPSS", "DPSP", "DPSD", "DPPS", "DPPP", "DPPD", "DPDS", "DPDP", "DPDD",
        "DDSS", "DDSP", "DDSD", "DDPS", "DDPP", "DDPD", "DDDS", "DDDP", "DDDD",
    };
    const int maxShellType = this->orbitalInfo_.getMaxShellType();

    // cutoff report for schwarz
    bool hasCutoffSchwarz = false;
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
            for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                    const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                    if (this->cutoffAll_schwartz_[shellTypeABCD] != 0) {
                        hasCutoffSchwarz = true;
                        break;
                    }
                }
            }
        }
    }
    if (hasCutoffSchwarz == true) {
        this->log_.info("schwarz cutoff report");
        this->log_.info(TlUtils::format("threshold: %e", this->CDAM_tau_));
        this->log_.info("type: alive / all (ratio)");
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
                for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                    const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                    for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                        const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                        
                        if (this->cutoffAll_schwartz_[shellTypeABCD] > 0) {
                            const double ratio = (double)this->cutoffAlive_schwartz_[shellTypeABCD]
                                / (double)this->cutoffAll_schwartz_[shellTypeABCD]
                                * 100.0;
                            this->log_.info(TlUtils::format(" %4s: %12ld / %12ld (%6.2f%%)",
                                                            typeStr4[shellTypeABCD],
                                                            this->cutoffAlive_schwartz_[shellTypeABCD],
                                                            this->cutoffAll_schwartz_[shellTypeABCD],
                                                            ratio));
                        }
                    }
                }
            }
        }
    }
}
