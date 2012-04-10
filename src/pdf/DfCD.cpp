#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <set>
#include "DfCD.h"
#include "DfEriEngine.h"
#include "TlOrbitalInfo.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

const int DfCD::MAX_SHELL_TYPE = 2 + 1;


DfCD::DfCD(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEriEngines_(NULL)
{
    this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    

    this->epsilon_ = std::sqrt(this->cutoffThreshold_);
    if ((*pPdfParam)["CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["CD_epsilon"].getDouble();
    }    


    // this->cutoffEpsilon1_ = this->cutoffThreshold_ * 0.01;
    // if ((*pPdfParam)["cutoff_epsilon1"].getStr().empty() != true) {
    //     this->cutoffEpsilon1_ = (*pPdfParam)["cutoff_epsilon1"].getDouble();
    // }    

    // this->cutoffEpsilon2_ = this->cutoffThreshold_;
    // if ((*pPdfParam)["cutoff_epsilon2"].getStr().empty() != true) {
    //     this->cutoffEpsilon2_ = (*pPdfParam)["cutoff_epsilon2"].getDouble();
    // }    

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

void DfCD::makeSuperMatrix()
{
    this->makeSuperMatrix_screening();
    //this->makeSuperMatrix_noScreening();
}

void DfCD::makeSuperMatrix_screening()
{
    //const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfPQs = this->numOfPQs_;

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    
    TlSparseSymmetricMatrix schwarzTable;
    std::vector<index_type> I2PQ; // I~ to (pq) index table; size of (I2PQ) is the number of I~.
    this->calcPQPQ(orbitalInfo, &schwarzTable, &I2PQ);
    const index_type numOfItilde = I2PQ.size();
    this->log_.info(TlUtils::format(" # of PQ dimension: %d", int(numOfPQs)));
    this->log_.info(TlUtils::format(" # of I~ dimension: %d", int(numOfItilde)));

    // 
    std::vector<index_type> PQ2I(numOfPQs, -1);
    for (index_type i = 0; i < numOfItilde; ++i) {
        const index_type pq_index = I2PQ[i];
        assert(pq_index < numOfPQs);
        PQ2I[pq_index] = i;
    }

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    TlSymmetricMatrix G(numOfItilde);
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
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    TlMatrix L = G.choleskyFactorization2(this->epsilon_);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", L.getNumOfCols()));
    L.save("L.mat");

    this->saveI2PQ(I2PQ);
}


void DfCD::calcPQPQ(const TlOrbitalInfoObject& orbitalInfo,
                    TlSparseSymmetricMatrix *pSchwarzTable,
                    std::vector<index_type> *pI2PQ)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(numOfAOs == orbitalInfo.getNumOfOrbitals());

    const double threshold = this->cutoffThreshold_;
    this->log_.info(TlUtils::format(" I~ threshold: %e", threshold));

    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
    this->log_.info(TlUtils::format(" pGTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(this->numOfPQs_);
    pSchwarzTable->clear();
    pSchwarzTable->resize(numOfAOs);

    for (index_type shellIndexP = 0; shellIndexP < numOfAOs; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;

        for (index_type shellIndexQ = 0; shellIndexQ <= shellIndexP; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                  shellIndexP,
                                                                  shellIndexQ,
                                                                  0.0);
            engine.calc(queryPQ, queryPQ, PQ, PQ);

            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    const int pq_index = p * maxStepsQ + q;
                    const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                    
                    const double value = std::fabs(engine.WORK[pqpq_index]);

                    // for schwartz
                    maxValue = std::max(maxValue, value);

                    // for I~ to pq table
                    if (value > threshold) {
                        pI2PQ->push_back(TlSymmetricMatrix::vtr_index(indexP, indexQ, numOfAOs));
                    }
                }
            }
            pSchwarzTable->set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
            
            shellIndexQ += maxStepsQ;
        }
        shellIndexP += maxStepsP;
    }
}

void DfCD::makeSuperMatrix_kernel2(const TlOrbitalInfo& orbitalInfo,
                                   const std::vector<DfTaskCtrl::Task4>& taskList,
                                   const std::vector<index_type>& PQ2I,
                                   TlSymmetricMatrix* pG)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        //this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
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
                   const std::vector<index_type>& PQ2I,
                   const DfEriEngine& engine,
                   TlSymmetricMatrix* pG)
{
    const index_type numOfAOs = this->m_nNumOfAOs;

    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;

        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            const index_type I_pq_index = TlSymmetricMatrix::vtr_index(indexP, indexQ, numOfAOs);
            assert(I_pq_index < PQ2I.size());
            const index_type indexPQ = PQ2I[I_pq_index];
            if (indexPQ == -1) {
                index += maxStepsR * maxStepsS;
                continue;
            }
            
            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;
                
                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    
                    const double value = engine.WORK[index];
                    const index_type maxIndexS = (indexP == indexR) ? indexQ : indexR;

                    const index_type I_rs_index = TlSymmetricMatrix::vtr_index(indexR, indexS, numOfAOs);
                    assert(I_rs_index < PQ2I.size());
                    const index_type indexRS = PQ2I[I_rs_index];
                    if (indexRS == -1) {
                        ++index;
                        continue;
                    }

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
                    ++index;
                 }
            }
        }
    }
}

void DfCD::saveI2PQ(const std::vector<int>& I2PQ) 
{
    std::string filepath = "I2PQ.vtr";
    std::ofstream ofs;
    ofs.open(filepath.c_str(), std::ofstream::out | std::ofstream::binary);

    const std::size_t size = I2PQ.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
    for (std::size_t i = 0; i < size; ++i) {
        int tmp = I2PQ[i];
        ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(int));
    }

    ofs.close();
}

std::vector<int> DfCD::getI2PQ()
{
    std::string filepath = "I2PQ.vtr";
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    std::vector<int> answer(size);
    int tmp = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int));
        answer[i] = tmp;
    }

    ifs.close();
    return answer;
}


void DfCD::makeSuperMatrix_noScreening()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type dim = numOfAOs * (numOfAOs +1) / 2;
    TlSymmetricMatrix G(dim);

#ifdef CHECK_LOOP
    this->check.resize(dim);
#endif // CHECK_LOOP
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
#ifdef CHECK_LOOP
    pDfTaskCtrl->setCutoffThreshold(-1.0);
#endif // CHECK_LOOP

    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->makeSuperMatrix_kernel(orbitalInfo,
                                     taskList,
                                     &G);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_, &taskList);
    }

    this->finalize(&G);
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

    G.save("G.mat");
#ifdef CHECK_LOOP
    this->check.save("check.mat");
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (std::fabs(this->check.get(i, j) - 1.0) > 1.0E-5) {
                std::cerr << TlUtils::format("count err: (%d, %d)=%f", i, j, this->check.get(i, j))
                          << std::endl;
            }
        }
    }
#endif // CHECK_LOOP
    
    TlMatrix L = G.choleskyFactorization2(this->epsilon_);
    L.save("L.mat");
}

void DfCD::makeSuperMatrix_kernel(const TlOrbitalInfo& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task4>& taskList,
                                  TlSymmetricMatrix* pG)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        //this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
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
            this->storeG(shellIndexP, maxStepsP,
                         shellIndexQ, maxStepsQ,
                         shellIndexR, maxStepsR,
                         shellIndexS, maxStepsS,
                         this->pEriEngines_[threadID], pG);
        }
    }
}

void DfCD::storeG(const index_type shellIndexP, const int maxStepsP,
                  const index_type shellIndexQ, const int maxStepsQ,
                  const index_type shellIndexR, const int maxStepsR,
                  const index_type shellIndexS, const int maxStepsS,
                  const DfEriEngine& engine,
                  TlSymmetricMatrix* pG)
{
    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;

        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            const std::size_t indexPQ = this->index(indexP, indexQ);
            
            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;
                
                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    
                    const double value = engine.WORK[index];
                    const index_type maxIndexS = (indexP == indexR) ? indexQ : indexR;
                    
                    const std::size_t indexRS = this->index(indexR, indexS);

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
                    ++index;
                 }
            }
        }
    }
}

void DfCD::makeSuperMatrix_exact()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfDims = numOfAOs * (numOfAOs +1) / 2;
    TlSymmetricMatrix G(numOfDims);
    //std::set<index_type> calcd; // 計算済みかどうかを格納
    
    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);
    
    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const int maxStepsP = 2 * shellTypeP + 1;
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeQ = MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            ShellArray::const_iterator qItEnd = shellArrayQ.end();

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);

            for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = MAX_SHELL_TYPE -1; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
        
                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;

                            //const index_type shellIndexPQ = this->index(shellIndexP, shellIndexQ);
                            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                                  shellIndexP, shellIndexQ,
                                                                                  0.0);
                            
                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;

                                    //const index_type shellIndexRS = this->index(shellIndexR, shellIndexS);
                                    const DfEriEngine::CGTO_Pair RS = engine.getCGTO_pair(orbitalInfo,
                                                                                          shellIndexR, shellIndexS,
                                                                                          0.0);
                                    
                                    engine.calc(queryPQ, queryRS, PQ, RS);

                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const index_type indexP = shellIndexP + i;
                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const index_type indexQ = shellIndexQ + j;
                                            const index_type indexPQ = this->index(indexP, indexQ);
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;
                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    const index_type indexRS = this->index(indexR, indexS);
                                                    
                                                    const double value = engine.WORK[index];

                                                    G.set(indexPQ, indexRS, value);
                                                    ++index;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    G.save("G_exact.mat");

    TlMatrix L = G.choleskyFactorization2(this->epsilon_);
    L.save("L_exact.mat");

    // check
    {
        TlMatrix Lt = L;
        Lt.transpose();
        
        TlMatrix LLt = L * Lt;
        LLt.save("LLt_exact.mat");
    }
}

DfCD::ShellArrayTable DfCD::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();

    index_type shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = orbitalInfo.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        shellArrayTable[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }

    return shellArrayTable;
}

DfCD::ShellPairArrayTable DfCD::getShellPairArrayTable(const ShellArrayTable& shellArrayTable)
{
    ShellPairArrayTable shellPairArrayTable(MAX_SHELL_TYPE * MAX_SHELL_TYPE);

    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
            const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
            ShellArray::const_iterator rItEnd = shellArrayR.end();

            const int shellPairType_PR = shellTypeP * MAX_SHELL_TYPE + shellTypeR;
            for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                const index_type indexP = *pIt;

                for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                    const index_type indexR = *rIt;

                    if (indexP >= indexR) {
                        ShellPair shellPair(indexP, indexR);
                        shellPairArrayTable[shellPairType_PR].push_back(shellPair);
                    }
                }
            }
        }
    }

    return shellPairArrayTable;
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
    return pDfTaskCtrl;
}

void DfCD::finalize(TlMatrix* pMtx)
{
    // do nothing
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
                                          const std::vector<int>& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlVector buf(this->numOfPQs_);
    for (index_type i = 0; i < numOfItilde; ++i) {
        const index_type index = I2PQ[i];
        buf[index] = L_col[i];
    }

    return TlSymmetricMatrix(buf, this->m_nNumOfAOs);
}

void DfCD::getJ(TlSymmetricMatrix* pJ)
{
    const index_type numOfAOs = this->m_nNumOfAOs;

    const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

    // cholesky vector
    TlMatrix L;
    L.load("L.mat");
    const index_type numOfCBs = L.getNumOfCols();

    const std::vector<int> I2PQ = this->getI2PQ();
    for (index_type I = 0; I < numOfCBs; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == numOfAOs);
        assert(LI.getNumOfCols() == numOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();
        
        *pJ += qi*LI;
    }
}


void DfCD::getK(const RUN_TYPE runType,
                TlSymmetricMatrix *pK)
{
    TlMatrix L;
    L.load("L.mat");
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    const std::vector<int> I2PQ = this->getI2PQ();
    for (index_type I = 0; I < numOfCBs; ++I) {
        TlSymmetricMatrix l = this->getCholeskyVector(L.getColVector(I), I2PQ);;
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
}


// void DfCD::getJ(TlSymmetricMatrix* pJ)
// {
//     const index_type numOfAOs = this->m_nNumOfAOs;

//     const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

//     // cholesky vector
//     TlMatrix L;
//     L.load("L.mat");
//     const index_type numOfCBs = L.getNumOfCols();
    
//     for (index_type I = 0; I < numOfCBs; ++I) {
//         TlSymmetricMatrix LI(numOfAOs);
//         for (index_type p = 0; p < numOfAOs; ++p) {
//             for (index_type q = 0; q <= p; ++q) {
//                 const index_type rs = p + (2 * numOfAOs - (q +1)) * q / 2;
//                 LI.set(p, q, L.get(rs, I));
//             }
//         }

//         TlMatrix QI = LI;
//         QI.dot(P);
//         double qi = QI.sum();

//         *pJ += qi*LI;
//     }
// }


// void DfCD::getK(const RUN_TYPE runType,
//                 TlSymmetricMatrix *pK)
// {
//     const index_type numOfAOs = this->m_nNumOfAOs;
    
//     TlMatrix L;
//     L.load("L.mat");
//     const index_type numOfCBs = L.getNumOfCols();
    
//     TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
//     const TlMatrix C = P.choleskyFactorization2();
    
//     for (index_type I = 0; I < numOfCBs; ++I) {
//         TlSymmetricMatrix l(numOfAOs);
//         for (index_type p = 0; p < numOfAOs; ++p) {
//             for (index_type q = 0; q <= p; ++q) {
//                 const index_type index = p + (2 * numOfAOs - (q +1)) * q / 2;
//                 l.set(p, q, L.get(index, I));
//             }
//         }
    
//         TlMatrix X = l * C;
//         TlMatrix Xt = X;
//         Xt.transpose();
        
//         TlSymmetricMatrix XX = X * Xt;
//         *pK += XX;
//     }
    
//     *pK *= -1.0;
// }
