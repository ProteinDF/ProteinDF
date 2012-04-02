#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <set>
#include "DfCD.h"
#include "DfEriEngine.h"
#include "TlOrbitalInfo.h"
#include "TlSymmetricMatrix.h"

const int DfCD::MAX_SHELL_TYPE = 2 + 1;


DfCD::DfCD(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEriEngines_(NULL)
{
    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
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
    
    const index_type numOfAOs = this->m_nNumOfAOs;
    const std::size_t dim = numOfAOs * (numOfAOs +1) / 2;
    TlSymmetricMatrix G(dim);

    this->check.resize(dim);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

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
    this->check.save("check.mat");
    
    // 正定値(固有値が正)であるかをチェック
    {
        TlMatrix eigvec;
        TlVector eigval;
        G.diagonal(&eigval, &eigvec);
        eigval.save("G_eigval.vec");
        eigvec.save("G_eigvec.mat");
    }

    TlMatrix L = G.choleskyFactorization2();
    L.save("L0.mat");

    // std::cerr << "dim=" << dim << std::endl;
    // std::cerr << "L size=(" << L.getNumOfRows()
    //           << ", " << L.getNumOfCols()
    //           << ")" << std::endl;
    // for (int i = 0; i < dim; ++i) {
    //     // for (int j = 0; j <= i; ++j) {
    //     //     const double v = L.get(i, j);
    //     //     L.set(i, j, v);
    //     // }
    //     for (int j = i +1; j < dim; ++j) {
    //         L.set(i, j, 0.0);
    //     }
    // }
    // L.save("L.mat");
    TlMatrix Lt = L;
    Lt.transpose();
    
    TlMatrix LL = L * Lt;
    LL.save("LL.mat");
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
                    
                    if ((indexP >= indexQ) && (maxIndexS >= indexS)) {
                        if ((this->index(shellIndexP, shellIndexQ) != this->index(shellIndexR, shellIndexS)) |
                            (indexPQ >= indexRS)) {
                            this->check.add(indexPQ, indexRS, 1.0);
                            pG->set(indexPQ, indexRS, value);
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
    std::set<index_type> calcd; // 計算済みかどうかを格納
    
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

                            const index_type shellIndexPQ = this->index(shellIndexP, shellIndexQ);
                            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                                  shellIndexP, shellIndexQ,
                                                                                  0.0);
                            
                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;

                                    const index_type shellIndexRS = this->index(shellIndexR, shellIndexS);
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

                                                    // check
                                                    if (((indexPQ == 115) && (indexRS == 111)) ||
                                                        ((indexPQ == 111) && (indexRS == 115))) {
                                                        std::cerr << TlUtils::format("CHECK <%d, %d> (%d, %d|%d, %d) [%d, %d|%d, %d]= %18.10f",
                                                                                     indexPQ, indexRS,
                                                                                     indexP, indexQ, indexR, indexS,
                                                                                     shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                                                                                     value)
                                                                  << std::endl;
                                                    }
                                                    // if ((shellIndexP == 9) && (shellIndexQ == 9) &&
                                                    //     (shellIndexR == 9) && (shellIndexS == 9)) {
                                                    //     std::cerr << TlUtils::format("CHECK <%2d, %2d> (%2d, %2d|%2d, %2d) [%2d, %2d|%2d, %2d]= %18.10f",
                                                    //                                  indexPQ, indexRS,
                                                    //                                  indexP, indexQ, indexR, indexS,
                                                    //                                  shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                                                    //                                  value)
                                                    //               << std::endl;
                                                    // }

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
    this->checkCholeskyFactorization(G);
    // 正定値(固有値が正)であるかをチェック
    {
        TlMatrix eigvec;
        TlVector eigval;
        G.diagonal(&eigval, &eigvec);
        eigval.save("G_exact_eigval.vec");
        eigvec.save("G_exact_eigvec.mat");
    }
}

void DfCD::checkCholeskyFactorization(const TlSymmetricMatrix& G)
{
    const index_type dim = G.getNumOfRows();
    
    FittingTbl ft(dim);
    {
        const index_type numOfAOs = this->m_nNumOfAOs;
        index_type index = 0;
        for (index_type i = 0; i < numOfAOs; ++i) {
            for (index_type j = 0; j <= i; ++j) {
                ft[index].p = i;
                ft[index].q = j;
                this->log_.debug(TlUtils::format("(%d, %d) -> %d", i, j, index));
                ++index;
            }
        }
        assert(index == dim);
    }

    
    std::vector<double> diagonals(dim);
    for (index_type i = 0; i < dim; ++i) {
        diagonals[i] = G.get(i, i);
    }

    for (index_type i = 0; i < dim; ++i) {
        const double max_i = diagonals[i];
        for (index_type j = 0; j < i; ++j) {
            const double max_j = diagonals[j];
            const double value = G.get(i, j);
            if (std::fabs(value) < 1.0E-16) {
                continue;
            }
            if (value > std::max(max_i, max_j)) {
                this->log_.debug(TlUtils::format("CD_check G((%d, %d), (%d, %d))=%16.10f > max(%16.10f, %16.10f)",
                                                 ft[i].p, ft[i].q,
                                                 ft[j].p, ft[j].q,
                                                 value,
                                                 max_i, max_j));
            }
        }
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
