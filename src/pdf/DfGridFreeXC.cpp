#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfGridFreeXC.h"
#include "DfFunctional_SVWN.h"

const int DfGridFreeXC::MAX_SHELL_TYPE = 2 + 1;

DfGridFreeXC::DfGridFreeXC(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfGridFreeXC::~DfGridFreeXC()
{
}


void DfGridFreeXC::buildFxc()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const TlSymmetricMatrix P = 0.5 * DfObject::getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    TlSymmetricMatrix S = DfObject::getSpqMatrix<TlSymmetricMatrix>();

    TlSymmetricMatrix M;
    this->getM_exact(P, &M);
    M.save("M_exact.mat");
    // this->getM(P, &M);
    // M.save("M.mat");

    // V version
    {
        TlMatrix V = DfObject::getXMatrix<TlMatrix>();
        TlMatrix tV = V;
        tV.transpose();
        {
            // check
            TlMatrix VSV = tV * S * V;
            VSV.save("VSV.mat"); // = I
        }
        
        TlSymmetricMatrix M_tilda = tV * M * V;
        M_tilda.save("M_tilda.mat");
        
        TlMatrix U;
        TlVector lamda;
        M_tilda.diagonal(&lamda, &U);
        U.save("U.mat");
        lamda.save("lamda.vct");
        
        TlSymmetricMatrix F_lamda = this->get_F_lamda(lamda);
        F_lamda.save("F_lamda.mat");
        
        TlMatrix tU = U;
        tU.transpose();
        
        TlSymmetricMatrix Fxc = S * V * U * F_lamda * tU * tV * S;
        DfObject::saveFxcMatrix(RUN_RKS, this->m_nIteration, Fxc);

        // check 
        // {
        //     TlSymmetricMatrix L(numOfAOs);
        //     for (int i = 0; i < numOfAOs; ++i) {
        //         L(i, i) = lamda[i];
        //     }
        
        //     TlSymmetricMatrix SVULUVS = S * V * U * L * tU * tV * S;
        //     SVULUVS.save("SVULUVS_M.mat");
        // }
    }

    // X version
    // {
    //     TlSymmetricMatrix X;
    //     {
    //         TlMatrix u;
    //         TlVector eigval;
    //         S.diagonal(&eigval, &u);
            
    //         TlSymmetricMatrix s(numOfAOs);
    //         for (int i = 0; i < numOfAOs; ++i) {
    //             s(i, i) = 1.0 / std::sqrt(eigval[i]);
    //         }
    //         TlMatrix tu = u;
    //         tu.transpose();
    //         X = u * s * tu;
    //     }
    //     X.save("X.mat");
    //     // check
    //     {
    //         TlMatrix XSX = X * S * X;
    //         XSX.save("XSX.mat");
    //     }
        
    //     TlSymmetricMatrix M_tilda = X * M * X;
    //     M_tilda.save("XMX.mat");
    //     TlMatrix U;
    //     TlVector lamda;
    //     M_tilda.diagonal(&lamda, &U);
    //     TlSymmetricMatrix F_lamda = this->get_F_lamda(lamda);

    //     TlSymmetricMatrix L(numOfAOs);
    //     for (int i = 0; i < numOfAOs; ++i) {
    //         L(i, i) = lamda[i];
    //     }

    //     TlMatrix tU = U;
    //     tU.transpose();

    //     TlSymmetricMatrix Fxc = X * S * U * F_lamda * tU * S * X;
    //     //TlSymmetricMatrix Fxc = S * X * U * F_lamda * tU * X * S;
    //     Fxc.save("Fxc_X.mat");

    //     TlSymmetricMatrix XSULUSX = X * S * U * L * tU * S * X;
    //     XSULUSX.save("XSULUSX_M.mat");

    //     TlSymmetricMatrix PF = P * Fxc;
    //     PF.save("PF.mat");
    // }
}


void DfGridFreeXC::getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    assert(pM != NULL);
    pM->resize(this->m_nNumOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    // pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    // pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    // pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_,
                                          &taskList,
                                          true);
    while (hasTask == true) {
        this->getM_part(orbitalInfo,
                        taskList,
                        P, pM);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_,
                                         &taskList);
    }
                    
    this->finalize(pM);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfGridFreeXC::createEngines()
{
    assert(this->pOvpEngines_ == NULL);
    
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->log_.info(TlUtils::format("create OpenMP ERI engine: %d", numOfThreads));
        this->pOvpEngines_ = new DfOverlapEngine[numOfThreads];
    }
#else
    this->pOvpEngines_ = new DfOverlapEngine[1];
#endif // _OPENMP
}


void DfGridFreeXC::destroyEngines()
{
    this->log_.info("delete OpenMP ERI engine");
    if (this->pOvpEngines_ != NULL) {
        delete[] this->pOvpEngines_;
    }
    this->pOvpEngines_ = NULL;
}


DfTaskCtrl* DfGridFreeXC::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfGridFreeXC::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}


TlSparseSymmetricMatrix DfGridFreeXC::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo)
{
    this->log_.info("make Schwartz cutoff table: start");
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    TlSparseSymmetricMatrix schwarz(maxShellIndex);

    DfOverlapEngine engine;
    // engine.setPrimitiveLevelThreshold(0.0);
    
    for (index_type shellIndexP = 0; shellIndexP < maxShellIndex; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;
        const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
        const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

        for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);
            const DfOverlapEngine::Query query(0, 0, 0, 0,
                                               shellTypeP, shellTypeQ,
                                               shellTypeP, shellTypeQ);

            engine.calc(query, posP, posQ, posP, posQ,
                        pgtosP, pgtosQ, pgtosP, pgtosQ);

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

    this->log_.info("make Schwartz cutoff table: end");
    return schwarz;
}


void DfGridFreeXC::getM_part(const TlOrbitalInfoObject& orbitalInfo,
                             const std::vector<DfTaskCtrl::Task4>& taskList,
                             const TlMatrixObject& P, TlMatrixObject* pM)
{
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        // this->pOvpEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for schedule(runtime)
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
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
            const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
            const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
            const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);
                        
            const DfOverlapEngine::Query query(0, 0, 0, 0,
                                               shellTypeP, shellTypeQ,
                                               shellTypeR, shellTypeS);

            this->pOvpEngines_[threadID].calc(query,
                                              posP, posQ, posR, posS,
                                              pgtosP, pgtosQ, pgtosR, pgtosS);
                        
            this->storeM(shellIndexP, maxStepsP,
                         shellIndexQ, maxStepsQ,
                         shellIndexR, maxStepsR,
                         shellIndexS, maxStepsS,
                         this->pOvpEngines_[threadID], P, pM);
        }
    }
        
}


void DfGridFreeXC::storeM(const index_type shellIndexP, const int maxStepsP,
                          const index_type shellIndexQ, const int maxStepsQ,
                          const index_type shellIndexR, const int maxStepsR,
                          const index_type shellIndexS, const int maxStepsS,
                          const DfOverlapEngine& engine,
                          const TlMatrixObject& P,
                          TlMatrixObject* pM)
{
    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;

        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            const double P_pq = P.get(indexP, indexQ);
            
            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;

                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    const double P_rs = P.get(indexR, indexS);
                    
                    const double value = engine.WORK[index];
                    const index_type maxIndexS = (indexP == indexR) ? indexQ : indexR;
                    
                    if ((indexP >= indexQ) && (maxIndexS >= indexS)) {
                        // Eq.1 : (indexP, indexQ) <= (indexR, indexS)
                        const double coefEq1 = (indexR != indexS) ? 2.0 : 1.0;
                        pM->add(indexP, indexQ, coefEq1 * P_rs * value);
                        
                        // Eq.2 : (indexR, indexS) <= (indexP, indexQ)
                        if ((shellIndexP != shellIndexR) || (shellIndexQ != shellIndexS) || (indexP == indexR)) {
                            if (((indexP + indexQ) != (indexR + indexS)) ||
                                ((indexP * indexQ) != (indexR * indexS))) {
                                // Eq.1の条件と重複しないようにするための措置
                                            
                                const double coefEq2 = (indexP != indexQ) ? 2.0 : 1.0;
                                pM->add(indexR, indexS, coefEq2 * P_pq * value);
                            }
                        }
                    }
                    ++index;
                }
            }
        }
    }
}


TlSymmetricMatrix DfGridFreeXC::get_F_lamda(const TlVector lamda)
{
    const int dim = lamda.getSize();
    TlSymmetricMatrix F_lamda(dim);

    DfFunctional_LDA* pFunc = new DfFunctional_SVWN();

    double fv_a = 0.0;
    double fv_b = 0.0;
    for (int i = 0; i < dim; ++i) {
        const double v = lamda.get(i);
        if (v > 1.0E-16) {
            pFunc->getDerivativeFunctional(v, v, &fv_a, &fv_b);
            F_lamda.set(i, i, fv_a);
        }
    }

    delete pFunc;
    pFunc = NULL;

    return F_lamda;
}


void DfGridFreeXC::getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    assert(pM != NULL);
    pM->resize(this->m_nNumOfAOs);
    TlMatrix tmpM(this->m_nNumOfAOs, this->m_nNumOfAOs);

    DfOverlapEngine engine;
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    // for debug
    // orbitalInfo.printCGTOs(std::cerr);

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
            
            for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = 0; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfOverlapEngine::Query query(0, 0, 0, 0,
                                                       shellTypeP, shellTypeQ,
                                                       shellTypeR, shellTypeS);

                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
                        const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;
                            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
                            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
                                const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
                                
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;
                                    const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
                                    const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);

                                    engine.calc(query,
                                                posP, posQ, posR, posS,
                                                pgtosP, pgtosQ, pgtosR, pgtosS);
                                    
                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const int indexP = shellIndexP + i;
                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const int indexQ = shellIndexQ + j;
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;
                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    
                                                    if (indexP >= indexQ) {
                                                        const double P_rs = P.get(indexR, indexS);
                                                        const double value = engine.WORK[index];
                                                        pM->add(indexP, indexQ, P_rs * value);
                                                    }

                                                    // for debug
                                                    {
                                                        const double P_rs = P.get(indexR, indexS);
                                                        const double value = engine.WORK[index];
                                                        // std::cerr << TlUtils::format("M(%d %d) P(%d %d)=% e (%d %d %d %d)=% e",
                                                        //                              indexP, indexQ, indexR, indexS, P_rs,
                                                        //                              indexP, indexQ, indexR, indexS, value)
                                                        //           << std::endl;
                                                        tmpM.add(indexP, indexQ, P_rs * value);
                                                    }
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

    tmpM.save("tmpM.mat");
}


// void DfGridFreeXC::getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
// {
//     assert(pM != NULL);
//     pM->resize(this->m_nNumOfAOs);

//     TlSymmetricMatrix S(this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2);
    
//     DfOverlapEngine engine;
    
//     const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
//                                     (*(this->pPdfParam_))["basis_sets"]);
//     const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
//     const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);

//     for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
//         const int maxStepsP = 2 * shellTypeP + 1;
//         const ShellArray shellArrayP = shellArrayTable[shellTypeP];
//         ShellArray::const_iterator pItEnd = shellArrayP.end();

//         for (int shellTypeQ = MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
//             const int maxStepsQ = 2 * shellTypeQ + 1;
//             const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
//             ShellArray::const_iterator qItEnd = shellArrayQ.end();
            
//             for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
//                 const int maxStepsR = 2 * shellTypeR + 1;
//                 const ShellArray shellArrayR = shellArrayTable[shellTypeR];
//                 ShellArray::const_iterator rItEnd = shellArrayR.end();

//                 for (int shellTypeS = 0; shellTypeS >= 0; --shellTypeS) {
//                     const int maxStepsS = 2 * shellTypeS + 1;
//                     const ShellArray shellArrayS = shellArrayTable[shellTypeS];
//                     ShellArray::const_iterator sItEnd = shellArrayS.end();

//                     const DfOverlapEngine::Query query(0, 0, 0, 0,
//                                                        shellTypeP, shellTypeQ,
//                                                        shellTypeR, shellTypeS);

//                     for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
//                         const index_type shellIndexP = *pIt;
//                         const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
//                         const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

//                         for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
//                             const index_type shellIndexQ = *qIt;
//                             const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
//                             const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

//                             for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
//                                 const index_type shellIndexR = *rIt;
//                                 const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
//                                 const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
                                
//                                 for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
//                                     const index_type shellIndexS = *sIt;
//                                     const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
//                                     const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);

//                                     engine.calc(query,
//                                                 posP, posQ, posR, posS,
//                                                 pgtosP, pgtosQ, pgtosR, pgtosS);
                                    
//                                     int index = 0;
//                                     for (int i = 0; i < maxStepsP; ++i) {
//                                         const int indexP = shellIndexP + i;
//                                         for (int j = 0; j < maxStepsQ; ++j) {
//                                             const int indexQ = shellIndexQ + j;
                                            
//                                             for (int k = 0; k < maxStepsR; ++k) {
//                                                 const int indexR = shellIndexR + k;
//                                                 for (int l = 0; l < maxStepsS; ++l) {
//                                                     const int indexS = shellIndexS + l;
                                                    
//                                                     if (indexP >= indexQ) {
//                                                         const double P_rs = P.get(indexR, indexS);
//                                                         const double value = engine.WORK[index];
//                                                         pM->add(indexP, indexQ, P_rs * value);
//                                                     }
                                                    
//                                                     ++index;
//                                                 }
//                                              }
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }


DfGridFreeXC::ShellArrayTable DfGridFreeXC::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
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


DfGridFreeXC::ShellPairArrayTable DfGridFreeXC::getShellPairArrayTable(const ShellArrayTable& shellArrayTable)
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
