#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfOverlapX.h"

const int DfOverlapX::MAX_SHELL_TYPE = 2 + 1;

DfOverlapX::DfOverlapX(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam)
{
#ifdef _OPENMP
    const int numOfThreads = omp_get_max_threads();
    this->pEngines_ = new DfOverlapEngine[numOfThreads];
#else
    this->pEngines_ = new DfOverlapEngine[1];
#endif // _OPENMP
}


DfOverlapX::~DfOverlapX()
{
    delete[] this->pEngines_;
    this->pEngines_ = NULL;
}


DfTaskCtrl* DfOverlapX::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfOverlapX::finalize(TlMatrix* pMtx)
{
    // do nothing
}


void DfOverlapX::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}


void DfOverlapX::finalize(TlVector* pVct)
{
    // do nothing
}


void DfOverlapX::getSpq(TlSymmetricMatrix* pSpq)
{
    assert(pSpq != NULL);
    const index_type numOfAOs = this->m_nNumOfAOs;
    pSpq->resize(numOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["model"]["coordinates"],
                                    (*(this->pPdfParam_))["model"]["basis_set"]);
    this->calcOverlap(orbitalInfo, pSpq);
    this->finalize(pSpq);
}


void DfOverlapX::getSab(TlSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pSab->resize(numOfAuxDens);

    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["model"]["coordinates"],
                                                    (*(this->pPdfParam_))["model"]["basis_set_auxD"]);
    this->calcOverlap(orbitalInfo_Density, pSab);
    this->finalize(pSab);
}



void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrixObject* pMatrix)
{
    DfTaskCtrl*  pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue(orbitalInfo,
                                       false,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo,
                               taskList,
                               pMatrix);
        
        hasTask = pTaskCtrl->getQueue(orbitalInfo,
                                      false,
                                      this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    
    delete pTaskCtrl;
    pTaskCtrl = NULL;
}


void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task2>& taskList,
                                  TlMatrixObject* pMatrix)
{
    // 第三中心点の固定値を設定
    static const TlPosition posR(0.0, 0.0, 0.0);
    static DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    static DfOverlapEngine::PGTOs pgtosR(1);
    pgtosR[0] = pgtoR;

    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosP = this->getPGTOs(orbitalInfo, shellIndexP);
            const DfOverlapEngine::PGTOs pgtosQ = this->getPGTOs(orbitalInfo, shellIndexQ);
            const DfOverlapEngine::Query query(0, 0, 0, shellTypeP, shellTypeQ, 0);
            
            this->pEngines_[threadID].calc(query, posP, posQ, posR, pgtosP, pgtosQ, pgtosR);
            
            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        pMatrix->add(globalShellIndexP, globalShellIndexQ, this->pEngines_[threadID].WORK[index]);
                    }
                    ++index;
                }
            }
        }
    }
}

void DfOverlapX::getForce(const TlSymmetricMatrix& W,
                          TlMatrix* pForce)
{
    assert(pForce != NULL);
    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["model"]["coordinates"],
                                    (*(this->pPdfParam_))["model"]["basis_set"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);

    
    for (int shellTypeP = DfOverlapX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        const index_type shellArraySizeP = shellArrayP.size();
        
        for (int shellTypeQ = DfOverlapX::MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            // const index_type shellArraySizeQ = shellArrayQ.size();
            
            for (index_type p = 0; p < shellArraySizeP; ++p) {
                const index_type shellIndexP = shellArrayP[p];

                this->getForce_partProc(orbitalInfo,
                                        shellTypeP, shellTypeQ,
                                        shellIndexP,
                                        shellArrayQ,
                                        W, pForce);
            }
        }
    }
}


void DfOverlapX::getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                                   const int shellTypeP, const int shellTypeQ,
                                   const index_type shellIndexP,
                                   const ShellArray& shellArrayQ,
                                   const TlSymmetricMatrix& W,
                                   TlMatrix* pForce)
{
    // 第三中心点の固定値を設定
    static const TlPosition posR(0.0, 0.0, 0.0);
    static const DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    static DfOverlapEngine::PGTOs pgtosR(1);
    pgtosR[0] = pgtoR;

    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const std::size_t shellArraySizeQ = shellArrayQ.size();
    const DfOverlapEngine::Query query(1, 0, 0, shellTypeP, shellTypeQ, 0);
    
    const DfOverlapEngine::PGTOs pgtosP = this->getPGTOs(orbitalInfo,
                                                         shellIndexP);
    const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
    const int atomIndexA = orbitalInfo.getAtomIndex(shellIndexP);

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t q = 0; q < shellArraySizeQ; ++q) {
            const index_type shellIndexQ = shellArrayQ[q];
            // 重複ペアの排除
//                     if (shellIndexP < shellIndexQ) {
//                         continue;
//                     }
                    
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosQ = this->getPGTOs(orbitalInfo,
                                                                 shellIndexQ);
            const int atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
                
            this->pEngines_[threadID].calc(query, posP, posQ, posR, pgtosP, pgtosQ, pgtosR);

            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, X, coef * dSdA);
                        pForce->add(atomIndexB, X, coef * dSdB);
                    }
                    ++index;
                }
            }
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, Y, coef * dSdA);
                        pForce->add(atomIndexB, Y, coef * dSdB);
                    }
                    ++index;
                }
            }
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, Z, coef * dSdA);
                        pForce->add(atomIndexB, Z, coef * dSdB);
                    }
                    ++index;
                }
            }
        }
    }
}


DfOverlapX::ShellArrayTable DfOverlapX::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
    const int maxShellIndex = orbitalInfo.getNumOfOrbitals();

    int shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = orbitalInfo.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        shellArrayTable[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }

    return shellArrayTable;
}


DfOverlapEngine::PGTOs DfOverlapX::getPGTOs(const TlOrbitalInfoObject& orbitalInfo,
                                            const int shellIndex)
{
    DfOverlapEngine::PGTOs pgtos;

    const int numOfContractions = orbitalInfo.getCgtoContraction(shellIndex);
    pgtos.resize(numOfContractions);
    for (int i = 0; i < numOfContractions; ++i) {
        const DfOverlapEngine::PGTO pgto(orbitalInfo.getCoefficient(shellIndex, i),
                                         orbitalInfo.getExponent(shellIndex, i));
        pgtos[i] = pgto;
    }
    
    return pgtos;
}


