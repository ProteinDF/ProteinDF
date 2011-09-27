#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfHpqX.h"

const int DfHpqX::MAX_SHELL_TYPE = 2 +1;

DfHpqX::DfHpqX(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["model"]["coordinates"],
                   (*pPdfParam)["model"]["basis_set"]) {
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->pEngines_ = new DfHpqEngine[numOfThreads];
    }
#else
    this->pEngines_ = new DfHpqEngine[1];
#endif // _OPENMP
        
    this->makeShellArrayTable();
}


DfHpqX::~DfHpqX()
{
    delete[] this->pEngines_;
    this->pEngines_ = NULL;
}


DfTaskCtrl* DfHpqX::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfHpqX::finalize(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2)
{
    // do nothing
}


void DfHpqX::getHpq(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2)
{
    const int numOfAOs = this->m_nNumOfAOs;
    
    // make coordinates
    const Fl_Geometry flGeom((*this->pPdfParam_)["model"]["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    const int numOfDummyAtoms = flGeom.getDummyatom();
    const int numOfRealAtoms = numOfAtoms - numOfDummyAtoms;
    std::vector<TlAtom> Cs(numOfRealAtoms);
    std::vector<TlAtom> Xs(numOfDummyAtoms);
    std::size_t realAtomIndex = 0;
    std::size_t dummyAtomIndex = 0;
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string atomName = flGeom.getAtom(i);
        const TlPosition p = flGeom.getCoordinate(i);
        const double charge = flGeom.getCharge(i);
        const TlAtom atom(atomName, p, charge);
        if (atomName == "X") {
            Xs[dummyAtomIndex] = atom;
            ++dummyAtomIndex;
        } else {
            Cs[realAtomIndex] = atom;
            ++realAtomIndex;
        }
    }
    assert(realAtomIndex == Cs.size());
    assert(dummyAtomIndex == Xs.size());
    
    DfHpqEngine engine;

    pHpq->resize(numOfAOs);
    pHpq2->resize(numOfAOs);

    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue(this->orbitalInfo_,
                                       false,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getHpq_part(this->orbitalInfo_,
                          taskList,
                          Cs, Xs,
                          pHpq, pHpq2);
        
        hasTask = pTaskCtrl->getQueue(this->orbitalInfo_,
                                      false,
                                      this->grainSize_, &taskList);
    } 

    if (this->chargeExtrapolateNumber_ > 0) {
        *pHpq2 /= static_cast<double>(this->chargeExtrapolateNumber_);
    }

    this->finalize(pHpq, pHpq2);
    pTaskCtrl->cutoffReport();
    
    delete pTaskCtrl;
}


void DfHpqX::getHpq_part(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         const std::vector<TlAtom>& Cs,
                         const std::vector<TlAtom>& Xs,
                         TlMatrixObject* pHpq,
                         TlMatrixObject* pHpq2)
{
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
            const DfHpqEngine::PGTOs pgtosP = this->getPGTOs(shellIndexP);
            const DfHpqEngine::PGTOs pgtosQ = this->getPGTOs(shellIndexQ);
            const DfHpqEngine::Query query(0, 0, shellTypeP, shellTypeQ); 

            this->pEngines_[threadID].calc(query, posP, posQ, pgtosP, pgtosQ, Cs, Xs);

            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;
                    
                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        const double coef = 1.0;
                        const double value = coef * (  this->pEngines_[threadID].WORK_KIN[index]
                                                       + this->pEngines_[threadID].WORK_NUC[index]);
                        pHpq->add(globalShellIndexP, globalShellIndexQ, value);
                        
                        const double valueFromDummyCharge = coef * this->pEngines_[threadID].WORK_NUCX[index];
                        pHpq2->add(globalShellIndexP, globalShellIndexQ, valueFromDummyCharge);
                    }
                    ++index;
                }
            }
        }
    }
}


void DfHpqX::getForce(const TlSymmetricMatrix& P,
                      TlMatrix* pForce)
{
    assert(pForce != NULL);
    assert(pForce->getNumOfRows() == this->m_nNumOfAtoms);
    assert(pForce->getNumOfCols() == 3);
    
    //DfHpqEngine engine;
    
    for (int shellTypeP = DfHpqX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray shellArrayP = this->shellArrayTable_[shellTypeP];
        const std::size_t shellArraySizeP = shellArrayP.size();
        
        for (int shellTypeQ = DfHpqX::MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const ShellArray shellArrayQ = this->shellArrayTable_[shellTypeQ];
            // const std::size_t shellArraySizeQ = shellArrayQ.size();
            
            for (std::size_t p = 0; p < shellArraySizeP; ++p) {
                const index_type shellIndexP = shellArrayP[p];

                this->getForce_partProc(this->orbitalInfo_,
                                        shellTypeP, shellTypeQ,
                                        shellIndexP,
                                        shellArrayQ,
                                        P, pForce);
            }
        }
    }

}


void DfHpqX::getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                               const int shellTypeP, const int shellTypeQ,
                               const index_type shellIndexP,
                               const ShellArray& shellArrayQ,
                               const TlMatrixObject& P,
                               TlMatrix* pForce)
{
    static const int BUFFER_SIZE_NUC = 3 * 5 * 5 * 5; // (xyz) * 5d * 5d * 5d
    const Fl_Geometry flGeom((*this->pPdfParam_)["model"]["coordinates"]);

    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    
    const DfHpqEngine::Query queryPQ10(1, 0, shellTypeP, shellTypeQ);
    const DfHpqEngine::Query queryPQ01(0, 1, shellTypeP, shellTypeQ);

    const TlPosition posP = orbitalInfo_.getPosition(shellIndexP);
    const index_type atomIndexA = orbitalInfo_.getAtomIndex(shellIndexP);

    const int numOfAtoms = this->m_nNumOfAtoms;
    const DfHpqEngine::PGTOs pgtosP = this->getPGTOs(shellIndexP);
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        double* p_dNuc_dA = new double[BUFFER_SIZE_NUC];
        double* p_dNuc_dB = new double[BUFFER_SIZE_NUC];
        
        const std::size_t shellArraySizeQ = shellArrayQ.size();
#pragma omp for schedule(runtime)
        for (std::size_t q = 0; q < shellArraySizeQ; ++q) {
            const index_type shellIndexQ = shellArrayQ[q];
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfHpqEngine::PGTOs pgtosQ = this->getPGTOs(shellIndexQ);
            const index_type atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
            
            // 重複ペアの排除
//                     if (shellIndexP < shellIndexQ) {
//                         continue;
//                     }
            
            // 運動エネルギー
            this->pEngines_[threadID].calcKineticPart(queryPQ10, posP, posQ, pgtosP, pgtosQ);
            {
                int index = 0;
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double dKin_dA = this->pEngines_[threadID].WORK_KIN[index];
                        const double dKin_dB = - dKin_dA;
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, X, coef * dKin_dA);
                            pForce->add(atomIndexB, X, coef * dKin_dB);
                        }
                        ++index;
                    }
                }
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double dKin_dA = this->pEngines_[threadID].WORK_KIN[index];
                        const double dKin_dB = - dKin_dA;
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, Y, coef * dKin_dA);
                            pForce->add(atomIndexB, Y, coef * dKin_dB);
                        }
                        ++index;
                    }
                }
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double dKin_dA = this->pEngines_[threadID].WORK_KIN[index];
                        const double dKin_dB = - dKin_dA;
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, Z, coef * dKin_dA);
                            pForce->add(atomIndexB, Z, coef * dKin_dB);
                        }
                        ++index;
                    }
                }
            }
            
            // 核-電子反発
            for (index_type atomIndexC = 0; atomIndexC < numOfAtoms; ++atomIndexC) {
                const std::string atomName = flGeom.getAtom(atomIndexC);
                const TlPosition p = flGeom.getCoordinate(atomIndexC);
                const double charge = flGeom.getCharge(atomIndexC);
                const TlAtom C(atomName, p, charge);
                std::vector<TlAtom> Cs(1);
                Cs[0] = C;
                
                this->pEngines_[threadID].calcNuclearAttractionPart(queryPQ10, posP, posQ, pgtosP, pgtosQ, Cs);
                std::copy(this->pEngines_[threadID].WORK_NUC,
                          this->pEngines_[threadID].WORK_NUC + BUFFER_SIZE_NUC, p_dNuc_dA);
                this->pEngines_[threadID].calcNuclearAttractionPart(queryPQ01, posP, posQ, pgtosP, pgtosQ, Cs);
                std::copy(this->pEngines_[threadID].WORK_NUC,
                          this->pEngines_[threadID].WORK_NUC + BUFFER_SIZE_NUC, p_dNuc_dB);
                
                int index = 0;
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double gradA = p_dNuc_dA[index];
                        const double gradB = p_dNuc_dB[index];
                        const double gradC = - (gradA + gradB);
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, X, coef * gradA);
                            pForce->add(atomIndexB, X, coef * gradB);
                            pForce->add(atomIndexC, X, coef * gradC);
                        }
                        ++index;
                    }
                }
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double gradA = p_dNuc_dA[index];
                        const double gradB = p_dNuc_dB[index];
                        const double gradC = - (gradA + gradB);
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, Y, coef * gradA);
                            pForce->add(atomIndexB, Y, coef * gradB);
                            pForce->add(atomIndexC, Y, coef * gradC);
                        }
                        ++index;
                    }
                }
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type indexP = shellIndexP + stepP;
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type indexQ = shellIndexQ + stepQ;
                        
                        double coef = P.get(indexP, indexQ);
                        const double gradA = p_dNuc_dA[index];
                        const double gradB = p_dNuc_dB[index];
                        const double gradC = - (gradA + gradB);
#pragma omp critical(DfHpqX__getForce_partProc)
                        {
                            pForce->add(atomIndexA, Z, coef * gradA);
                            pForce->add(atomIndexB, Z, coef * gradB);
                            pForce->add(atomIndexC, Z, coef * gradC);
                        }
                        ++index;
                    }
                }
            }
        }

        delete[] p_dNuc_dA;
        p_dNuc_dA = NULL;
        delete[] p_dNuc_dB;
        p_dNuc_dB = NULL;
    }
}


void DfHpqX::makeShellArrayTable()
{
    this->shellArrayTable_.resize(MAX_SHELL_TYPE);
    const index_type maxShellIndex = this->orbitalInfo_.getNumOfOrbitals();

    int shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = this->orbitalInfo_.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        this->shellArrayTable_[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }
}


DfHpqEngine::PGTOs DfHpqX::getPGTOs(const index_type shellIndex)
{
    DfHpqEngine::PGTOs pgtos;

    const int numOfContractions = this->orbitalInfo_.getCgtoContraction(shellIndex);
    pgtos.resize(numOfContractions);
    for (int i = 0; i < numOfContractions; ++i) {
        const DfHpqEngine::PGTO pgto(this->orbitalInfo_.getCoefficient(shellIndex, i),
                                     this->orbitalInfo_.getExponent(shellIndex, i));
        pgtos[i] = pgto;
    }
    
    return pgtos;
}


