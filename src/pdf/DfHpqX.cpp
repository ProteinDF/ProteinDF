// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfHpqX.h"
#include "TlOrbitalInfoObject.h"

DfHpqX::DfHpqX(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["coordinates"],
                   (*pPdfParam)["basis_sets"]) {
    this->makeShellArrayTable();
}


DfHpqX::~DfHpqX()
{
}


void DfHpqX::createEngines()
{
    this->pEngines_ = new DfHpqEngine[this->numOfThreads_];
}


void DfHpqX::destroyEngines()
{
    if (this->pEngines_ != NULL) {
        delete[] this->pEngines_;
        this->pEngines_ = NULL;
    }
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

void DfHpqX::finalize(std::vector<double>* pValues)
{
    // do nothing
}

void DfHpqX::getHpq(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2)
{
    const int numOfAOs = this->m_nNumOfAOs;
    
    // make coordinates
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    const int numOfDummyAtoms = flGeom.getNumOfDummyAtoms();
    const int numOfRealAtoms = numOfAtoms - numOfDummyAtoms;
    std::vector<TlAtom> Cs(numOfRealAtoms);
    std::vector<TlAtom> Xs(numOfDummyAtoms);
    std::size_t realAtomIndex = 0;
    std::size_t dummyAtomIndex = 0;
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string atomName = flGeom.getAtomSymbol(i);
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
    
    pHpq->resize(numOfAOs);
    pHpq2->resize(numOfAOs);

    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                        true,
                                        this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getHpq_part(this->orbitalInfo_,
                          taskList,
                          Cs, Xs,
                          pHpq, pHpq2);
        
        hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                       true,
                                       this->grainSize_, &taskList);
    } 

    if (this->chargeExtrapolateNumber_ > 0) {
        *pHpq2 /= static_cast<double>(this->chargeExtrapolateNumber_);
    }

    this->finalize(pHpq, pHpq2);

    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
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

    this->createEngines();
    
    const int maxShellType = TlOrbitalInfoObject::getMaxShellType();
    for (int shellTypeP = maxShellType -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray shellArrayP = this->shellArrayTable_[shellTypeP];
        const std::size_t shellArraySizeP = shellArrayP.size();
        
        for (int shellTypeQ = maxShellType -1; shellTypeQ >= 0; --shellTypeQ) {
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

    this->destroyEngines();
}


void DfHpqX::getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                               const int shellTypeP, const int shellTypeQ,
                               const index_type shellIndexP,
                               const ShellArray& shellArrayQ,
                               const TlMatrixObject& P,
                               TlMatrix* pForce)
{
    static const int BUFFER_SIZE_NUC = 3 * 5 * 5 * 5; // (xyz) * 5d * 5d * 5d
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);

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
                const std::string atomName = flGeom.getAtomSymbol(atomIndexC);
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
    this->shellArrayTable_.resize(TlOrbitalInfoObject::getMaxShellType());
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
    const int numOfContractions = this->orbitalInfo_.getCgtoContraction(shellIndex);
    DfHpqEngine::PGTOs pgtos(numOfContractions);

    for (int i = 0; i < numOfContractions; ++i) {
        const DfHpqEngine::PGTO pgto(this->orbitalInfo_.getCoefficient(shellIndex, i),
                                     this->orbitalInfo_.getExponent(shellIndex, i));
        pgtos[i] = pgto;
    }
    
    return pgtos;
}

std::vector<double> DfHpqX::getESP(const TlMatrixObject& P,
                                   const std::vector<TlPosition>& grids)
{
    std::vector<double> values(grids.size());

    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                        true,
                                        this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getESP_part(this->orbitalInfo_,
                          taskList,
                          P,
                          grids,
                          &values);
        
        hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                       true,
                                       this->grainSize_, &taskList);
    } 

    this->finalize(&values);

    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();

    return values;
}

void DfHpqX::getESP_part(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         const TlMatrixObject& P,
                         const std::vector<TlPosition>& grids,
                         std::vector<double>* pValues)
{
    assert(pValues != NULL);
    const std::size_t numOfGrids = grids.size();
    pValues->resize(numOfGrids);

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

            if (shellIndexP != shellIndexQ) {
                for (std::size_t r = 0; r < numOfGrids; ++r) {
                    const TlPosition& posR = grids[r];
                    this->pEngines_[threadID].calc(query, posP, posQ, pgtosP, pgtosQ, posR);

                    double esp = 0.0;
                    int index = 0;
                    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                        const index_type globalShellIndexP = shellIndexP + stepP;
                        for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                            const index_type globalShellIndexQ = shellIndexQ + stepQ;
                            
                            esp += 2.0 * P.get(globalShellIndexP, globalShellIndexQ) * this->pEngines_[threadID].WORK_NUC[index];
                            ++index;
                        }
                    }

#pragma omp critical(DfHpqX__getESP_part1)
                    {
                        (*pValues)[r] += esp;
                    }
                }

            } else {
                for (std::size_t r = 0; r < numOfGrids; ++r) {
                    const TlPosition& posR = grids[r];
                    this->pEngines_[threadID].calc(query, posP, posQ, pgtosP, pgtosQ, posR);

                    double esp = 0.0;
                    int index = 0;
                    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                        const index_type globalShellIndexP = shellIndexP + stepP;                    
                        for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                            const index_type globalShellIndexQ = shellIndexQ + stepQ;
                            
                            if (globalShellIndexP >= globalShellIndexQ) {
                                double coef = (globalShellIndexP != globalShellIndexQ) ? 2.0 : 1.0;
                                esp += coef * P.get(globalShellIndexP, globalShellIndexQ) * this->pEngines_[threadID].WORK_NUC[index];
                            }
                            ++index;
                        }
                    }

#pragma omp critical(DfHpqX__getESP_part2)
                    {
                        (*pValues)[r] += esp;
                    }
                }
            }
        }
    }
}

