#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfCD.h"
#include "DfEriEngine.h"
#include "TlOrbitalInfo.h"
#include "TlSymmetricMatrix.h"

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
    TlSymmetricMatrix T(dim);
    
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
                                     &T);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_, &taskList);
    }

    this->finalize(&T);
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();

    T.save("T.mat");
}


void DfCD::makeSuperMatrix_kernel(const TlOrbitalInfo& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task4>& taskList,
                                  TlSymmetricMatrix* pT)
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
            this->storeT(shellIndexP, maxStepsP,
                         shellIndexQ, maxStepsQ,
                         shellIndexR, maxStepsR,
                         shellIndexS, maxStepsS,
                         this->pEriEngines_[threadID], pT);
        }
    }
}

void DfCD::storeT(const index_type shellIndexP, const int maxStepsP,
                  const index_type shellIndexQ, const int maxStepsQ,
                  const index_type shellIndexR, const int maxStepsR,
                  const index_type shellIndexS, const int maxStepsS,
                  const DfEriEngine& engine,
                  TlSymmetricMatrix* pT)
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
                    
                    if ((indexP >= indexQ) && (maxIndexS >= indexS)) {
                        const std::size_t indexRS = this->index(indexR, indexS);
                        pT->add(indexPQ, indexRS, value);
                    }
                    ++index;
                }
            }
        }
    }
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
