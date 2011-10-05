#include "DfTaskCtrl.h"
#include "TlMath.h"
#include "TlOrbitalInfo.h"
#include "TlFmt.h"

DfTaskCtrl::DfTaskCtrl(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    TlOrbitalInfo orbitalInfo((*pPdfParam)["coordinates"],
                              (*pPdfParam)["basis_sets"]);
    this->maxShellType_ = orbitalInfo.getMaxShellType();


    this->lengthScaleParameter_ = 1.0;
    if ((*pPdfParam)["length_scale_parameter"].getStr() != "") {
        this->lengthScaleParameter_ = (*pPdfParam)["length_scale_parameter"].getDouble();
    }
    
    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    

    this->cutoffEpsilon1_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon1"].getStr().empty() != true) {
        this->cutoffEpsilon1_ = (*pPdfParam)["cutoff_epsilon1"].getDouble();
    }    

    this->cutoffEpsilon2_ = this->cutoffThreshold_;
    if ((*pPdfParam)["cutoff_epsilon2"].getStr().empty() != true) {
        this->cutoffEpsilon2_ = (*pPdfParam)["cutoff_epsilon2"].getDouble();
    }    

    this->cutoffEpsilon3_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon3"].getStr().empty() != true) {
        this->cutoffEpsilon3_ = (*pPdfParam)["cutoff_epsilon3"].getDouble();
    }    
}


DfTaskCtrl::~DfTaskCtrl()
{
}


void DfTaskCtrl::clearCutoffStats(const TlOrbitalInfoObject& orbitalInfo)
{
    const int maxShellType = orbitalInfo.getMaxShellType();
    const int numOfShellPairType = maxShellType* maxShellType;
    const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;

    this->cutoffAll_E1_.clear();
    this->cutoffAlive_E1_.clear();
    this->cutoffAll_E1_.resize(numOfShellPairType, 0);
    this->cutoffAlive_E1_.resize(numOfShellPairType, 0);

    this->cutoffAll_E2_.clear();
    this->cutoffAlive_E2_.clear();
    this->cutoffAll_E2_.resize(numOfShellPairType, 0);
    this->cutoffAlive_E2_.resize(numOfShellPairType, 0);

    this->cutoffAll_schwarz_.clear();
    this->cutoffAlive_schwarz_.clear();
    this->cutoffAll_schwarz_.resize(numOfShellQuartetType, 0);
    this->cutoffAlive_schwarz_.resize(numOfShellQuartetType, 0);

}


DfTaskCtrl::ShellArrayTable DfTaskCtrl::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    const int maxShellType = orbitalInfo.getMaxShellType();
    ShellArrayTable shellArrayTable(maxShellType);
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    for (int i = 0; i < maxShellType; ++i) {
        shellArrayTable[i].reserve(maxShellIndex);
    }

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


bool DfTaskCtrl::getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const bool isCutoffByDistibution,
                          const int maxGrainSize,
                          std::vector<Task2>* pTaskList,                          
                          bool initialize)
{
    assert(pTaskList != NULL);

    const int maxShellType = orbitalInfo.getMaxShellType();
    static ShellArrayTable shellArrayTable;
    static int shellTypeP = maxShellType -1;
    static int shellTypeQ = maxShellType -1;
    static std::size_t shellArrayIndexP = 0;
    static std::size_t shellArrayIndexQ = 0;
    
    pTaskList->clear();
    pTaskList->reserve(maxGrainSize);

    if (initialize == true) {
        this->clearCutoffStats(orbitalInfo);
        shellArrayTable = this->makeShellArrayTable(orbitalInfo);
        shellTypeP = maxShellType -1;
        shellTypeQ = maxShellType -1;
        shellArrayIndexP = 0;
        shellArrayIndexQ = 0;

        return true;
    }

    int grainSize = 0;
    Task2 task;
    for ( ; shellTypeP >= 0; ) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        const size_t shellArraySizeP = shellArrayP.size();
        
        for ( ; shellTypeQ >= 0; ) {
            ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            
            for ( ; shellArrayIndexP < shellArraySizeP; ) {
                const index_type shellIndexP = shellArrayP[shellArrayIndexP];
                task.shellIndex1 = shellIndexP;

                if (isCutoffByDistibution == true) {
                    shellArrayQ =
                        this->selectShellArrayByDistribution(shellArrayQ,
                                                             shellIndexP, orbitalInfo);
                }
                ShellArray::iterator qItEnd = std::upper_bound(shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
                const std::size_t shellArraySizeQ = std::distance(shellArrayQ.begin(), qItEnd);
                
                for ( ; shellArrayIndexQ < shellArraySizeQ; ) {
                    const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
                    task.shellIndex2 = shellIndexQ;

                    pTaskList->push_back(task);
                    ++shellArrayIndexQ;
                    ++grainSize;

                    if (grainSize >= maxGrainSize) {
                        return true;
                    }
                }
                shellArrayIndexQ = 0;
                ++shellArrayIndexP;
            }
            shellArrayIndexP = 0;
            --shellTypeQ;
        }
        shellTypeQ = maxShellType -1;
        --shellTypeP;
    }

    return (pTaskList->empty() != true);
}


bool DfTaskCtrl::getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                           const TlSparseSymmetricMatrix& schwarzTable,
                           const int maxGrainSize,
                           std::vector<Task4>* pTaskList,
                           bool initialize)
{
    assert(pTaskList != NULL);

    const int maxShellType = orbitalInfo.getMaxShellType();
    static ShellArrayTable shellArrayTable;
    static ShellPairArrayTable shellPairArrayTable;
    static int shellTypeP = maxShellType -1;
    static int shellTypeQ = maxShellType -1;
    static int shellTypeR = maxShellType -1;
    static int shellTypeS = maxShellType -1;
    static std::size_t prIndex = 0;
    static std::size_t shellArrayIndexQ = 0;
    static std::size_t shellArrayIndexS = 0;
    
    pTaskList->clear();
    pTaskList->reserve(maxGrainSize);

    if (initialize == true) {
        this->clearCutoffStats(orbitalInfo);
        shellArrayTable = this->makeShellArrayTable(orbitalInfo);
        shellPairArrayTable = this->getShellPairArrayTable(orbitalInfo,
                                                           shellArrayTable);
        shellPairArrayTable = this->selectShellPairArrayTableByDensity(shellPairArrayTable,
                                                                       orbitalInfo);
        shellTypeP = maxShellType -1;
        shellTypeQ = maxShellType -1;
        shellTypeR = maxShellType -1;
        shellTypeS = maxShellType -1;
        prIndex = 0;
        shellArrayIndexQ = 0;
        shellArrayIndexS = 0;

        return true;
    }

    int grainSize = 0;
    DfTaskCtrl::Task4 task;
    for ( ; shellTypeP >= 0; ) {
        //const int maxStepsP = 2 * shellTypeP + 1;

        for ( ; shellTypeR >= 0; ) {
            //const int maxStepsR = 2 * shellTypeR + 1;
            
            const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
            const ShellPairArray& shellPairArray_PR = shellPairArrayTable[shellPairType_PR];
            const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();
            
            for ( ; prIndex < numOfShellPairArray_PR; ) {
                const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
                const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
                task.shellIndex1 = shellIndexP;
                task.shellIndex3 = shellIndexR;
                
                for ( ; shellTypeQ >= 0; ) {
                    //const int maxStepsQ = 2 * shellTypeQ + 1;
                    const ShellArray shellArrayQ =
                        this->selectShellArrayByDistribution(shellArrayTable[shellTypeQ],
                                                             shellIndexP,
                                                             orbitalInfo);
                    ShellArray::const_iterator qItEnd = std::upper_bound(shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
                    const std::size_t shellArraySizeQ = std::distance(shellArrayQ.begin(), qItEnd);
                    
                    for ( ; shellTypeS >= 0; ) {
                        //const int maxStepsS = 2 * shellTypeS + 1;
                        const ShellArray shellArrayS =
                            this->selectShellArrayByDistribution(shellArrayTable[shellTypeS],
                                                                 shellIndexR,
                                                                 orbitalInfo);

                        for ( ; shellArrayIndexQ < shellArraySizeQ; ) {
                            const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
                            task.shellIndex2 = shellIndexQ;
                            
                            const index_type maxShellIndexS = (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
                            ShellArray::const_iterator sItEnd = std::upper_bound(shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
                            const std::size_t shellArraySizeS = std::distance(shellArrayS.begin(), sItEnd);
                            for ( ; shellArrayIndexS < shellArraySizeS; ) {
                                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                                task.shellIndex4 = shellIndexS;

                                // schwarz cutoff
                                const int shellQuartetType = this->getShellQuartetType(orbitalInfo,
                                                                                       shellTypeP, shellTypeQ,
                                                                                       shellTypeR, shellTypeS);
                                const bool isAlive = this->isAliveBySchwarzCutoff(shellIndexP, shellIndexQ,
                                                                                  shellIndexR, shellIndexS,
                                                                                  shellQuartetType,
                                                                                  schwarzTable,
                                                                                  this->cutoffThreshold_);
                                ++shellArrayIndexS;
                                if (isAlive == true) {
                                    pTaskList->push_back(task);
                                    ++grainSize;
                                    
                                    if (grainSize >= maxGrainSize) {
                                        return true;
                                    }
                                }
                            }
                            shellArrayIndexS = 0;
                            ++shellArrayIndexQ;
                        }
                        shellArrayIndexQ = 0;
                        --shellTypeS;
                    }
                    shellTypeS = maxShellType -1;
                    --shellTypeQ;
                }
                shellTypeQ = maxShellType -1;
                ++prIndex;
            }
            prIndex = 0;
            --shellTypeR;
        }
        shellTypeR = maxShellType -1;
        --shellTypeP;
    }

    if  (pTaskList->empty() != true) {
        return true;
    }

    return false;
}

bool DfTaskCtrl::getQueue4_K(const TlOrbitalInfoObject& orbitalInfo,
                             const TlSparseSymmetricMatrix& schwarzTable,
                             const TlMatrixObject& P,
                             const std::vector<TlMatrixObject::index_type>& rowIndexes,
                             const std::vector<TlMatrixObject::index_type>& colIndexes,
                             const int maxGrainSize,
                             std::vector<Task4>* pTaskList,
                             bool initialize)
{
    assert(pTaskList != NULL);

    const int maxShellType = orbitalInfo.getMaxShellType();
    static ShellArrayTable shellArrayTable;
    static ShellPairArrayTable shellPairArrayTable;
    static int shellTypeP = maxShellType -1;
    static int shellTypeQ = maxShellType -1;
    static int shellTypeR = maxShellType -1;
    static int shellTypeS = maxShellType -1;
    static std::size_t prIndex = 0;
    static std::size_t shellArrayIndexQ = 0;
    static std::size_t shellArrayIndexS = 0;
    
    pTaskList->clear();
    pTaskList->reserve(maxGrainSize);

    if (initialize == true) {
        this->clearCutoffStats(orbitalInfo);
        shellArrayTable = this->makeShellArrayTable(orbitalInfo);
        shellPairArrayTable = this->getShellPairArrayTable(orbitalInfo,
                                                           shellArrayTable);
        shellPairArrayTable = this->selectShellPairArrayTableByDensity(shellPairArrayTable,
                                                                       orbitalInfo);
        shellTypeP = maxShellType -1;
        shellTypeQ = maxShellType -1;
        shellTypeR = maxShellType -1;
        shellTypeS = maxShellType -1;
        prIndex = 0;
        shellArrayIndexQ = 0;
        shellArrayIndexS = 0;

        return true;
    }

    //int grainSize = 0;
    DfTaskCtrl::Task4 task;

    const index_type numOfRowIndexes = rowIndexes.size();
    const index_type numOfColIndexes = colIndexes.size();
    for (index_type index_p = 0; index_p < numOfRowIndexes; ++index_p) {
        const index_type shellIndexP = rowIndexes[index_p];
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = shellTypeP * 2 + 1;

        task.shellIndex1 = shellIndexP;

        for (index_type index_r = 0; index_r < numOfColIndexes; ++index_r) {
            const index_type shellIndexR = colIndexes[index_r];
            const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
            const int maxStepsR = shellTypeR * 2 + 1;

            const double max_P_pr = std::max<double>(1.0E-20,
                                                     this->getMaxValue(P,
                                                                       shellIndexP, maxStepsP,
                                                                       shellIndexR, maxStepsR));
            task.shellIndex3 = shellIndexR;

            for (index_type index_q = 0; index_q < numOfRowIndexes; ++index_q) {
                const index_type shellIndexQ = rowIndexes[index_q];
                if (shellIndexQ > shellIndexP) {
                    continue;
                }
                const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
                const int maxStepsQ = shellTypeQ * 2 + 1;

                const double max_P_qr = this->getMaxValue(P,
                                                          shellIndexQ, maxStepsQ,
                                                          shellIndexR, maxStepsR);
                const double max_P_pqr = std::max(max_P_pr, max_P_qr);
                task.shellIndex2 = shellIndexQ;
                
                const index_type maxShellIndexS = (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
                for (index_type index_s = 0; index_s < numOfColIndexes; ++index_s) {
                    const index_type shellIndexS = colIndexes[index_s];
                    if (shellIndexS > maxShellIndexS) {
                        continue;
                    }
                    const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
                    const int maxStepsS = shellTypeS * 2 + 1;

                    const double max_P_ps = this->getMaxValue(P,
                                                              shellIndexP, maxStepsP,
                                                              shellIndexS, maxStepsS);
                    const double max_P_qs = this->getMaxValue(P,
                                                              shellIndexQ, maxStepsQ,
                                                              shellIndexS, maxStepsS);
                    // schwarz cutoff
                    const double max_P_pqs = std::max(max_P_ps, max_P_qs);
                    const double max_P_pqrs = std::max(max_P_pqr, max_P_pqs);
                    const double cutoffThreshold = this->cutoffThreshold_ / max_P_pqrs;
                    const int shellQuartetType = this->getShellQuartetType(orbitalInfo,
                                                                           shellTypeP, shellTypeQ,
                                                                           shellTypeR, shellTypeS);
                    const bool isAlive = this->isAliveBySchwarzCutoff(shellIndexP, shellIndexQ,
                                                                      shellIndexR, shellIndexS,
                                                                      shellQuartetType,
                                                                      schwarzTable,
                                                                      cutoffThreshold);
                    if (isAlive == true) {
                        task.shellIndex4 = shellIndexS;
                        pTaskList->push_back(task);
                    }
                }
            }
        }
    }
    
    if  (pTaskList->empty() != true) {
        return true;
    }

    return false;
}


bool DfTaskCtrl::getQueue4_K0(const TlOrbitalInfoObject& orbitalInfo,
                              const TlSparseSymmetricMatrix& schwarzTable,
                              const TlMatrixObject& P,
                              const std::vector<TlMatrixObject::index_type>& rowIndexes,
                              const std::vector<TlMatrixObject::index_type>& colIndexes,
                              const int maxGrainSize,
                              std::vector<Task4>* pTaskList,
                              bool initialize)
{
    assert(pTaskList != NULL);

    const int maxShellType = orbitalInfo.getMaxShellType();
    static ShellArrayTable shellArrayTable;
    static ShellPairArrayTable shellPairArrayTable;
    static int shellTypeP = maxShellType -1;
    static int shellTypeQ = maxShellType -1;
    static int shellTypeR = maxShellType -1;
    static int shellTypeS = maxShellType -1;
    static std::size_t prIndex = 0;
    static std::size_t shellArrayIndexQ = 0;
    static std::size_t shellArrayIndexS = 0;
    
    pTaskList->clear();
    pTaskList->reserve(maxGrainSize);

    if (initialize == true) {
        this->clearCutoffStats(orbitalInfo);
        shellArrayTable = this->makeShellArrayTable(orbitalInfo);
        shellPairArrayTable = this->getShellPairArrayTable(orbitalInfo,
                                                           shellArrayTable);
        shellPairArrayTable = this->selectShellPairArrayTableByDensity(shellPairArrayTable,
                                                                       orbitalInfo);
        shellTypeP = maxShellType -1;
        shellTypeQ = maxShellType -1;
        shellTypeR = maxShellType -1;
        shellTypeS = maxShellType -1;
        prIndex = 0;
        shellArrayIndexQ = 0;
        shellArrayIndexS = 0;

        return true;
    }

    int grainSize = 0;
    DfTaskCtrl::Task4 task;
    for ( ; shellTypeP >= 0; ) {
        const int maxStepsP = 2 * shellTypeP + 1;

        for ( ; shellTypeR >= 0; ) {
            const int maxStepsR = 2 * shellTypeR + 1;
            
            const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
            const ShellPairArray& shellPairArray_PR = shellPairArrayTable[shellPairType_PR];
            const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();
            
            for ( ; prIndex < numOfShellPairArray_PR; ) {
                const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
                const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
                task.shellIndex1 = shellIndexP;
                task.shellIndex3 = shellIndexR;

                // if ((std::binary_search(rowIndexes.begin(), rowIndexes.end(), shellIndexP) != true) ||
                //     (std::binary_search(colIndexes.begin(), colIndexes.end(), shellIndexR) != true)) {
                //     continue;
                // }
                const double max_P_pr = std::max<double>(1.0E-20,
                                                         this->getMaxValue(P,
                                                                           shellIndexP, maxStepsP,
                                                                           shellIndexR, maxStepsR));
                
                for ( ; shellTypeQ >= 0; ) {
                    const int maxStepsQ = 2 * shellTypeQ + 1;
                    const ShellArray shellArrayQ =
                        this->selectShellArrayByDistribution(shellArrayTable[shellTypeQ],
                                                             shellIndexP,
                                                             orbitalInfo);
                    ShellArray::const_iterator qItEnd = std::upper_bound(shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
                    const std::size_t shellArraySizeQ = std::distance(shellArrayQ.begin(), qItEnd);
                    
                    for ( ; shellTypeS >= 0; ) {
                        const int maxStepsS = 2 * shellTypeS + 1;
                        const ShellArray shellArrayS =
                            this->selectShellArrayByDistribution(shellArrayTable[shellTypeS],
                                                                 shellIndexR,
                                                                 orbitalInfo);

                        for ( ; shellArrayIndexQ < shellArraySizeQ; ) {
                            const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
                            task.shellIndex2 = shellIndexQ;
                            
                            // if (std::binary_search(rowIndexes.begin(), rowIndexes.end(), shellIndexQ) != true) {
                            //     continue;
                            // }
                            const double max_P_qr = this->getMaxValue(P,
                                                                      shellIndexQ, maxStepsQ,
                                                                      shellIndexR, maxStepsR);
                            const double max_P_pqr = std::max(max_P_pr, max_P_qr);
                            
                            const index_type maxShellIndexS = (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
                            ShellArray::const_iterator sItEnd = std::upper_bound(shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
                            const std::size_t shellArraySizeS = std::distance(shellArrayS.begin(), sItEnd);
                            for ( ; shellArrayIndexS < shellArraySizeS; ) {
                                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                                task.shellIndex4 = shellIndexS;

                                // if (std::binary_search(colIndexes.begin(), colIndexes.end(), shellIndexS) != true) {
                                //     continue;
                                // }
                                const double max_P_ps = this->getMaxValue(P,
                                                                          shellIndexP, maxStepsP,
                                                                          shellIndexS, maxStepsS);
                                const double max_P_qs = this->getMaxValue(P,
                                                                          shellIndexQ, maxStepsQ,
                                                                          shellIndexS, maxStepsS);
                                
                                // schwarz cutoff
                                const double max_P_pqs = std::max(max_P_ps, max_P_qs);
                                const double max_P_pqrs = std::max(max_P_pqr, max_P_pqs);
                                const double cutoffThreshold = this->cutoffThreshold_ / max_P_pqrs;
                                const int shellQuartetType = this->getShellQuartetType(orbitalInfo,
                                                                                       shellTypeP, shellTypeQ,
                                                                                       shellTypeR, shellTypeS);
                                const bool isAlive = this->isAliveBySchwarzCutoff(shellIndexP, shellIndexQ,
                                                                                  shellIndexR, shellIndexS,
                                                                                  shellQuartetType,
                                                                                  schwarzTable,
                                                                                  cutoffThreshold);
                                ++shellArrayIndexS;
                                if (isAlive == true) {
                                    pTaskList->push_back(task);
                                    ++grainSize;
                                    
                                    if (grainSize >= maxGrainSize) {
                                        return true;
                                    }
                                }
                            }
                            shellArrayIndexS = 0;
                            ++shellArrayIndexQ;
                        }
                        shellArrayIndexQ = 0;
                        --shellTypeS;
                    }
                    shellTypeS = maxShellType -1;
                    --shellTypeQ;
                }
                shellTypeQ = maxShellType -1;
                ++prIndex;
            }
            prIndex = 0;
            --shellTypeR;
        }
        shellTypeR = maxShellType -1;
        --shellTypeP;
    }

    if  (pTaskList->empty() != true) {
        return true;
    }

    return false;
}

double DfTaskCtrl::getMaxValue(const TlMatrixObject& P,
                               const TlMatrixObject::index_type row, const TlMatrixObject::index_type d_row, 
                               const TlMatrixObject::index_type col, const TlMatrixObject::index_type d_col)
{
    const TlMatrixObject::index_type maxRow = row + d_row;
    const TlMatrixObject::index_type maxCol = col + d_col;
    double answer = 0.0;
    for (TlMatrixObject::index_type r = row; r < maxRow; ++r) {
        for (TlMatrixObject::index_type c = col; c < maxCol; ++c) {
            answer = std::max(answer, std::fabs(P.getLocal(r, c)));
        }
    }

    return answer;
}

bool DfTaskCtrl::getQueue_Force4(const TlOrbitalInfoObject& orbitalInfo,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const int maxGrainSize,
                                 std::vector<Task4>* pTaskList,
                                 bool initialize)
{
    assert(pTaskList != NULL);

    const int maxShellType = orbitalInfo.getMaxShellType();
    static ShellArrayTable shellArrayTable;
    static ShellPairArrayTable shellPairArrayTable;
    static int shellTypeP = maxShellType -1;
    static int shellTypeQ = maxShellType -1;
    static int shellTypeR = maxShellType -1;
    static int shellTypeS = maxShellType -1;
    static std::size_t prIndex = 0;
    static std::size_t shellArrayIndexQ = 0;
    static std::size_t shellArrayIndexS = 0;
    
    pTaskList->clear();
    pTaskList->reserve(maxGrainSize);

    if (initialize == true) {
        this->clearCutoffStats(orbitalInfo);
        shellArrayTable = this->makeShellArrayTable(orbitalInfo);
        shellPairArrayTable = this->getShellPairArrayTable(orbitalInfo,
                                                           shellArrayTable);
        shellPairArrayTable = this->selectShellPairArrayTableByDensity(shellPairArrayTable,
                                                                       orbitalInfo);
        shellTypeP = maxShellType -1;
        shellTypeQ = maxShellType -1;
        shellTypeR = maxShellType -1;
        shellTypeS = maxShellType -1;
        prIndex = 0;
        shellArrayIndexQ = 0;
        shellArrayIndexS = 0;

        return true;
    }

    int grainSize = 0;
    DfTaskCtrl::Task4 task;
    for ( ; shellTypeP >= 0; ) {

        for ( ; shellTypeR >= 0; ) {
            const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
            const ShellPairArray& shellPairArray_PR = shellPairArrayTable[shellPairType_PR];
            const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();
            
            for ( ; prIndex < numOfShellPairArray_PR; ) {
                const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
                const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
                task.shellIndex1 = shellIndexP;
                task.shellIndex3 = shellIndexR;
                
                for ( ; shellTypeQ >= 0; ) {
                    //const int maxStepsQ = 2 * shellTypeQ + 1;
                    const ShellArray shellArrayQ =
                        this->selectShellArrayByDistribution(shellArrayTable[shellTypeQ],
                                                             shellIndexP,
                                                             orbitalInfo);
                    ShellArray::const_iterator qItEnd = std::upper_bound(shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
                    const index_type shellArraySizeQ = std::distance(shellArrayQ.begin(), qItEnd);
                    
                    for ( ; shellTypeS >= 0; ) {
                        //const int maxStepsS = 2 * shellTypeS + 1;
                        const ShellArray shellArrayS =
                            this->selectShellArrayByDistribution(shellArrayTable[shellTypeS],
                                                                 shellIndexR,
                                                                 orbitalInfo);
                        ShellArray::const_iterator sItEnd = std::upper_bound(shellArrayS.begin(), shellArrayS.end(), shellIndexR);
                        const std::size_t shellArraySizeS = std::distance(shellArrayS.begin(), sItEnd);

                        for ( ; shellArrayIndexQ < shellArraySizeQ; ) {
                            const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
                            task.shellIndex2 = shellIndexQ;
                            
                            for ( ; shellArrayIndexS < shellArraySizeS; ) {
                                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                                task.shellIndex4 = shellIndexS;

                                // schwarz cutoff
                                const int shellQuartetType = this->getShellQuartetType(orbitalInfo,
                                                                                       shellTypeP, shellTypeQ,
                                                                                       shellTypeR, shellTypeS);
                                const bool isAlive = this->isAliveBySchwarzCutoff(shellIndexP, shellIndexQ,
                                                                                  shellIndexP, shellIndexQ,
                                                                                  shellQuartetType,
                                                                                  schwarzTable,
                                                                                  this->cutoffThreshold_);
                                ++shellArrayIndexS;
                                if (isAlive == true) {
                                    pTaskList->push_back(task);
                                    ++grainSize;
                                    
                                    if (grainSize >= maxGrainSize) {
                                        return true;
                                    }
                                }
                            }
                            shellArrayIndexS = 0;
                            ++shellArrayIndexQ;
                        }
                        shellArrayIndexQ = 0;
                        --shellTypeS;
                    }
                    shellTypeS = maxShellType -1;
                    --shellTypeQ;
                }
                shellTypeQ = maxShellType -1;
                ++prIndex;
            }
            prIndex = 0;
            --shellTypeR;
        }
        shellTypeR = maxShellType -1;
        --shellTypeP;
    }

    if  (pTaskList->empty() != true) {
        return true;
    }

    return false;
}


// J. Chem. Phys.,105,2726 (1996)
// eq.32
DfTaskCtrl::ShellArray DfTaskCtrl::selectShellArrayByDistribution(const ShellArray& inShellArray,
                                                                  const index_type companionShellIndex,
                                                                  const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArray answer;
    answer.reserve(inShellArray.size());
    
    const TlPosition posB = orbitalInfo.getPosition(companionShellIndex);
    const int shellTypeB = orbitalInfo.getShellType(companionShellIndex);
    // orbitalInfoのPGTOリストは指数が小さい順にソートされているため、
    // 最初(index=0)の指数のみをチェックすれば良い。
    const double exponentB = orbitalInfo.getExponent(companionShellIndex, 0);
    
    // check
    const int maxShellType = orbitalInfo.getMaxShellType();
    static const double INV_EQ32_COEF = 1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
    const double threshold = this->cutoffEpsilon2_ * INV_EQ32_COEF;
    ShellArray::const_iterator itEnd = inShellArray.end();
    for (ShellArray::const_iterator it = inShellArray.begin(); it != itEnd; ++it) {
        const int shellPairType = orbitalInfo.getShellType(*it) * maxShellType + shellTypeB;
        const double distance2 = posB.squareDistanceFrom(orbitalInfo.getPosition(*it));
        const double exponentA = orbitalInfo.getExponent(*it, 0);

        const double zetaP = exponentA + exponentB;
        const double zeta = exponentA * exponentB / zetaP;

        const double exponent = - zeta * distance2;
        const double coef = 1.0 / (std::pow(zetaP, 1.25));

        if (coef * std::exp(exponent) >= threshold) {
            answer.push_back(*it);

#pragma omp atomic
            ++(this->cutoffAlive_E2_[shellPairType]);
        }

#pragma omp atomic
        ++(this->cutoffAll_E2_[shellPairType]);
    }

    // swap technique
    ShellArray(answer).swap(answer);
    
    return answer;
}


DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::getShellPairArrayTable(const TlOrbitalInfoObject& orbitalInfo,
                                                                   const ShellArrayTable& shellArrayTable)
{
    const int maxShellType = orbitalInfo.getMaxShellType();
    ShellPairArrayTable shellPairArrayTable(maxShellType * maxShellType);

    for (int shellTypeP = maxShellType -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = maxShellType -1; shellTypeR >= 0; --shellTypeR) {
            const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
            ShellArray::const_iterator rItEnd = shellArrayR.end();

            const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
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


// J. Chem. Phys.,105,2726 (1996)
// eq.31
// 1/r cutoff
DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::selectShellPairArrayTableByDensity(
    const ShellPairArrayTable& inShellPairArrayTable,
    const TlOrbitalInfoObject& orbitalInfo)
{
    const int maxShellType = this->maxShellType_;
    assert(inShellPairArrayTable.size() == (maxShellType * maxShellType));
    
    const double cutoffThreshold = this->cutoffEpsilon1_;
    const double CONTRIBUTE_COEF = 2.0 * std::pow(TlMath::PI(), 2.5);
    
    TlFmt& FmT = TlFmt::getInstance();
    static const double coef[6] = {
        -0.017450254,
         0.132520568,
        -0.047915444,
         0.792267596,
        -0.583721015,
         0.697593555
    };
    // l = 1.0のときのgamma値
    static const double gamma1[6] = {
        0.03,
        0.081722098,
        0.222616709,
        0.606423483,
        1.651939972,
        4.5
    };
    std::vector<double> gamma(6);
    {
        const double inv_lprime = 1.0 / this->lengthScaleParameter_;
        const double ll = inv_lprime * inv_lprime;
        for (int i = 0; i < 6; ++i) {
            gamma[i] = gamma1[i] * ll;
        }
    }

    const int maxShellPairType = maxShellType * maxShellType;
    ShellPairArrayTable answer(maxShellPairType);
    for (int shellPairType = 0; shellPairType < maxShellPairType; ++shellPairType) {
        const ShellPairArray& shellPairArray = inShellPairArrayTable[shellPairType];
        const std::size_t shellPairArraySize = shellPairArray.size();
        ShellPairArray tmp;
        tmp.reserve(shellPairArraySize);
        
        for (std::size_t shellPairIndex = 0; shellPairIndex < shellPairArraySize; ++shellPairIndex) {
            const index_type shellIndexA = shellPairArray[shellPairIndex].shellIndex1;
            const index_type shellIndexB = shellPairArray[shellPairIndex].shellIndex2;

            const TlPosition posA = orbitalInfo.getPosition(shellIndexA);
            const int numOfContractionsA = orbitalInfo.getCgtoContraction(shellIndexA);
            const TlPosition posB = orbitalInfo.getPosition(shellIndexB);
            const int numOfContractionsB = orbitalInfo.getCgtoContraction(shellIndexB);
            const double AB2 = posB.squareDistanceFrom(posA);
        
            double judge = 0.0;
            for (int pgtoIndexA = 0; pgtoIndexA < numOfContractionsA; ++pgtoIndexA) {
                const double coefA = orbitalInfo.getCoefficient(shellIndexA, pgtoIndexA);
                const double zetaA = orbitalInfo.getExponent(shellIndexA, pgtoIndexA);
                
                for (int pgtoIndexB = 0; pgtoIndexB < numOfContractionsB; ++pgtoIndexB) {
                    const double coefB = orbitalInfo.getCoefficient(shellIndexB, pgtoIndexB);
                    const double zetaB = orbitalInfo.getExponent(shellIndexB, pgtoIndexB);
                    
                    const double coefAB = coefA * coefB;
                    const double zetaAB = zetaA * zetaB;
                    const double zetaA_B = zetaA + zetaB;
                    for (int i = 0; i < 6; ++i) {
                        const double zetaAgamma = zetaA * gamma[i];
                        const double zetaBgamma = zetaB * gamma[i];
                        const double param = zetaAB + zetaAgamma + zetaBgamma;
                        const double term1 = CONTRIBUTE_COEF / (std::sqrt(zetaA_B) * param);
                        const double term2 = std::exp(-zetaAB * gamma[i] * AB2 / param);
                        const double T = zetaAB * zetaAB * AB2 / (zetaA_B * param);
                        double term3 = 0.0;
                        FmT.getFmT(0, T, &term3);
                        judge += coefAB * coef[i] * term1 * term2 * term3;
                    }
                }
            }

            if (std::fabs(judge) > cutoffThreshold) {
                tmp.push_back(shellPairArray[shellPairIndex]);

#pragma omp atomic
                ++(this->cutoffAlive_E1_[shellPairType]);
            }

#pragma omp atomic
            ++(this->cutoffAll_E1_[shellPairType]);
        }

        // swap technique
        ShellPairArray(tmp).swap(tmp);

        answer[shellPairType] = tmp;
    }
    
    return answer;
}


TlSparseSymmetricMatrix DfTaskCtrl::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo,
                                                     DfEriEngine* pEngine)
{
    assert(pEngine != NULL);
    
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
            pEngine->calc(queryPQ, queryPQ, PQ, PQ);

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


int DfTaskCtrl::getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo,
                                    const int shellTypeP,
                                    const int shellTypeQ,
                                    const int shellTypeR,
                                    const int shellTypeS)
{
    const int maxShellType = orbitalInfo.getMaxShellType();
    // schwarz cutoff
    const int shellQuartetType =
        ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;

    return shellQuartetType;
}


bool DfTaskCtrl::isAliveBySchwarzCutoff(const index_type shellIndexP,
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
        ++(this->cutoffAlive_schwarz_[shellQuartetType]);
    }

#pragma omp atomic
    ++(this->cutoffAll_schwarz_[shellQuartetType]);

    return answer;
}


void DfTaskCtrl::cutoffReport()
{
    const int maxShellType = this->maxShellType_;
    static const char typeStr2[][3] = {
        "SS", "SP", "SD",
        "PS", "PP", "PD",
        "DS", "DP", "DD"
    };
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
    
    // cutoff report for Epsilon1
    bool hasCutoff1 = false;
    for (int shellTypeA = 0; ((hasCutoff1 == false) && (shellTypeA < maxShellType)); ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellPairType = shellTypeA * maxShellType + shellTypeB;
            if (this->cutoffAll_E1_[shellPairType] != 0) {
                hasCutoff1 = true;
                break;
            }
        }
    }
    if (hasCutoff1 == true) {
        this->logger(TlUtils::format(" density cutoff report\n"));
        this->logger(TlUtils::format(" type: alive / all (ratio )\n"));
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellPairType = shellTypeA * maxShellType + shellTypeB;
                
                if (this->cutoffAll_E1_[shellPairType] > 0) {
                    const double ratio = (double)this->cutoffAlive_E1_[shellPairType]
                        / (double)this->cutoffAll_E1_[shellPairType]
                        * 100.0;
                    this->logger(TlUtils::format(" %2s: %12ld / %12ld (%6.2f%%)\n",
                                                 typeStr2[shellPairType],
                                                 this->cutoffAlive_E1_[shellPairType],
                                                 this->cutoffAll_E1_[shellPairType],
                                                 ratio));
                }
            }
        }
        this->logger("\n");
    }

    // cutoff report for Epsilon2
    bool hasCutoff2 = false;
    for (int shellTypeA = 0; ((hasCutoff2 == false) && (shellTypeA < maxShellType)); ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellPairType = shellTypeA * maxShellType + shellTypeB;
            if (this->cutoffAll_E2_[shellPairType] != 0) {
                hasCutoff2 = true;
                break;
            }
        }
    }
    if (hasCutoff2 == true) {
        this->logger(TlUtils::format(" distribute cutoff report\n"));
        this->logger(TlUtils::format(" type: alive / all (ratio )\n"));
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellPairType = shellTypeA * maxShellType + shellTypeB;
                
                if (this->cutoffAll_E2_[shellPairType] > 0) {
                    const double ratio = (double)this->cutoffAlive_E2_[shellPairType]
                        / (double)this->cutoffAll_E2_[shellPairType]
                        * 100.0;
                    this->logger(TlUtils::format(" %2s: %12ld / %12ld (%6.2f%%)\n",
                                                 typeStr2[shellPairType],
                                                 this->cutoffAlive_E2_[shellPairType],
                                                 this->cutoffAll_E2_[shellPairType],
                                                 ratio));
                }
            }
        }
        this->logger("\n");
    }

    // cutoff report for schwarz
    bool hasCutoffSchwarz = false;
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
            for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                    const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                    if (this->cutoffAll_schwarz_[shellTypeABCD] != 0) {
                        hasCutoffSchwarz = true;
                        break;
                    }
                }
            }
        }
    }
    if (hasCutoffSchwarz == true) {
        this->logger(TlUtils::format(" schwarz cutoff report\n"));
        this->logger(TlUtils::format(" type: alive / all (ratio )\n"));
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
                for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                    const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                    for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                        const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                        
                        if (this->cutoffAll_schwarz_[shellTypeABCD] > 0) {
                            const double ratio = (double)this->cutoffAlive_schwarz_[shellTypeABCD]
                                / (double)this->cutoffAll_schwarz_[shellTypeABCD]
                                * 100.0;
                            this->logger(TlUtils::format(" %4s: %12ld / %12ld (%6.2f%%)\n",
                                                         typeStr4[shellTypeABCD],
                                                         this->cutoffAlive_schwarz_[shellTypeABCD],
                                                         this->cutoffAll_schwarz_[shellTypeABCD],
                                                         ratio));
                        }
                    }
                }
            }
        }
    }
}
