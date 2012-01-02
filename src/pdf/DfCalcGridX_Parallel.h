#ifndef DFCALCGRIDX_PARALLEL_H
#define DFCALCGRIDX_PARALLEL_H

#include <set>
#include <cassert>
#include "DfCalcGridX.h"
#include "GridDataManager.h"
#include "TlCommunicate.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlPartialSymmetricMatrix.h"
#include "TlSparseVectorMatrix.h"

class DfCalcGridX_Parallel : public DfCalcGridX {
public:
    DfCalcGridX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCalcGridX_Parallel();

public:
    // template<class DfFunctionalType>
    // double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
    //                                    DfFunctionalType* pFunctional,
    //                                    TlSymmetricMatrix* pF_A);
    // template<class DfFunctionalType>
    // double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
    //                                    const TlSymmetricMatrix& P_B,
    //                                    DfFunctionalType* pFunctional,
    //                                    TlSymmetricMatrix* pF_A,
    //                                    TlSymmetricMatrix* pF_B);
public:
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pF_A,
                                       TlSymmetricMatrix* pF_B);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pF_A,
                                       TlSymmetricMatrix* pF_B);
    
    
protected:
    virtual void calcRho_LDA(const TlSymmetricMatrix& P_A);
    virtual void calcRho_LDA(const TlSymmetricMatrix& P_A,
                             const TlSymmetricMatrix& P_B);
    virtual void calcRho_GGA(const TlSymmetricMatrix& P_A);
    virtual void calcRho_GGA(const TlSymmetricMatrix& P_A,
                             const TlSymmetricMatrix& P_B);

    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlSymmetricMatrix* pF_A,
                    TlSymmetricMatrix* pF_B);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlSymmetricMatrix* pF_A,
                    TlSymmetricMatrix* pF_B);
    
    TlMatrix distributeGridMatrix(const int iteration);
    void gatherGridMatrix(const TlMatrix& gridMat);

// -----------------------------------------------------------------------------
    
public:
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       DfFunctional_LDA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B,
                                       DfFunctional_LDA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A,
                                       TlDistributeSymmetricMatrix* pF_B);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       DfFunctional_GGA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B,
                                       DfFunctional_GGA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A,
                                       TlDistributeSymmetricMatrix* pF_B);
   
protected:
    template<class SymmetricMatrixType, class DfFunctionalClass>
    double calcXCIntegForFockAndEnergy_atomParallel(const SymmetricMatrixType& P,
                                                    DfFunctionalClass* pFunctional,
                                                    SymmetricMatrixType* pF);

    template<class SymmetricMatrixType, class DfFunctionalClass>
    double calcXCIntegForFockAndEnergy_atomParallel(const SymmetricMatrixType& PA,
                                                    const SymmetricMatrixType& PB,
                                                    DfFunctionalClass* pFunctional,
                                                    SymmetricMatrixType* pFA,
                                                    SymmetricMatrixType* pFB);

    template<class SymmetricMatrixType, class DfFunctionalClass>
    double calcXCIntegForFockAndEnergy_MasterSlave(const SymmetricMatrixType& P,
                                                   DfFunctionalClass* pFunctional,
                                                   SymmetricMatrixType* pF);

    template<class SymmetricMatrixType, class DfFunctionalClass>
    double calcXCIntegForFockAndEnergy_MasterSlave(const SymmetricMatrixType& PA,
                                                   const SymmetricMatrixType& PB,
                                                   DfFunctionalClass* pFunctional,
                                                   SymmetricMatrixType* pFA,
                                                   SymmetricMatrixType* pFB);


    virtual void defineCutOffValues(const TlSymmetricMatrix& P);

    virtual void defineCutOffValues(const TlSymmetricMatrix& PA,
                                    const TlSymmetricMatrix& PB);


    virtual void backupGridData();
    virtual void flushGridData();
    
protected:

    void defineCutOffValues(const TlDistributeSymmetricMatrix& P);
    void defineCutOffValues(const TlDistributeSymmetricMatrix& PA,
                            const TlDistributeSymmetricMatrix& PB);

    void broadcastGridInfo(std::vector<int>* pGridCounts,
                           std::vector<double>* pGridInfo);
    
    void calcPhys(const TlDistributeSymmetricMatrix& P,
                  const std::vector<int>& AO_list,
                  const std::vector<double>& gridInfo,
                  std::vector<double>* pRhoAs);

    /// 分散型密度行列とグリッド情報から、各グリッド上の電子密度・勾配を算出する。
    void calcPhys(const TlDistributeSymmetricMatrix& P,
                  const std::vector<double>& gridInfo,
                  std::vector<double>* pRhoAs,
                  std::vector<double>* pGradRhoAXs,
                  std::vector<double>* pGradRhoAYs,
                  std::vector<double>* pGradRhoAZs);
    void calcPhys_Master(const TlDistributeSymmetricMatrix& P);
    void calcPhys_Slave(const TlDistributeSymmetricMatrix& P,
                        const std::vector<double>& gridInfo,
                        std::vector<double>* pRhoAs,
                        std::vector<double>* pGradRhoAXs,
                        std::vector<double>* pGradRhoAYs,
                        std::vector<double>* pGradRhoAZs);
        
    void calcPhysX(const TlMatrixObject& P,
                   const index_type* pAO_List,
                   const std::size_t AO_ListSize,
                   const std::vector<double>& gridInfo,
                   std::vector<double>* pRhoAs,
                   std::vector<double>* pGradRhoAXs,
                   std::vector<double>* pGradRhoAYs,
                   std::vector<double>* pGradRhoAZs);
    bool getQueueX(index_type* pStartShell_I,
                   index_type* pStartShell_J,
                   index_type* pEndShell_I,
                   index_type* pEndShell_J,
                   double* pProgress,
                   bool initialize = false);
    
    std::vector<int> getScaLapackLocalAOs(const TlDistributeSymmetricMatrix& P) const;

    void getRhoAtGridPoint(const TlDistributeSymmetricMatrix& P,
                           const std::vector<WFGrid>& phis,
                           double* pRhoA);
    void getRhoAtGridPoint(const TlDistributeSymmetricMatrix& P,
                           const std::vector<WFGrid>& phis,
                           const std::vector<WFGrid>& gradPhiXs,
                           const std::vector<WFGrid>& gradPhiYs,
                           const std::vector<WFGrid>& gradPhiZs,
                           double* pRhoA,
                           double* pGradRhoAX,
                           double* pGradRhoAY,
                           double* pGradRhoAZ);

    void addPreviousPhys(std::vector<double>* pRhoA);
    void addPreviousPhys(std::vector<double>* pRhoA,
                         std::vector<double>* pRhoB);
    void addPreviousPhys(std::vector<double>* pRhoA,
                         std::vector<double>* pGradRhoAX,
                         std::vector<double>* pGradRhoAY,
                         std::vector<double>* pGradRhoAZ);
    void addPreviousPhys(std::vector<double>* pRhoA,
                         std::vector<double>* pGradRhoAX,
                         std::vector<double>* pGradRhoAY,
                         std::vector<double>* pGradRhoAZ,
                         std::vector<double>* pRhoB,
                         std::vector<double>* pGradRhoBX,
                         std::vector<double>* pGradRhoBY,
                         std::vector<double>* pGradRhoBZ);

    void saveCurrentPhys(const std::vector<int>& gridCounts,
                         const std::vector<double>& rhoA);
    void saveCurrentPhys(const std::vector<int>& gridCounts,
                         const std::vector<double>& rhoA,
                         const std::vector<double>& rhoB);
    void saveCurrentPhys(const std::vector<int>& gridCounts,
                         const std::vector<double>& rhoA,
                         const std::vector<double>& gradRhoAX,
                         const std::vector<double>& gradRhoAY,
                         const std::vector<double>& gradRhoAZ);
    void saveCurrentPhys(const std::vector<int>& gridCounts,
                         const std::vector<double>& rhoA,
                         const std::vector<double>& gradRhoAX,
                         const std::vector<double>& gradRhoAY,
                         const std::vector<double>& gradRhoAZ,
                         const std::vector<double>& rhoB,
                         const std::vector<double>& gradRhoBX,
                         const std::vector<double>& gradRhoBY,
                         const std::vector<double>& gradRhoBZ);

    double buildK(const std::vector<int>& gridCounts,
                  const std::vector<double>& gridInfo,
                  const std::vector<double>& rhoAs,
                  DfFunctional_LDA* pFunctional,
                  TlSparseSymmetricMatrix* pF);
    double buildK(const std::vector<int>& gridCounts,
                  const std::vector<double>& gridInfo,
                  const std::vector<double>& rhoAs,
                  const std::vector<double>& rhoBs,
                  DfFunctional_LDA* pFunctional,
                  TlSparseSymmetricMatrix* pFA,
                  TlSparseSymmetricMatrix* pFB);
    double buildK(const std::vector<int>& gridCounts,
                  const std::vector<double>& gridInfo,
                  const std::vector<double>& rhoAs,
                  const std::vector<double>& gradRhoAXs,
                  const std::vector<double>& gradRhoAYs,
                  const std::vector<double>& gradRhoAZs,
                  DfFunctional_GGA* pFunctional,
                  TlSparseSymmetricMatrix* pF);

    double buildK_rev2(const std::vector<int>& gridCounts,
                       const std::vector<double>& gridInfo,
                       const std::vector<double>& rhoAs,
                       const std::vector<double>& gradRhoAXs,
                       const std::vector<double>& gradRhoAYs,
                       const std::vector<double>& gradRhoAZs,
                       DfFunctional_GGA* pFunctional,
                       TlDistributeSymmetricMatrix* pF);
    void buildK_Master(const std::vector<double>& gridInfo,
                       const std::vector<double>& rhoAs,
                       const std::vector<double>& gradRhoAXs,
                       const std::vector<double>& gradRhoAYs,
                       const std::vector<double>& gradRhoAZs,
                       DfFunctional_GGA* pFunctional);
    void buildK_Slave(const std::vector<int>& gridCounts,
                      const std::vector<double>& gridInfo,
                      const std::vector<double>& rhoAs,
                      const std::vector<double>& gradRhoAXs,
                      const std::vector<double>& gradRhoAYs,
                      const std::vector<double>& gradRhoAZs,
                      DfFunctional_GGA* pFunctional,
                      std::list<TlPartialSymmetricMatrix>* pF);
    
    double buildK(const std::vector<int>& gridCounts,
                  const std::vector<double>& gridInfo,
                  const std::vector<double>& rhoAs,
                  const std::vector<double>& gradRhoAXs,
                  const std::vector<double>& gradRhoAYs,
                  const std::vector<double>& gradRhoAZs,
                  const std::vector<double>& rhoBs,
                  const std::vector<double>& gradRhoBXs,
                  const std::vector<double>& gradRhoBYs,
                  const std::vector<double>& gradRhoBZs,
                  DfFunctional_GGA* pFunctional,
                  TlSparseSymmetricMatrix* pFA,
                  TlSparseSymmetricMatrix* pFB);

protected:
    double calcXCIntegForFockAndEnergy1(const TlDistributeSymmetricMatrix& P,
                                        DfFunctional_LDA* pFunctional,
                                        TlDistributeSymmetricMatrix* pF);
    double calcXCIntegForFockAndEnergy2(const TlDistributeSymmetricMatrix& P,
                                        DfFunctional_LDA* pFunctional,
                                        TlDistributeSymmetricMatrix* pF);

    /// XC項のエネルギーおよび行列要素の計算(GGA版)
    ///
    /// 密度行列の通信をバックグラウンドで行う。
    double calcXC_BG(const TlDistributeSymmetricMatrix& P,
                     DfFunctional_GGA* pFunctional,
                     TlDistributeSymmetricMatrix* pF);

    /// XC項のエネルギーおよび行列要素の計算(GGA版)
    ///
    /// 密度行列のlocal行列を利用したDivide&Conquer法により並列計算を行う。
    double calcXC_DC(const TlDistributeSymmetricMatrix& P,
                     DfFunctional_GGA* pFunctional,
                     TlDistributeSymmetricMatrix* pF);

    void calcPhys_new(const TlDistributeMatrix& P,
                      const std::vector<double>& gridInfo,
                      std::vector<double>* pRhoVals);
    void calcPhys_new(const TlDistributeMatrix& P,
                      const std::vector<double>& gridInfo,
                      std::vector<double>* pRhoVals,
                      std::vector<double>* pGradRhoXVals,
                      std::vector<double>* pGradRhoYVals,
                      std::vector<double>* pGradRhoZVals);
    TlMatrix getWaveFunctionCoef(const std::vector<index_type>& AO_indexes,
                                 const TlPosition& gridPosition);
    
    void getWaveFunctionCoef(const std::vector<index_type>& AO_indexes,
                             const TlPosition& gridPosition,
                             TlMatrix* pWF,
                             TlMatrix* pGradWF_X,
                             TlMatrix* pGradWF_Y,
                             TlMatrix* pGradWF_Z);
    
protected:
    // tag for MPI
    enum {
        TAG_REQUEST_JOB = 9001,
        TAG_ASSIGN_JOB = 9002,
        TAG_TERMINATE_SLAVE = 9003,
        TAG_TERMINATE_OK = 9004
    };


    struct DensityCalcJob {
    public:
        DensityCalcJob(index_type globalDim,
                       index_type start_row, index_type start_col,
                       index_type end_row, index_type end_col)
            : startRow(start_row), startCol(start_col), endRow(end_row), endCol(end_col),
              partialMatrix(globalDim, start_row, start_col,
                            std::max((end_row - start_row), (end_col - start_col))) {
        }
    public:
        index_type startRow;
        index_type startCol;
        index_type endRow;
        index_type endCol;
        TlPartialSymmetricMatrix partialMatrix;
    };

    virtual void saveGridMatrix(const TlMatrix& gridMat);

    // NEW IMPLIMENT -----------------------------------------------------------
protected:
    // virtual TlMatrix getGridMatrix();
    
    // bool getTask_DC(TlMatrix* pOutMat,
    //                 bool initialize = false);

    // virtual void finalizeGridMatrix(const TlMatrix& gridMat);
    
protected:
    int assignAtomRange_;
    int assignAoRange_;
    int assignJobsPerProc_;
    std::size_t densityMatrixCacheMemSize_;

    int calcMode_;
};


/// 原子毎に並列化を行う
/// 各プロセス毎の処理を担当するため、
/// エネルギーや行列の集計が必要
template<class SymmetricMatrixType, class DfFunctionalClass>
double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy_atomParallel(const SymmetricMatrixType& P,
                                                                      DfFunctionalClass* pFunctional,
                                                                      SymmetricMatrixType* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    const int nNumOfAtoms = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;

    // 原子ごとに並列化---------------------------------------------------
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    const int nRange = nNumOfAtoms;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り用
    const int nLocalStart = nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nLocalStart + nInterval), nNumOfAtoms);

    double dEnergy = DfCalcGridX::calcXCIntegForFockAndEnergy(nLocalStart, nLocalEnd,
                                                              P, pFunctional,
                                                              pF);
    // -------------------------------------------------------------------

    return dEnergy;
}


/// 原子毎に並列化を行う
/// 各プロセス毎の処理を担当するため、
/// エネルギーや行列の集計が必要
template<class SymmetricMatrixType, class DfFunctionalClass>
double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy_atomParallel(const SymmetricMatrixType& PA,
                                                                      const SymmetricMatrixType& PB,
                                                                      DfFunctionalClass* pFunctional,
                                                                      SymmetricMatrixType* pFA,
                                                                      SymmetricMatrixType* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    const int nNumOfAtoms = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;

    // 原子ごとに並列化---------------------------------------------------
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    const int nRange = nNumOfAtoms;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り用
    const int nLocalStart = nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nLocalStart + nInterval), nNumOfAtoms);

    double dEnergy = DfCalcGridX::calcXCIntegForFockAndEnergy(nLocalStart, nLocalEnd,
                                                              PA, PB, pFunctional,
                                                              pFA, pFB);
    // -------------------------------------------------------------------

    return dEnergy;
}


template<class SymmetricMatrixType, class DfFunctionalClass>
double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy_MasterSlave(const SymmetricMatrixType& P,
                                                                     DfFunctionalClass* pFunctional,
                                                                     SymmetricMatrixType* pF)
{
    assert(pFunctional != NULL);
    assert(pF != NULL);

    double energy = 0.0;
    const int numOfAtoms = this->numOfRealAtoms_;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    ////const int nRank = rComm.getRank();

    // master-slave message ID
    enum {
        SLAVE_REQUEST_JOB = 0,
        MASTER_ASSIGN_JOB = 1,
        MASTER_MAKE_SLAVE_FINISH_JOB = 2
    };

    if (rComm.isMaster() == true) {
        // for Master
        int numOfFinishJobProc = 0;
        int currentAtomIndex = 0;
        int progress = 0;
        const int assignAtomRange = this->assignAtomRange_;

        while (numOfFinishJobProc < (nProc -1)) {
            int msgFromSlave = 0;
            int slaveProc = 0;
            rComm.receiveDataFromAnySource(msgFromSlave, &slaveProc);
            //std::cout << TlUtils::format("[%d] recv %d from %d", nRank, msgFromSlave, slaveProc) << std::endl;
            assert(msgFromSlave == SLAVE_REQUEST_JOB);

            //std::cerr << TlUtils::format("[0] progress %d / %d", currentAtomIndex, numOfAtoms) << std::endl;
            if (currentAtomIndex < numOfAtoms) {
                const int msgToSlave = MASTER_ASSIGN_JOB;
                //std::cout << TlUtils::format("[%d] send %d to %d", nRank, msgToSlave, slaveProc) << std::endl;
                rComm.sendData(msgToSlave, slaveProc);

                std::vector<int> atomRange(2);
                atomRange[0] = currentAtomIndex;
                atomRange[1] = std::min(currentAtomIndex + assignAtomRange, this->numOfRealAtoms_);
                rComm.sendData(atomRange, slaveProc);
                currentAtomIndex += assignAtomRange;

                const int currentProgress = (currentAtomIndex / numOfAtoms * 100) / 10;
                if (currentProgress > progress) {
                    this->loggerTime(TlUtils::format(" %3d%% proceed.", currentProgress));
                    progress = currentProgress;
                }
                
                if (this->isDebugOut_ == true) {
                    this->logger(TlUtils::format(" [0]: send (%d, %d) to %d.\n",
                                                 atomRange[0], atomRange[1], slaveProc));
                }
            } else {
                const int msgToSlave = MASTER_MAKE_SLAVE_FINISH_JOB;
                rComm.sendData(msgToSlave, slaveProc);
                ++numOfFinishJobProc;
            }
        }

    } else {
        // for Slave
        const int masterProc = 0;

        while (true) {
            int msgToMaster = SLAVE_REQUEST_JOB;
            //std::cerr << TlUtils::format("[%d] send %d to master", nRank, msgToMaster) << std::endl;
            assert(msgToMaster == SLAVE_REQUEST_JOB);
            rComm.sendData(msgToMaster, masterProc);

            int msgFromMaster = 0;
            rComm.receiveData(msgFromMaster, masterProc);
            if (msgFromMaster == MASTER_MAKE_SLAVE_FINISH_JOB) {
                //std::cerr << TlUtils::format("[%d] break.", nRank) << std::endl;
                break;
            } else {
                std::vector<int> atomRange;
                rComm.receiveData(atomRange, masterProc);
                assert(atomRange.size() == 2);
                const int localStart = atomRange[0];
                const int localEnd = atomRange[1];
                energy += DfCalcGridX::calcXCIntegForFockAndEnergy(localStart, localEnd,
                                                                   P, pFunctional, pF);
                //std::cerr << TlUtils::format("[%d] end job. from %d to %d.", nRank, localStart, localEnd) << std::endl;
            }
        }
    }

    rComm.barrier();

    return energy;
}


template<class SymmetricMatrixType, class DfFunctionalClass>
double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy_MasterSlave(const SymmetricMatrixType& PA,
                                                                     const SymmetricMatrixType& PB,
                                                                     DfFunctionalClass* pFunctional,
                                                                     SymmetricMatrixType* pFA,
                                                                     SymmetricMatrixType* pFB)
{
    assert(pFunctional != NULL);
    assert(pFA != NULL);
    assert(pFB != NULL);

    double energy = 0.0;
    const int numOfAtoms = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    //const int nRank = rComm.getRank();

    // master-slave message ID
    enum {
        SLAVE_REQUEST_JOB,
        MASTER_ASSIGN_JOB,
        MASTER_MAKE_SLAVE_FINISH_JOB
    };

    if (rComm.isMaster() == true) {
        // for Master
        int numOfFinishJobProc = 0;
        int currentAtomIndex = 0;
        const int assignAtomRange = this->assignAtomRange_;

        while (numOfFinishJobProc < (nProc -1)) {
            int msgFromSlave = 0;
            int slaveProc = 0;
            rComm.receiveDataFromAnySource(msgFromSlave, &slaveProc);
            assert(msgFromSlave == SLAVE_REQUEST_JOB);

            if (currentAtomIndex < numOfAtoms) {
                const int msgToSlave = MASTER_ASSIGN_JOB;
                rComm.sendData(msgToSlave, slaveProc);

                std::vector<int> atomRange(2);
                atomRange[0] = currentAtomIndex;
                atomRange[1] = std::min(currentAtomIndex + assignAtomRange, this->numOfRealAtoms_);
                rComm.sendData(atomRange, slaveProc);
                currentAtomIndex += assignAtomRange;
            } else {
                const int msgToSlave = MASTER_MAKE_SLAVE_FINISH_JOB;
                rComm.sendData(msgToSlave, slaveProc);
                ++numOfFinishJobProc;
            }
        }

    } else {
        // for Slave
        const int masterProc = 0;

        while (true) {
            int msgToMaster = SLAVE_REQUEST_JOB;
            rComm.sendData(msgToMaster, masterProc);

            int msgFromMaster = 0;
            rComm.receiveData(msgFromMaster, masterProc);
            if (msgFromMaster == MASTER_MAKE_SLAVE_FINISH_JOB) {
                break;
            } else {
                std::vector<int> atomRange;
                rComm.receiveData(atomRange, masterProc);
                assert(atomRange.size() == 2);
                const int localStart = atomRange[0];
                const int localEnd = atomRange[1];
                energy += DfCalcGridX::calcXCIntegForFockAndEnergy(localStart, localEnd,
                                                                   PA, PB, pFunctional,
                                                                   pFA, pFB);
            }
        }
    }

    return energy;
}


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// template<class DfFunctionalType>
// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
//                                                          DfFunctionalType* pFunctional,
//                                                          TlSymmetricMatrix* pF_A)
// {
//     assert(pFunctional != NULL);
//     assert(pF_A != NULL);

    
//     TlCommunicate& rComm = TlCommunicate::getInstance();





    
//     TlMatrix localGridMat; // このプロセスが一時的に保持するグリッド行列
//     index_type numOfRowsOfLocalGridMat = 0;

//     double energy = 0.0;
//     TlMatrix tmpGridMat;   // while 1cycleで用いるグリッド行列
//     this->getTask_DC(NULL, true); // initialize
//     while (this->getTask_DC(&tmpGridMat)) {
//         energy += DfCalcGridX::calcXCIntegForFockAndEnergy(P_A,
//                                                            pFunctional,
//                                                            pF_A,
//                                                            &tmpGridMat);

//         // merge
//         const index_type numOfRowsOfTmpMat = tmpGridMat.getNumOfRows();
//         const index_type numOfColsOfTmpMat = tmpGridMat.getNumOfCols();
//         localGridMat.resize(numOfRowsOfLocalGridMat + numOfRowsOfTmpMat,
//                             tmpGridMat.getNumOfCols());
//         localGridMat.setBlockMatrix(numOfRowsOfLocalGridMat, 0,
//                                     tmpGridMat);
//         numOfRowsOfLocalGridMat += tmpGridMat.getNumOfRows();
//     }

//     this->finalizeGridMatrix(localGridMat);
//     rComm.allReduce_SUM(*pF_A);
//     pF_A->save("FxcA.mat");
//     rComm.allReduce_SUM(energy);

//     return energy;
// }

// template<class DfFunctionalType>
// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
//                                                          const TlSymmetricMatrix& P_B,
//                                                          DfFunctionalType* pFunctional,
//                                                          TlSymmetricMatrix* pF_A,
//                                                          TlSymmetricMatrix* pF_B)
// {
//     this->log_.debug("DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy() in.");
//     assert(pFunctional != NULL);
//     assert(pF_A != NULL);
//     assert(pF_B != NULL);
//     TlCommunicate& rComm = TlCommunicate::getInstance();
    
//     double energy = 0.0;

//     TlMatrix localGridMat; // このプロセスが一時的に保持するグリッド行列
//     index_type numOfRowsOfLocalGridMat = 0;

//     TlMatrix tmpGridMat;   // while 1cycleで用いるグリッド行列
//     this->getTask_DC(NULL, true); // initialize
//     while (this->getTask_DC(&tmpGridMat)) {
//         energy += DfCalcGridX::calcXCIntegForFockAndEnergy(P_A, P_B,
//                                                            pFunctional,
//                                                            pF_A, pF_B,
//                                                            &tmpGridMat);

//         // merge
//         const index_type numOfRowsOfTmpMat = tmpGridMat.getNumOfRows();
//         const index_type numOfColsOfTmpMat = tmpGridMat.getNumOfCols();
//         localGridMat.resize(numOfRowsOfLocalGridMat + numOfRowsOfTmpMat,
//                             tmpGridMat.getNumOfCols());
//         localGridMat.setBlockMatrix(numOfRowsOfLocalGridMat, 0,
//                                     tmpGridMat);
//         numOfRowsOfLocalGridMat += tmpGridMat.getNumOfRows();
//     }

//     this->finalizeGridMatrix(localGridMat);
//     rComm.allReduce_SUM(*pF_A);
//     pF_A->save("FxcA.mat");
//     if (pF_B != NULL) {
//         rComm.allReduce_SUM(*pF_B);
//     }
//     rComm.allReduce_SUM(energy);

//     return energy;
// }

#endif // DFCALCGRIDX_PARALLEL_H
