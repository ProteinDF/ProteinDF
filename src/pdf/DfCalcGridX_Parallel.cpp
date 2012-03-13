#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <unistd.h>
#include <cassert>
#include <vector>
#include <map>
#include <queue>
#include <cassert>

#include "DfCalcGridX_Parallel.h"
#include "GridDataManager.h"
#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlFileMatrix.h"

//#define USE_FILE_MATRIX
#define JOB_PROTOCOL_SIZE (4)

DfCalcGridX_Parallel::DfCalcGridX_Parallel(TlSerializeData* pPdfParam)
    : DfCalcGridX(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;

    this->assignAtomRange_ = 1;
    if (pdfParam["xc-ms-atom-range"].getStr().empty() != true) {
        this->assignAtomRange_ = pdfParam["xc-ms-atom-range"].getInt();
    }

    this->assignAoRange_ = 64;
    if (pdfParam["xc_ms_ao_range"].getStr().empty() != true) {
        this->assignAoRange_ = pdfParam["xc_ms_ao_range"].getInt();
    }

    this->assignJobsPerProc_ = 2;
    if (pdfParam["xc_ms_job_per_proc"].getStr().empty() != true) {
        this->assignJobsPerProc_ = pdfParam["xc_ms_job_per_proc"].getInt();
    }
    
    this->densityMatrixCacheMemSize_ = 100 * 1024UL * 1024UL; // 100MB

    this->calcMode_ = pdfParam["grid_calcmode"].getInt();
}

DfCalcGridX_Parallel::~DfCalcGridX_Parallel()
{
}


// void DfCalcGridX_Parallel::saveGridMatrix(const TlMatrix& gridMat)
// {
//     if (this->m_bIsUpdateXC == true) {
//         TlCommunicate& rComm = TlCommunicate::getInstance();
        
//         const index_type numOfGrids = gridMat.getNumOfRows();
//         int startCol = 0;
//         int endCol = 0;
//         if (this->functionalType_ == DfXCFunctional::LDA) {
//             startCol = GM_LDA_RHO_ALPHA;
//             endCol = GM_LDA_RHO_ALPHA + 1;
//             if (this->m_nMethodType == METHOD_UKS) {
//                 endCol = GM_LDA_RHO_BETA + 1;
//             }
//         } else {
//             assert(this->functionalType_ == DfXCFunctional::GGA);
//             startCol = GM_GGA_RHO_ALPHA;
//             endCol = GM_GGA_GRAD_RHO_Z_ALPHA + 1;
//             if (this->m_nMethodType == METHOD_UKS) {
//                 endCol = GM_GGA_GRAD_RHO_Z_BETA + 1;
//             }
//         }
//         TlMatrix commMat = gridMat.getBlockMatrix(0, numOfGrids,
//                                                   startCol, endCol - startCol);
        
//         rComm.allReduce_SUM(commMat);
        
//         if (rComm.isMaster() == true) {
//             TlMatrix newMat = gridMat;
//             newMat.setBlockMatrix(0, startCol, commMat);
//             DfCalcGridX::saveGridMatrix(this->m_nIteration, newMat);
//         }
//     }
// }


// void DfCalcGridX_Parallel::backupGridData()
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     rComm.barrier();
// }


// void DfCalcGridX_Parallel::flushGridData()
// {
//     this->logger(" flush grid data.\n");
//     typedef std::map<GridDataManager::ChunkType, std::map<int, std::vector<double> > > physicalValuesType;
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int nProc = rComm.getNumOfProc();

//     enum {
//         FINISH_MSG = -1
//     };

//     if (rComm.isMaster() == true) {
//         GridDataManager gdm(this->getGridDataFilePath());

//         // for own data
//         physicalValuesType::const_iterator pEnd = this->physicalValues_.end();
//         for (physicalValuesType::const_iterator p = this->physicalValues_.begin(); p != pEnd; ++p) {
//             const GridDataManager::ChunkType chunkType = p->first;
//             std::map<int, std::vector<double> >::const_iterator qEnd = p->second.end();
//             for (std::map<int, std::vector<double> >::const_iterator q = p->second.begin(); q != qEnd; ++q) {
//                 gdm.setData(q->first, chunkType, q->second);
//             }
//         }
//         this->physicalValues_.clear();

//         // message from slave
//         int finishProcCount = 0;
//         while (finishProcCount < (nProc -1)) {
//             int msg = 0;
//             int src = 0;
//             rComm.receiveDataFromAnySource(msg, &src);

//             int ack_msg = 0; // ACKnowledgement mssage
//             rComm.sendData(ack_msg, src);

//             if (msg >= 0) {
//                 const GridDataManager::ChunkType chunkType = (GridDataManager::ChunkType)msg;
//                 int atom = 0;
//                 std::vector<double> data;
//                 rComm.receiveData(atom, src);
//                 rComm.receiveData(data, src);
//                 gdm.setData(atom, chunkType, data);
//             } else if (msg == FINISH_MSG) {
//                 ++finishProcCount;
//                 continue;
//             }
//         }
//     } else {
//         const int root = 0;
//         int ack_msg = 0;

//         physicalValuesType::const_iterator pEnd = this->physicalValues_.end();
//         for (physicalValuesType::const_iterator p = this->physicalValues_.begin(); p != pEnd; ++p) {
//             const int chunkType = (int)p->first;
//             std::map<int, std::vector<double> >::const_iterator qEnd = p->second.end();
//             for (std::map<int, std::vector<double> >::const_iterator q = p->second.begin(); q != qEnd; ++q) {
//                 const int atom = q->first;
//                 const std::vector<double>& data = q->second;
//                 rComm.sendData(chunkType, root);
//                 rComm.receiveData(ack_msg, root);
//                 assert(ack_msg == 0);

//                 rComm.sendData(atom, root);
//                 rComm.sendData(data, root);
//             }
//         }

//         const int finishMsg = FINISH_MSG;
//         rComm.sendData(finishMsg, root);
//         rComm.receiveData(ack_msg, root);
//         assert(ack_msg == 0);

//         this->physicalValues_.clear();
//     }
// }


void DfCalcGridX_Parallel::defineCutOffValues(const TlSymmetricMatrix& P)
{
    // TODO: 並列化
    DfCalcGridX::defineCutOffValues(P);
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlSymmetricMatrix& PA,
                                              const TlSymmetricMatrix& PB)
{
    // TODO: 並列化
    DfCalcGridX::defineCutOffValues(PA, PB);
}



// =============================================================================
// ScaLAPACK版
//
// memo:
// - grid dataはmasterノードのみ保存する。
// =============================================================================
// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P,
//                                                          DfFunctional_LDA* pFunctional,
//                                                          TlDistributeSymmetricMatrix* pF)
// {
    // if (this->enableExperimentalCode_ != true) {
    //     return this->calcXCIntegForFockAndEnergy1(P, pFunctional, pF);
    // } else {
    //     this->logger(" !!! experimental code !!!\n");
    //     return this->calcXCIntegForFockAndEnergy2(P, pFunctional, pF);
    // }
// }


// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy1(const TlDistributeSymmetricMatrix& P,
//                                                           DfFunctional_LDA* pFunctional,
//                                                           TlDistributeSymmetricMatrix* pF)
// {
//     assert(pFunctional != NULL);
//     assert(pF != NULL);

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     this->physicalValues_.clear();
//     this->defineCutOffValues(P);

//     // gird情報(座標, 重み)をMasterから全ノードに配布
//     std::vector<int> gridCounts;
//     std::vector<double> gridInfo;
//     this->broadcastGridInfo(&gridCounts, &gridInfo);
    
//     // 各ノードはローカルScaLAPACK行列に対し, 各グリッド点の電子密度, 勾配を計算
//     const std::vector<int> AO_list = this->getScaLapackLocalAOs(P);
    
//     std::vector<double> rhoAs;
//     this->calcPhys(P, AO_list, gridInfo,
//                    &rhoAs);
//     for (std::size_t i = 0; i < rhoAs.size(); ++i) {
//         rhoAs[i] *= 0.5;
//     }
    
//     // 各グリッド点の電子密度, 勾配をMasterノードに集計/再配布
//     if ((this->m_bIsUpdateXC == true) && (this->m_nIteration > 1)) {
//         if (rComm.isMaster() == true) {
//             // 前回の電子密度を加算
//             this->addPreviousPhys(&rhoAs);
//         }
//     }
//     rComm.allReduce_SUM(rhoAs);
//     if (this->m_bIsUpdateXC == true) {
//         if (rComm.isMaster() == true) {
//             // 今回の電子密度を保存
//             this->saveCurrentPhys(gridCounts, rhoAs);
//         }
//     }
    
//     // 各グリッド点の電子密度, 勾配からXC行列を作成する
//     TlSparseSymmetricMatrix tmpF(this->m_nNumOfAOs);
//     double energy = this->buildK(gridCounts, gridInfo,
//                                  rhoAs,
//                                  pFunctional, &tmpF);
//     pF->resize(this->m_nNumOfAOs);
//     pF->mergeSparseMatrix(tmpF);
//     rComm.allReduce_SUM(energy);
    
//     return energy;
// }


// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy2(const TlDistributeSymmetricMatrix& P,
//                                                           DfFunctional_LDA* pFunctional,
//                                                           TlDistributeSymmetricMatrix* pF)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(pFunctional != NULL);
//     assert(pF != NULL);

//     this->logger(" numerical integral using distribute matrix.\n");

//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     this->physicalValues_.clear();
//     this->defineCutOffValues(P);

//     // gird情報(座標, 重み)をMasterから全ノードに配布
//     this->loggerTime(" broadcast grid information.");
//     std::vector<int> gridCounts;
//     std::vector<double> gridInfo;
//     this->broadcastGridInfo(&gridCounts, &gridInfo);
    
//     // 各ノードはローカルScaLAPACK行列に対し, 各グリッド点の電子密度, 勾配を計算
//     this->loggerTime(" calculate population of each grid.");
//     std::vector<double> rhoAs;
//     this->calcPhys_new(TlDistributeMatrix(P), gridInfo,
//                        &rhoAs);
//     for (std::size_t i = 0; i < rhoAs.size(); ++i) {
//         rhoAs[i] *= 0.5;
//     }
    
//     // 各グリッド点の電子密度, 勾配をMasterノードに集計/再配布
//     this->loggerTime(" load previous parameters.");
//     if ((this->m_bIsUpdateXC == true) && (this->m_nIteration > 1)) {
//         if (rComm.isMaster() == true) {
//             // 前回の電子密度を加算
//             this->addPreviousPhys(&rhoAs);
//         }
//     }
//     this->loggerTime(" reduce and broadcast the populations.");
//     rComm.allReduce_SUM(rhoAs);
//     this->loggerTime(" save current parameters.");
//     if (this->m_bIsUpdateXC == true) {
//         if (rComm.isMaster() == true) {
//             // 今回の電子密度を保存
//             this->saveCurrentPhys(gridCounts, rhoAs);
//         }
//     }
    
//     // 各グリッド点の電子密度, 勾配からXC行列を作成する
//     this->loggerTime(" build the KS matrix of XC term");
//     TlSparseSymmetricMatrix tmpF(this->m_nNumOfAOs);
//     double energy = this->buildK(gridCounts, gridInfo,
//                                  rhoAs,
//                                  pFunctional, &tmpF);
//     pF->resize(this->m_nNumOfAOs);
//     pF->mergeSparseMatrix(tmpF);
//     rComm.allReduce_SUM(energy);
//     this->loggerTime("\n");
    
//     return energy;
// }


// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
//                                                          const TlDistributeSymmetricMatrix& P_B,
//                                                          DfFunctional_LDA* pFunctional,
//                                                          TlDistributeSymmetricMatrix* pF_A,
//                                                          TlDistributeSymmetricMatrix* pF_B)
// {
    // assert(pFunctional != NULL);
    // assert(pF_A != NULL);
    // assert(pF_B != NULL);

    // TlCommunicate& rComm = TlCommunicate::getInstance();
    // if (rComm.isMaster() == true) {
    //     DfCalcGridX::backupGridData();
    // }
    // this->physicalValues_.clear();
    // this->defineCutOffValues(P_A, P_B);

    // // gird情報(座標, 重み)をMasterから全ノードに配布
    // std::vector<int> gridCounts;
    // std::vector<double> gridInfo;
    // this->broadcastGridInfo(&gridCounts, &gridInfo);
    
    // // 各ノードはローカルScaLAPACK行列に対し, 各グリッド点の電子密度, 勾配を計算
    // const std::vector<int> AO_list = this->getScaLapackLocalAOs(P_A); // PBも同じ
    // std::vector<double> rhoAs;
    // std::vector<double> rhoBs;
    // this->calcPhys(P_A, AO_list, gridInfo,
    //                &rhoAs);
    // this->calcPhys(P_B, AO_list, gridInfo,
    //                &rhoBs);
    
    // // 各グリッド点の電子密度, 勾配をMasterノードに集計/再配布
    // if ((this->m_bIsUpdateXC == true) && (this->m_nIteration > 1)) {
    //     if (rComm.isMaster() == true) {
    //         // 前回の電子密度を加算
    //         this->addPreviousPhys(&rhoAs, &rhoBs);
    //     }
    // }
    // rComm.allReduce_SUM(rhoAs);
    // if (this->m_bIsUpdateXC == true) {
    //     if (rComm.isMaster() == true) {
    //         // 今回の電子密度を保存
    //         this->saveCurrentPhys(gridCounts, rhoAs, rhoBs);
    //     }
    // }
    
    // // 各グリッド点の電子密度, 勾配からXC行列を作成する
    // TlSparseSymmetricMatrix tmpFA(this->m_nNumOfAOs);
    // TlSparseSymmetricMatrix tmpFB(this->m_nNumOfAOs);
    // double energy = this->buildK(gridCounts, gridInfo,
    //                              rhoAs, rhoBs,
    //                              pFunctional, &tmpFA, &tmpFB);
    // pF_A->resize(this->m_nNumOfAOs);
    // pF_B->resize(this->m_nNumOfAOs);
    // pF_A->mergeSparseMatrix(tmpFA);
    // pF_B->mergeSparseMatrix(tmpFB);
    // rComm.allReduce_SUM(energy);
    
//     return 0.0;
// }


// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P,
//                                                          DfFunctional_GGA* pFunctional,
//                                                          TlDistributeSymmetricMatrix* pF)
// {
    // this->logger(" numerical integral using distribute matrix.\n");

    // double answer = 0.0;
    // if (this->calcMode_ == 2) {
    //     this->logger(" local density matrix is used.\n");
    //     answer = this->calcXC_DC(P, pFunctional, pF);
    // } else {
    //     this->logger(" global density matrix is transfered in the background.\n");
    //     answer = this->calcXC_BG(P, pFunctional, pF);
    // }

    // return answer;
//     return 0.0;
// }


// double DfCalcGridX_Parallel::calcXC_BG(const TlDistributeSymmetricMatrix& P,
//                                        DfFunctional_GGA* pFunctional,
//                                        TlDistributeSymmetricMatrix* pF)
// {
//     assert(pFunctional != NULL);
//     assert(pF != NULL);

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     this->physicalValues_.clear();
//     this->defineCutOffValues(P);

//     // gird情報(座標, 重み)をMasterから全ノードに配布
//     this->loggerTime(" broadcast grid information.");
//     std::vector<int> gridCounts;
//     std::vector<double> gridInfo;
//     this->broadcastGridInfo(&gridCounts, &gridInfo);
    
//     // 各ノードはローカルScaLAPACK行列に対し, 各グリッド点の電子密度, 勾配を計算
//     this->loggerTime(" calculate population of each grid.");
//     //const std::vector<int> AO_list = this->getScaLapackLocalAOs(P);
//     std::vector<double> rhoAs;
//     std::vector<double> gradRhoAXs;
//     std::vector<double> gradRhoAYs;
//     std::vector<double> gradRhoAZs;
//     // this->calcPhys(P, AO_list, gridInfo,
//     //                &rhoAs, &gradRhoAXs, &gradRhoAYs, &gradRhoAZs);
//     this->calcPhys(P, gridInfo,
//                    &rhoAs, &gradRhoAXs, &gradRhoAYs, &gradRhoAZs);

//     const std::size_t numOfAllGrids = rhoAs.size();
// #pragma omp parallel for
//     for (std::size_t i = 0; i < numOfAllGrids; ++i) {
//         rhoAs[i] *= 0.5;
//         gradRhoAXs[i] *= 0.5;
//         gradRhoAYs[i] *= 0.5;
//         gradRhoAZs[i] *= 0.5;
//     }
    
//     // 各グリッド点の電子密度, 勾配をMasterノードに集計/再配布
//     this->loggerTime(" load previous parameters.");
//     if ((this->m_bIsUpdateXC == true) && (this->m_nIteration > 1)) {
//         if (rComm.isMaster() == true) {
//             // 前回の電子密度を加算
//             this->addPreviousPhys(&rhoAs, &gradRhoAXs, &gradRhoAYs, &gradRhoAZs);
//         }
//     }
//     this->loggerTime(" reduce and broadcast the populations.");
//     {
//         std::vector<double> tmp(numOfAllGrids * 4);
//         std::copy(rhoAs.begin(), rhoAs.end(), tmp.begin());
//         std::copy(gradRhoAXs.begin(), gradRhoAXs.end(), tmp.begin() + numOfAllGrids);
//         std::copy(gradRhoAYs.begin(), gradRhoAYs.end(), tmp.begin() + numOfAllGrids *2);
//         std::copy(gradRhoAZs.begin(), gradRhoAZs.end(), tmp.begin() + numOfAllGrids *3);
//         rComm.allReduce_SUM((double*)&(tmp[0]), numOfAllGrids * 4);
//         std::copy(tmp.begin(),                     tmp.begin() + numOfAllGrids,     rhoAs.begin());
//         std::copy(tmp.begin() + numOfAllGrids,     tmp.begin() + numOfAllGrids * 2, gradRhoAXs.begin());
//         std::copy(tmp.begin() + numOfAllGrids * 2, tmp.begin() + numOfAllGrids * 3, gradRhoAYs.begin());
//         std::copy(tmp.begin() + numOfAllGrids * 3, tmp.begin() + numOfAllGrids * 4, gradRhoAZs.begin());
//     }
//     this->loggerTime(" save current parameters.");
//     if (this->m_bIsUpdateXC == true) {
//         if (rComm.isMaster() == true) {
//             // 今回の電子密度を保存
//             this->saveCurrentPhys(gridCounts, rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs);
//         }
//     }
    
//     // 各グリッド点の電子密度, 勾配からXC行列を作成する
// //     this->loggerTime(" build the KS matrix of XC term.\n");
// //     TlSparseSymmetricMatrix tmpF(this->m_nNumOfAOs);
// //     double energy = this->buildK(gridCounts, gridInfo,
// //                                  rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
// //                                  pFunctional, &tmpF);
// //     pF->resize(this->m_nNumOfAOs);

// //     this->loggerTime(" reduce the KS matrix of XC term.\n");
// //     pF->mergeSparseMatrix(tmpF);
// //     rComm.allReduce_SUM(energy);
//     double energy = this->buildK_rev2(gridCounts, gridInfo,
//                                       rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
//                                       pFunctional, pF);

//     this->loggerTime(" finish.\n");
//     return energy;
// }


// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& PA,
//                                                          const TlDistributeSymmetricMatrix& PB,
//                                                          DfFunctional_GGA* pFunctional,
//                                                          TlDistributeSymmetricMatrix* pFA,
//                                                          TlDistributeSymmetricMatrix* pFB)
// {
    // assert(pFunctional != NULL);
    // assert(pFA != NULL);
    // assert(pFB != NULL);

    // TlCommunicate& rComm = TlCommunicate::getInstance();
    // if (rComm.isMaster() == true) {
    //     DfCalcGridX::backupGridData();
    // }
    // this->physicalValues_.clear();
    // this->defineCutOffValues(PA, PB);

    // // gird情報(座標, 重み)をMasterから全ノードに配布
    // std::vector<int> gridCounts;
    // std::vector<double> gridInfo;
    // this->broadcastGridInfo(&gridCounts, &gridInfo);
    
    // // 各ノードはローカルScaLAPACK行列に対し, 各グリッド点の電子密度, 勾配を計算
    // std::vector<double> rhoAs;
    // std::vector<double> gradRhoAXs;
    // std::vector<double> gradRhoAYs;
    // std::vector<double> gradRhoAZs;
    // std::vector<double> rhoBs;
    // std::vector<double> gradRhoBXs;
    // std::vector<double> gradRhoBYs;
    // std::vector<double> gradRhoBZs;
    // this->calcPhys(PA, gridInfo,
    //                &rhoAs, &gradRhoAXs, &gradRhoAYs, &gradRhoAZs);
    // this->calcPhys(PB, gridInfo,
    //                &rhoBs, &gradRhoBXs, &gradRhoBYs, &gradRhoBZs);
    
    // // 各グリッド点の電子密度, 勾配をMasterノードに集計/再配布
    // if ((this->m_bIsUpdateXC == true) && (this->m_nIteration > 1)) {
    //     if (rComm.isMaster() == true) {
    //         // 前回の電子密度を加算
    //         this->addPreviousPhys(&rhoAs, &gradRhoAXs, &gradRhoAYs, &gradRhoAZs,
    //                               &rhoBs, &gradRhoBXs, &gradRhoBYs, &gradRhoBZs);
    //     }
    // }
    // rComm.allReduce_SUM(rhoAs);
    // rComm.allReduce_SUM(gradRhoAXs);
    // rComm.allReduce_SUM(gradRhoAYs);
    // rComm.allReduce_SUM(gradRhoAZs);
    // rComm.allReduce_SUM(rhoBs);
    // rComm.allReduce_SUM(gradRhoBXs);
    // rComm.allReduce_SUM(gradRhoBYs);
    // rComm.allReduce_SUM(gradRhoBZs);
    // if (this->m_bIsUpdateXC == true) {
    //     if (rComm.isMaster() == true) {
    //         // 今回の電子密度を保存
    //         this->saveCurrentPhys(gridCounts,
    //                               rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
    //                               rhoBs, gradRhoBXs, gradRhoBYs, gradRhoBZs);
    //     }
    // }
    
    // // 各グリッド点の電子密度, 勾配からXC行列を作成する
    // TlSparseSymmetricMatrix tmpFA(this->m_nNumOfAOs);
    // TlSparseSymmetricMatrix tmpFB(this->m_nNumOfAOs);
    // double energy = this->buildK(gridCounts, gridInfo,
    //                              rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
    //                              rhoBs, gradRhoBXs, gradRhoBYs, gradRhoBZs,
    //                              pFunctional, &tmpFA, &tmpFB);
    // pFA->resize(this->m_nNumOfAOs);
    // pFB->resize(this->m_nNumOfAOs);
    // pFA->mergeSparseMatrix(tmpFA);
    // pFB->mergeSparseMatrix(tmpFB);
    // rComm.allReduce_SUM(energy);
    
    // return energy;
//     return 0.0;
// }


void DfCalcGridX_Parallel::defineCutOffValues(const TlDistributeSymmetricMatrix& P)
{
    const double maxValueOfP = std::max(P.getMaxAbsoluteElement(), 1.0E-16);
    if (maxValueOfP < 1.0) {
        this->m_densityCutOffValueA /= maxValueOfP;
    }
    this->logger(TlUtils::format(" density cutoff value = %e\n",
                                 this->m_densityCutOffValueA));
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.broadcast(this->m_densityCutOffValueA);
}


void DfCalcGridX_Parallel::defineCutOffValues(const TlDistributeSymmetricMatrix& PA,
                                              const TlDistributeSymmetricMatrix& PB)
{
    const double maxValueOfPA = std::max(PA.getMaxAbsoluteElement(), 1.0E-16);
    const double maxValueOfPB = std::max(PB.getMaxAbsoluteElement(), 1.0E-16);
    if (maxValueOfPA < 1.0) {
        this->m_densityCutOffValueA /= maxValueOfPA;
    }
    if (maxValueOfPB < 1.0) {
        this->m_densityCutOffValueB /= maxValueOfPB;
    }
    this->logger(TlUtils::format(" density cutoff value(alpha) = %e\n",
                                 this->m_densityCutOffValueA));
    this->logger(TlUtils::format(" density cutoff value(beta ) = %e\n",
                                 this->m_densityCutOffValueB));

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.broadcast(this->m_densityCutOffValueA);
    rComm.broadcast(this->m_densityCutOffValueB);
}


// void DfCalcGridX_Parallel::broadcastGridInfo(std::vector<int>* pGridCounts,
//                                              std::vector<double>* pGridInfo)
// {
//     assert(pGridCounts != NULL);
//     assert(pGridInfo != NULL);
    
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int numOfRealAtoms = this->numOfRealAtoms_;
    
//     if (rComm.isMaster() == true) {
//         GridDataManager gdm(this->getGridDataFilePath());
//         pGridCounts->clear();
//         pGridCounts->resize(numOfRealAtoms);

//         std::size_t allGrid = 0;
//         for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//             const std::vector<double> coordX = gdm.getData(atom, GridDataManager::COORD_X);
//             const std::vector<double> coordY = gdm.getData(atom, GridDataManager::COORD_Y);
//             const std::vector<double> coordZ = gdm.getData(atom, GridDataManager::COORD_Z);
//             const std::vector<double> weight = gdm.getData(atom, GridDataManager::GRID_WEIGHT);
//             assert(coordX.size() == coordY.size());
//             assert(coordY.size() == coordZ.size());
//             assert(coordZ.size() == weight.size());
            
//             const int numOfGrids = coordX.size();
//             (*pGridCounts)[atom] = numOfGrids;

//             pGridInfo->resize(pGridInfo->size() + numOfGrids*4);
//             for (int grid = 0; grid < numOfGrids; ++grid) {
//                 (*pGridInfo)[allGrid++] = coordX[grid];
//                 (*pGridInfo)[allGrid++] = coordY[grid];
//                 (*pGridInfo)[allGrid++] = coordZ[grid];
//                 (*pGridInfo)[allGrid++] = weight[grid];
//             }
//         }
//         assert(allGrid == pGridInfo->size());
//     }

//     rComm.broadcast(*pGridCounts);
//     rComm.broadcast(*pGridInfo);
// }


// void DfCalcGridX_Parallel::calcPhys(const TlDistributeSymmetricMatrix& P,
//                                     const std::vector<int>& AO_list,
//                                     const std::vector<double>& gridInfo,
//                                     std::vector<double>* pRhoAs)
// {
//     const int numOfAllGrids = gridInfo.size() / 4;
//     pRhoAs->clear();
//     pRhoAs->resize(numOfAllGrids, 0.0);

//     for (int grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const TlPosition position(x, y, z);

//         // calc phi table
//         std::vector<WFGrid> phis;
//         this->getPhiTable(position, AO_list, &phis);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());

//         // get rho at grid point
//         double rhoA = 0.0;
//         this->getRhoAtGridPoint(P, phis, &rhoA);
//         (*pRhoAs)[grid] += rhoA;
//     }
// }


// void DfCalcGridX_Parallel::calcPhys_new(const TlDistributeMatrix& P,
//                                         const std::vector<double>& gridInfo,
//                                         std::vector<double>* pRhoVals)
// {
//     std::vector<index_type> rowIndexes = P.getRowIndexTable();
//     std::vector<index_type> colIndexes = P.getColIndexTable();
//     const std::size_t numOfRows = rowIndexes.size();
//     const std::size_t numOfCols = colIndexes.size();

//     const TlMatrix localP = P.getLocalMatrix();
    
//     const std::size_t numOfGrids = gridInfo.size() / 4;
//     pRhoVals->resize(numOfGrids);
//     for (int grid = 0; grid < numOfGrids; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const TlPosition r(x, y, z);

//         TlMatrix wfCoef_row = this->getWaveFunctionCoef(rowIndexes, r);
//         TlMatrix wfCoef_col = this->getWaveFunctionCoef(colIndexes, r);
//         wfCoef_col.transpose();
//         TlMatrix coefMatrix = wfCoef_row * wfCoef_col;
//         assert(coefMatrix.getNumOfRows() == numOfRows);
//         assert(coefMatrix.getNumOfCols() == numOfCols);
//         coefMatrix.dot(localP);

//         (*pRhoVals)[grid] = coefMatrix.sum();
//     }    
// }


// void DfCalcGridX_Parallel::calcPhys_new(const TlDistributeMatrix& P,
//                                         const std::vector<double>& gridInfo,
//                                         std::vector<double>* pRhoVals,
//                                         std::vector<double>* pGradRhoXVals,
//                                         std::vector<double>* pGradRhoYVals,
//                                         std::vector<double>* pGradRhoZVals)
// {
//     std::vector<index_type> rowIndexes = P.getRowIndexTable();
//     std::vector<index_type> colIndexes = P.getColIndexTable();
//     const std::size_t numOfRows = rowIndexes.size();
//     const std::size_t numOfCols = colIndexes.size();

//     const TlMatrix localP = P.getLocalMatrix();
//     assert(localP.getNumOfRows() == numOfRows);
//     assert(localP.getNumOfCols() == numOfCols);
    
//     const std::size_t numOfGrids = gridInfo.size() / 4;
//     pRhoVals->resize(numOfGrids);
//     pGradRhoXVals->resize(numOfGrids);
//     pGradRhoYVals->resize(numOfGrids);
//     pGradRhoZVals->resize(numOfGrids);
//     for (int grid = 0; grid < numOfGrids; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const TlPosition r(x, y, z);

//         TlMatrix wfCoef_row = this->getWaveFunctionCoef(rowIndexes, r);
//         TlMatrix wfCoef_col;
//         TlMatrix wfCoef_gradX;
//         TlMatrix wfCoef_gradY;
//         TlMatrix wfCoef_gradZ;
//         this->getWaveFunctionCoef(colIndexes, r,
//                                   &wfCoef_col,
//                                   &wfCoef_gradX, &wfCoef_gradY, &wfCoef_gradZ);
//         wfCoef_col = this->getWaveFunctionCoef(colIndexes, r);
//         // rho
//         {
//             wfCoef_col.transpose();
//             TlMatrix coefMatrix = wfCoef_row * wfCoef_col;
//             coefMatrix.dot(localP);
//             (*pRhoVals)[grid] = coefMatrix.sum();
//         }
//         // gradRho
//         {
//             wfCoef_gradX.transpose();
//             TlMatrix coefMatrix = wfCoef_row * wfCoef_gradX;
//             coefMatrix.dot(localP);
//             (*pGradRhoXVals)[grid] = 2.0 * coefMatrix.sum();
//         }
//         {
//             wfCoef_gradY.transpose();
//             TlMatrix coefMatrix = wfCoef_row * wfCoef_gradY;
//             coefMatrix.dot(localP);
//             (*pGradRhoYVals)[grid] = 2.0 * coefMatrix.sum();
//         }
//         {
//             wfCoef_gradZ.transpose();
//             TlMatrix coefMatrix = wfCoef_row * wfCoef_gradZ;
//             coefMatrix.dot(localP);
//             (*pGradRhoZVals)[grid] = 2.0 * coefMatrix.sum();
//         }
//     }    
// }


// TlMatrix DfCalcGridX_Parallel::getWaveFunctionCoef(const std::vector<index_type>& AO_indexes,
//                                                    const TlPosition& gridPosition)
// {
//     const std::size_t size = AO_indexes.size();
//     TlMatrix wfCoef(size, 1);
//     for (int aoIndex = 0; aoIndex < size; ++aoIndex) {
//         const index_type orb = AO_indexes[aoIndex];

//         const TlPosition r = gridPosition - this->m_tlOrbInfo.getPosition(orb);
//         const double r2 = r.squareDistanceFrom();
//         const int basisType = this->m_tlOrbInfo.getBasisType(orb);
//         const double prefactor = this->getPrefactor(basisType, r);

//         double value = 0.0;
//         const int numOfPGTOs = this->m_tlOrbInfo.getCgtoContraction(orb);
//         for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
//             const double alpha = this->m_tlOrbInfo.getExponent(orb, pgtoIndex);
//             const double exponent = alpha * r2;
//             const double coef = this->m_tlOrbInfo.getCoefficient(orb, pgtoIndex);
//             value += coef * std::exp(- exponent);
//         }

//         wfCoef.set(aoIndex, 0, prefactor * value);
//     }
//     return wfCoef;
// }


// void DfCalcGridX_Parallel::getWaveFunctionCoef(const std::vector<index_type>& AO_indexes,
//                                                const TlPosition& gridPosition,
//                                                TlMatrix* pWF,
//                                                TlMatrix* pGradWF_X,
//                                                TlMatrix* pGradWF_Y,
//                                                TlMatrix* pGradWF_Z)
// {
//     const std::size_t size = AO_indexes.size();
//     pWF->resize(size, 1);
//     pGradWF_X->resize(size, 1);
//     pGradWF_Y->resize(size, 1);
//     pGradWF_Z->resize(size, 1);

//     double prefactorX = 0.0;
//     double prefactorY = 0.0;
//     double prefactorZ = 0.0;
//     for (int aoIndex = 0; aoIndex < size; ++aoIndex) {
//         const index_type orb = AO_indexes[aoIndex];

//         const TlPosition r = gridPosition - this->m_tlOrbInfo.getPosition(orb);
//         const double r2 = r.squareDistanceFrom();
//         const int basisType = this->m_tlOrbInfo.getBasisType(orb);
//         const double prefactor = this->getPrefactor(basisType, r);

//         double value = 0.0;
//         double gradX = 0.0;
//         double gradY = 0.0;
//         double gradZ = 0.0;
//         const int numOfPGTOs = this->m_tlOrbInfo.getCgtoContraction(orb);
//         for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
//             const double alpha = this->m_tlOrbInfo.getExponent(orb, pgtoIndex);
//             const double exponent = alpha * r2;
//             const double e = std::exp(- exponent);
//             const double coef = this->m_tlOrbInfo.getCoefficient(orb, pgtoIndex);
//             this->getPrefactorForDerivative(basisType, alpha, r,
//                                             &prefactorX, &prefactorY, &prefactorZ);

//             value += coef * e;
//             gradX += coef * prefactorX * e;
//             gradY += coef * prefactorY * e;
//             gradZ += coef * prefactorZ * e;
//         }
//         pWF->set(aoIndex, 0, prefactor* value);
//         pGradWF_X->set(aoIndex, 0, gradX);
//         pGradWF_Y->set(aoIndex, 0, gradY);
//         pGradWF_Z->set(aoIndex, 0, gradZ);
//     }
// }


// void DfCalcGridX_Parallel::calcPhys(const TlDistributeSymmetricMatrix& P,
//                                     const std::vector<double>& gridInfo,
//                                     std::vector<double>* pRhoAs,
//                                     std::vector<double>* pGradRhoAXs,
//                                     std::vector<double>* pGradRhoAYs,
//                                     std::vector<double>* pGradRhoAZs)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(rComm.checkNonBlockingCommunications() == true);
//     // const int numOfProcs = rComm.getNumOfProc();

//     const int numOfAllGrids = gridInfo.size() / 4;
//     pRhoAs->clear();
//     pRhoAs->resize(numOfAllGrids, 0.0);
//     pGradRhoAXs->clear();
//     pGradRhoAXs->resize(numOfAllGrids, 0.0);
//     pGradRhoAYs->clear();
//     pGradRhoAYs->resize(numOfAllGrids, 0.0);
//     pGradRhoAZs->clear();
//     pGradRhoAZs->resize(numOfAllGrids, 0.0);

//     // timer
//     TlTime totalTimer;
//     TlTime calcTimer;
//     calcTimer.stop();
    
//     this->loggerTime(" calculate data on grids");
//     this->logger(TlUtils::format("  assign AO range = %d\n", this->assignAoRange_));
//     this->logger(TlUtils::format("  assign jobs per proc = %d\n", this->assignJobsPerProc_));
//     //const static int JOB_PROTOCOL_SIZE = 4; // index=0: shell index
//     if (rComm.isMaster() == true) {
//         this->calcPhys_Master(P);
//     } else {
//         this->calcPhys_Slave(P, gridInfo,
//                              pRhoAs, pGradRhoAXs, pGradRhoAYs, pGradRhoAZs);
//     }

//     // 後始末(master, slave 共通)
//     P.getPartialMatrixX((TlPartialSymmetricMatrix*)NULL, true);

//     // timer
// //     for (int proc = 0; proc < numOfProcs; ++proc) {
// //         if (proc == rank) {
// //             const std::string timeStr = TlUtils::format("  [%5d] Total: %.2f sec. Calc: %.2f sec.\n",
// //                                                         proc,
// //                                                         totalTimer.getElapseTime(),
// //                                                         calcTimer.getElapseTime());
// //             this->loggerP(timeStr);
// //         }
// //         rComm.barrier();
// //     }
    
    
//     assert(rComm.checkNonBlockingCommunications() == true);
// }


// void DfCalcGridX_Parallel::calcPhys_Master(const TlDistributeSymmetricMatrix& P)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int numOfProcs = rComm.getNumOfProc();
//     const int jobsPerProc = this->assignJobsPerProc_; // 1プロセスあたりの同時ジョブ数
    
//     this->getQueueX(NULL, NULL, NULL, NULL,
//                     NULL, true); // initialize only

//     // slaveからのjobのリクエスト
//     int cmdRequestJob = 0;
//     bool isWaitingCmdRequestJob = false;

//     // ジョブキュー
//     typedef std::list<int> AssignJobQueueType;
//     AssignJobQueueType assignJobQueue;

//     //
//     std::vector<bool> isWaitingJobList(numOfProcs);
//     std::vector<index_type*> jobList(numOfProcs);
//     for (int i = 0; i < numOfProcs; ++i) {
//         jobList[i] = new int[JOB_PROTOCOL_SIZE];
//     }

//     int terminateJobProc = 0;
//     bool isTerminate = false;
//     bool isBreakLoop = false;
//     int prevProgress = -1;
    
//     while (isBreakLoop != true) {
//         int src = 0;
//         // phase1
//         if (isWaitingCmdRequestJob != true) {
//             rComm.iReceiveDataFromAnySource(cmdRequestJob, TAG_REQUEST_JOB);
//             isWaitingCmdRequestJob = true;
//             this->log_.debug(" [0] RECV[**] request job\n");
//         }
        
//         // phase2
//         if ((isWaitingCmdRequestJob == true) &&
//             (rComm.test(cmdRequestJob, &src) == true)) {
//             rComm.wait(cmdRequestJob);
//             isWaitingCmdRequestJob = false;
//             this->log_.debug(TlUtils::format(" [0] RECV[OK] request job from [%d].\n",
//                                              src));
            
//             assignJobQueue.push_back(src);
//         }
        
//         // phase3
//         AssignJobQueueType::iterator itEnd = assignJobQueue.end();
//         for (AssignJobQueueType::iterator it = assignJobQueue.begin(); it != itEnd; ) {
//             const int slave = *it;
            
//             index_type startShell_I = 0;
//             index_type startShell_J = 0;
//             index_type endShell_I = 0;
//             index_type endShell_J = 0;
//             double progress = 0.0;
//             const bool isContinued = this->getQueueX(&startShell_I, &startShell_J,
//                                                      &endShell_I, &endShell_J,
//                                                      &progress);
            
//             // for progress
//             {
//                 const int currentProgress = int(progress / 10.0);
//                 if (currentProgress > prevProgress) {
//                     this->loggerTime(TlUtils::format("  %3d%% done.", currentProgress * 10));
//                     prevProgress = currentProgress;
//                 }
//             }
            
//             if (isContinued == true) {
//                 // assign job
//                 if (isWaitingJobList[slave] == true) {
//                     rComm.wait(jobList[slave]);
//                     isWaitingJobList[slave] = false;
//                     this->log_.debug(TlUtils::format(" [0] SEND[OK] assign job to [%d].\n",
//                                                      slave));
//                 }
//                 jobList[slave][0] = startShell_I;
//                 jobList[slave][1] = startShell_J;
//                 jobList[slave][2] = endShell_I;
//                 jobList[slave][3] = endShell_J;
//                 rComm.iSendDataX(jobList[slave], JOB_PROTOCOL_SIZE, slave, TAG_ASSIGN_JOB);
//                 isWaitingJobList[slave] = true;
//                 this->log_.debug(TlUtils::format(" [0] SEND[**] assign job to [%d]: (%d, %d) -> (%d, %d).\n",
//                                                  slave,
//                                                  startShell_I, startShell_J, endShell_I, endShell_J));
//             } else {
//                 // finish job
//                 ++terminateJobProc;
//                 if (terminateJobProc >= (numOfProcs -1) * jobsPerProc) {
//                     this->log_.debug(" [0] SEND[OK] terminate to others.\n");
//                     int cmdTerminateSlave = -1;
//                     for (int i = 0; i < numOfProcs -1; ++i) {
//                         rComm.sendData(cmdTerminateSlave, i +1, TAG_TERMINATE_SLAVE);
//                     }
//                     isTerminate = true;
//                 }
//             }
            
//             it = assignJobQueue.erase(it);
//         }
        
//         P.getPartialMatrixX(NULL);
        
//         // terminate
//         if ((isTerminate == true) && (assignJobQueue.empty() == true)) {
//             isBreakLoop = true;
//         }
//     } // end while

//     // 後始末
//     for (int i = 0; i < numOfProcs; ++i) {
//         if (isWaitingJobList[i] == true) {
//             rComm.cancel(jobList[i]); 
//             isWaitingJobList[i] = false;
//         }
            
//         delete[] jobList[i];
//         jobList[i] = NULL;
//     }

//     if (isWaitingCmdRequestJob == true) {
//         rComm.cancel(cmdRequestJob);
//         isWaitingCmdRequestJob = false;
//     }
// }


// void DfCalcGridX_Parallel::calcPhys_Slave(const TlDistributeSymmetricMatrix& P,
//                                           const std::vector<double>& gridInfo,
//                                           std::vector<double>* pRhoAs,
//                                           std::vector<double>* pGradRhoAXs,
//                                           std::vector<double>* pGradRhoAYs,
//                                           std::vector<double>* pGradRhoAZs)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int root = 0;

//     // masterからのジョブの割り当て
//     const int jobsPerProc = this->assignJobsPerProc_;
//     std::vector<std::vector<int> > cmdAssignJobs(jobsPerProc);
//     std::vector<bool> isWaitingCmdAssignJobs(jobsPerProc);
//     for (int i = 0; i < jobsPerProc; ++i) {
//         cmdAssignJobs[i] = std::vector<int>(JOB_PROTOCOL_SIZE, 0);
//         isWaitingCmdAssignJobs[i] = false;
//     }
    
//     // masterからの終了メッセージ
//     int cmdTerminateSlave = 0;
//     bool isWaitingCmdTerminateSlave = false;
        
//     // ジョブのリスト
//     typedef std::list<DensityCalcJob> JobListType;
//     JobListType jobList;
        
//     // 初期メッセージ
//     std::vector<int> cmdRequestJob(jobsPerProc, 0);
//     std::vector<bool> isWaitingCmdRequestJob(jobsPerProc);
//     for (int i = 0; i < jobsPerProc; ++i) {
//         rComm.iSendData(cmdRequestJob[i], root, TAG_REQUEST_JOB);
//         isWaitingCmdRequestJob[i] = true;
//     }

//     bool isTerminate = false;
//     bool isBreakLoop = false;

// //     TlTime time_all;
// //     time_all.start();
// //     TlTime time_calc;
// //     time_calc.stop();
//     std::size_t numOfMisHits = 0;
//     std::size_t numOfChances = 0;
    
// #ifdef _OPENMP
//     // OpenMPが有効な場合、以下のforループはnestされる
//     const int ompNestedBackup = omp_get_nested();
//     omp_set_nested(1);
// #endif // _OPENMP

//     while (isBreakLoop != true) {
//         P.getPartialMatrixX(NULL);

//         // phase1
//         for (int i = 0; i < jobsPerProc; ++i) {
//             if (isWaitingCmdAssignJobs[i] != true) {
//                 rComm.iReceiveDataX(&(cmdAssignJobs[i][0]), JOB_PROTOCOL_SIZE, root, TAG_ASSIGN_JOB);
//                 isWaitingCmdAssignJobs[i] = true;
//                 this->log_.debug(TlUtils::format(" [%d] RECV[**] assign job[%d] from [0]\n",
//                                                  rComm.getRank(), i));
//             }
//         }
//         if (isWaitingCmdTerminateSlave != true) {
//             rComm.iReceiveData(cmdTerminateSlave, root, TAG_TERMINATE_SLAVE);
//             isWaitingCmdTerminateSlave = true;
//             this->log_.debug(TlUtils::format(" [%d] RECV[**] terminate from [0]\n",
//                                              rComm.getRank()));
//         }

//         // phase2
//         for (int i = 0; i < jobsPerProc; ++i) {
//             if ((isWaitingCmdAssignJobs[i] == true) &&
//                 (rComm.test(&(cmdAssignJobs[i][0])) == true)) {
//                 rComm.wait(&(cmdAssignJobs[i][0]));
//                 isWaitingCmdAssignJobs[i] = false;
//                 this->log_.debug(TlUtils::format(" [%d] RECV[OK] assign job[%d] from [0]\n",
//                                                  rComm.getRank(), i));
                
//                 DensityCalcJob job(this->m_nNumOfAOs,
//                                    cmdAssignJobs[i][0],
//                                    cmdAssignJobs[i][1],
//                                    cmdAssignJobs[i][2],
//                                    cmdAssignJobs[i][3]);
//                 jobList.push_back(job);
//                 continue;
//             }
//         }
//         if ((isWaitingCmdTerminateSlave == true) &&
//             (rComm.test(cmdTerminateSlave) == true)) {
//             this->log_.debug(TlUtils::format(" [%d] RECV[**] terminate from [0]\n",
//                                              rComm.getRank()));
//             rComm.wait(cmdTerminateSlave);
//             isWaitingCmdTerminateSlave = false;
            
//             isTerminate = true;
//         }

//         // phase3
//         JobListType::iterator itEnd = jobList.end();
//         for (JobListType::iterator it = jobList.begin(); it != itEnd; ) {
//             // 計算
//             const index_type start_I = it->startRow;
//             const index_type start_J = it->startCol;
//             const index_type end_I = it->endRow;
//             const index_type end_J = it->endCol;
            
//             ++numOfChances;
//             if (P.getPartialMatrixX(&(it->partialMatrix)) != true) {
//                 ++it;
//                 ++numOfMisHits;
//                 continue;
//             }
            
//             index_type* pAOs = NULL;
//             int numOfAOs = 0;
//             if (start_I != start_J) {
//                 index_type range_i = end_I - start_I;
//                 index_type range_j = end_J - start_J;
//                 numOfAOs = range_i + range_j;
//                 pAOs = new index_type[numOfAOs];
//                 for (index_type i = 0; i < range_i; ++i) {
//                     pAOs[i] = start_I + i;
//                 }
//                 for (index_type j = 0; j < range_j; ++j) {
//                     pAOs[range_i + j] = start_J + j;
//                 }
//             } else {
//                 numOfAOs = std::max(end_I, end_J) - start_I;
//                 pAOs = new index_type[numOfAOs];
//                 for (int i = 0; i < numOfAOs; ++i) {
//                     pAOs[i] = start_I + i;
//                 }
//             }
            
//             this->log_.debug(TlUtils::format(" [%d] calc...\n",
//                                              rComm.getRank()));
//             //time_calc.start();
//             this->calcPhysX(it->partialMatrix, pAOs, numOfAOs, gridInfo,
//                             pRhoAs, pGradRhoAXs, pGradRhoAYs, pGradRhoAZs);
//             //time_calc.stop();
                        
//             delete[] pAOs;
//             pAOs = NULL;
            
//             for (int i = 0; i < jobsPerProc; ++i) {
//                 if ((isWaitingCmdRequestJob[i] == true) &&
//                     (rComm.test(cmdRequestJob[i]) == true)) {
//                     rComm.wait(cmdRequestJob[i]);
//                     isWaitingCmdRequestJob[i] = false; 
//                     this->log_.debug(TlUtils::format(" [%d] SEND[OK] request job to [0]\n",
//                                                      rComm.getRank()));
                    
//                     rComm.iSendData(cmdRequestJob[i], root, TAG_REQUEST_JOB);
//                     isWaitingCmdRequestJob[i] = true;
//                     this->log_.debug(TlUtils::format(" [%d] SEND[**] request job to [0]\n",
//                                                      rComm.getRank()));
//                     break; // 1つだけ処理すればよい
//                 }
//             }
            
//             it = jobList.erase(it);
//         }
                
//         // phase4
//         if ((isTerminate == true) &&(jobList.empty() == true)) {
//             this->log_.debug(TlUtils::format(" [%d] SEND[--] terminate OK to [0]\n",
//                                              rComm.getRank()));
//             isBreakLoop = true;
//         }
//     } // end of while
        
// #ifdef _OPENMP
//     omp_set_nested(ompNestedBackup);
// #endif // _OPENMP
    
//     // statics
// //     time_all.stop();
// //     const std::size_t numOfHits = numOfChances - numOfMisHits;
// //     const double distributionEffectRatio = (double(numOfHits) / double(numOfChances)) * 100.0;
// //     this->loggerP(TlUtils::format(" [%5d] calc. time: %8.2e sec./%8.2e sec. (%6.2f%%)(dist. effect: %5.2f%%)\n",
// //                                   rComm.getRank(),
// //                                   time_calc.getElapseTime(), time_all.getElapseTime(),
// //                                   time_calc.getElapseTime() / time_all.getElapseTime() * 100.0,
// //                                   distributionEffectRatio));
    
//     // 後始末
//     for (int i = 0; i < jobsPerProc; ++i) {
//         if (isWaitingCmdAssignJobs[i] == true) {
//             rComm.cancel(&(cmdAssignJobs[i][0]));
//             isWaitingCmdAssignJobs[i] = false;
//         }
//     }
        
//     if (isWaitingCmdTerminateSlave == true) {
//         rComm.cancel(cmdTerminateSlave);
//         isWaitingCmdTerminateSlave = false;
//     }

//     for (int i = 0; i < jobsPerProc; ++i) {
//         if (isWaitingCmdRequestJob[i] == true) {
//             rComm.cancel(cmdRequestJob[i]); // 送信
//             isWaitingCmdRequestJob[i] = false;
//         }
//     }
// }


// bool DfCalcGridX_Parallel::getQueueX(index_type* pStartShell_I,
//                                      index_type* pStartShell_J,
//                                      index_type* pEndShell_I,
//                                      index_type* pEndShell_J,
//                                      double* pProgress,
//                                      bool initialize)
// {
//     const index_type numOfAOs = this->m_nNumOfAOs;
//     const int range = this->assignAoRange_;

//     const static std::size_t totalElements = numOfAOs * (numOfAOs +1) / 2;
//     static index_type currentI = 0;
//     static index_type currentJ = 0;
//     static std::size_t calcdElements = 0;
//     if (initialize == true) {
//         currentI = 0;
//         currentJ = 0;
//         calcdElements = 0;
//         return true;
//     }

    
//     while (currentI < numOfAOs) {
//         *pStartShell_I = currentI;
//         *pEndShell_I = std::min(currentI + range, numOfAOs);
        
//         while (currentJ <= currentI) {
//             *pStartShell_J = currentJ;
//             *pEndShell_J = std::min(currentJ + range, numOfAOs);

//             // for progress
//             if (pProgress != NULL) {
//                 const int width_row = *pEndShell_I - *pStartShell_I;
//                 const int width_col = *pEndShell_J - *pStartShell_J;
//                 std::size_t elements;
//                 if (currentI != currentJ) {
//                     elements = width_row * width_col;
//                 } else {
//                     elements = width_row * (width_col +1) / 2;
//                 }
//                 calcdElements += elements;
//                 *pProgress = double(calcdElements) / double(totalElements) * 100.0;
//             }
            
//             currentJ += range;

//             return true;
//         }
//         currentJ = 0;
//         currentI += range;
//     }

//     return false;
// }


// void DfCalcGridX_Parallel::calcPhysX(const TlMatrixObject& P,
//                                      const index_type* pAOs,
//                                      const std::size_t numOfAOs,
//                                      const std::vector<double>& gridInfo,
//                                      std::vector<double>* pRhoAs,
//                                      std::vector<double>* pGradRhoAXs,
//                                      std::vector<double>* pGradRhoAYs,
//                                      std::vector<double>* pGradRhoAZs)
// {
//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
    
// #pragma omp parallel for schedule(runtime)
//     for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const TlPosition position(x, y, z);

//         // calc phi table
//         std::vector<WFGrid> phis;
//         std::vector<WFGrid> gradPhiXs;
//         std::vector<WFGrid> gradPhiYs;
//         std::vector<WFGrid> gradPhiZs;
//         this->getPhiTable(position, pAOs, numOfAOs,
//                           &phis, &gradPhiXs, &gradPhiYs, &gradPhiZs);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
//         std::sort(gradPhiXs.begin(), gradPhiXs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiYs.begin(), gradPhiYs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiZs.begin(), gradPhiZs.end(), WFGrid_sort_functional());

//         // get rho at grid point
//         double rhoA, gradRhoAX, gradRhoAY, gradRhoAZ;
//         DfCalcGridX::getRhoAtGridPoint(P, phis, gradPhiXs, gradPhiYs, gradPhiZs,
//                                        &rhoA, &gradRhoAX, &gradRhoAY, &gradRhoAZ);
//         (*pRhoAs)[grid] += rhoA;
//         (*pGradRhoAXs)[grid] += gradRhoAX;
//         (*pGradRhoAYs)[grid] += gradRhoAY;
//         (*pGradRhoAZs)[grid] += gradRhoAZ;
//     }
// }


// std::vector<int>
// DfCalcGridX_Parallel::getScaLapackLocalAOs(const TlDistributeSymmetricMatrix& P) const
// {
//     std::vector<int> rowIndexTable = P.getRowIndexTable();
//     std::vector<int> colIndexTable = P.getColIndexTable();
//     std::vector<int> answer(rowIndexTable.size() + colIndexTable.size());

//     std::merge(rowIndexTable.begin(), rowIndexTable.end(),
//                colIndexTable.begin(), colIndexTable.end(),
//                answer.begin());
//     std::vector<int>::iterator endIt = std::unique(answer.begin(), answer.end());
//     answer.erase(endIt, answer.end());
    
//     return answer;
// }


// void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlDistributeSymmetricMatrix& P,
//                                              const std::vector<WFGrid>& phis,
//                                              double* pRhoA)
// {
//     const double densityCutOffValue = this->m_densityCutOffValueA;
//     double rho = 0.0;

//     // 密度行列のループ(index = p)の最大を求める
//     double value = std::sqrt(densityCutOffValue);

//     // 探査空間を狭めるための変数
//     int max_q  = phis.size();

//     std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phis.begin(), phis.end(),
//                                                                 WFGrid(0, value),
//                                                                 WFGrid_sort_functional());
//     const int max_p = std::distance(phis.begin(), pEnd);
//     for (int p = 0; p < max_p; ++p) {
//         const int AO_p = phis[p].index;
//         const double phi_p = phis[p].value;
//         const double cutValue = std::fabs(densityCutOffValue / phi_p);

//         rho += P.getLocal(AO_p, AO_p) * phi_p * phi_p;

//         std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(phis.begin() + p +1,
//                                                                     phis.begin() + max_q,
//                                                                     WFGrid(0, cutValue),
//                                                                     WFGrid_sort_functional());
//         max_q = std::distance(phis.begin(), qEnd);
//         for (int q = p + 1; q < max_q; ++q) {
//             const int AO_q = phis[q].index;
//             const double phi_q = phis[q].value;
//             rho += 2.0 * P.getLocal(AO_p, AO_q) * phi_p * phi_q;
//         }
//     }

//     *pRhoA = rho;
// }


// void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlDistributeSymmetricMatrix& P,
//                                              const std::vector<WFGrid>& phis,
//                                              const std::vector<WFGrid>& gradPhiXs,
//                                              const std::vector<WFGrid>& gradPhiYs,
//                                              const std::vector<WFGrid>& gradPhiZs,
//                                              double* pRhoA,
//                                              double* pGradRhoAX,
//                                              double* pGradRhoAY,
//                                              double* pGradRhoAZ)
// {
//     const double densityCutOffValue = this->m_densityCutOffValueA;
//     double rho = 0.0;
//     double gradRhoX = 0.0;
//     double gradRhoY = 0.0;
//     double gradRhoZ = 0.0;

//     // 密度行列のループ(index = p)の最大を求める
//     double value = std::sqrt(densityCutOffValue);
//     if (gradPhiXs.size() > 0) {
//         value = std::min(value, std::fabs(gradPhiXs[0].value));
//     }
//     if (gradPhiYs.size() > 0) {
//         value = std::min(value, std::fabs(gradPhiYs[0].value));
//     }
//     if (gradPhiZs.size() > 0) {
//         value = std::min(value, std::fabs(gradPhiZs[0].value));
//     }

//     // 探査空間を狭めるための変数
//     int max_q  = phis.size();
//     int max_qx = gradPhiXs.size();
//     int max_qy = gradPhiYs.size();
//     int max_qz = gradPhiZs.size();

//     std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(phis.begin(), phis.end(),
//                                                                 WFGrid(0, value),
//                                                                 WFGrid_sort_functional());
//     const int max_p = std::distance(phis.begin(), pEnd);
//     for (int p = 0; p < max_p; ++p) {
//         const int AO_p = phis[p].index;
//         const double phi_p = phis[p].value;
//         const double cutValue = std::fabs(densityCutOffValue / phi_p);

//         rho += P.getLocal(AO_p, AO_p) * phi_p * phi_p;

//         std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(phis.begin() + p +1,
//                                                                     phis.begin() + max_q,
//                                                                     WFGrid(0, cutValue),
//                                                                     WFGrid_sort_functional());
//         max_q = std::distance(phis.begin(), qEnd);
//         for (int q = p + 1; q < max_q; ++q) {
//             const int AO_q = phis[q].index;
//             const double phi_q = phis[q].value;
//             rho += 2.0 * P.getLocal(AO_p, AO_q) * phi_p * phi_q;
//         }

//         std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(gradPhiXs.begin(),
//                                                                      gradPhiXs.begin() + max_qx,
//                                                                      WFGrid(0, cutValue),
//                                                                      WFGrid_sort_functional());
//         max_qx = std::distance(gradPhiXs.begin(), qxEnd);
//         for (int qx = 0; qx < max_qx; ++qx) {
//             const std::size_t AO_q = gradPhiXs[qx].index;
//             const double phi_q = gradPhiXs[qx].value;
//             gradRhoX += 2.0 * P.getLocal(AO_p, AO_q) * phi_p * phi_q;
//         }

//         std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(gradPhiYs.begin(),
//                                                                      gradPhiYs.begin() + max_qy,
//                                                                      WFGrid(0, cutValue),
//                                                                      WFGrid_sort_functional());
//         max_qy = std::distance(gradPhiYs.begin(), qyEnd);
//         for (int qy = 0; qy < max_qy; ++qy) {
//             const int AO_q = gradPhiYs[qy].index;
//             const double phi_q = gradPhiYs[qy].value;
//             gradRhoY += 2.0 * P.getLocal(AO_p, AO_q) * phi_p * phi_q;
//         }

//         std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(gradPhiZs.begin(),
//                                                                      gradPhiZs.begin() + max_qz,
//                                                                      WFGrid(0, cutValue),
//                                                                      WFGrid_sort_functional());
//         max_qz = std::distance(gradPhiZs.begin(), qzEnd);
//         for (int qz = 0; qz < max_qz; ++qz) {
//             const std::size_t AO_q = gradPhiZs[qz].index;
//             const double phi_q = gradPhiZs[qz].value;
//             gradRhoZ += 2.0 * P.getLocal(AO_p, AO_q) * phi_p * phi_q;
//         }
//     }

//     *pRhoA = rho;
//     *pGradRhoAX = gradRhoX;
//     *pGradRhoAY = gradRhoY;
//     *pGradRhoAZ = gradRhoZ;
// }


// void DfCalcGridX_Parallel::addPreviousPhys(std::vector<double>* pRhoAs)
// {
//     assert(pRhoAs != NULL);

//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const std::vector<double> rhoAs = gdm.getData(atom, GridDataManager::DENSITY);
        
//         std::transform(rhoAs.begin(), rhoAs.end(),
//                        pRhoAs->begin() + atomStart,
//                        pRhoAs->begin() + atomStart, std::plus<double>());
//         atomStart += rhoAs.size();
//     }        
// }


// void DfCalcGridX_Parallel::addPreviousPhys(std::vector<double>* pRhoAs,
//                                            std::vector<double>* pRhoBs)
// {
//     assert(pRhoAs != NULL);
//     assert(pRhoBs != NULL);

//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const std::vector<double> rhoAs = gdm.getData(atom, GridDataManager::DENSITY_ALPHA);
//         const std::vector<double> rhoBs = gdm.getData(atom, GridDataManager::DENSITY_BETA);
        
//         std::transform(rhoAs.begin(), rhoAs.end(),
//                        pRhoAs->begin() + atomStart,
//                        pRhoAs->begin() + atomStart, std::plus<double>());
//         std::transform(rhoBs.begin(), rhoBs.end(),
//                        pRhoBs->begin() + atomStart,
//                        pRhoBs->begin() + atomStart, std::plus<double>());
//         atomStart += rhoAs.size();
//     }        
// }


// void DfCalcGridX_Parallel::addPreviousPhys(std::vector<double>* pRhoAs,
//                                            std::vector<double>* pGradRhoAXs,
//                                            std::vector<double>* pGradRhoAYs,
//                                            std::vector<double>* pGradRhoAZs)
// {
//     assert(pRhoAs != NULL);
//     assert(pGradRhoAXs != NULL);
//     assert(pGradRhoAYs != NULL);
//     assert(pGradRhoAZs != NULL);

//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const std::vector<double> rhoAs = gdm.getData(atom, GridDataManager::DENSITY);
//         const std::vector<double> gradRhoAXs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_X);
//         const std::vector<double> gradRhoAYs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Y);
//         const std::vector<double> gradRhoAZs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Z);
        
//         std::transform(rhoAs.begin(), rhoAs.end(),
//                        pRhoAs->begin() + atomStart,
//                        pRhoAs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAXs.begin(), gradRhoAXs.end(),
//                        pGradRhoAXs->begin() + atomStart,
//                        pGradRhoAXs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAYs.begin(), gradRhoAYs.end(),
//                        pGradRhoAYs->begin() + atomStart,
//                        pGradRhoAYs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAZs.begin(), gradRhoAZs.end(),
//                        pGradRhoAZs->begin() + atomStart,
//                        pGradRhoAZs->begin() + atomStart, std::plus<double>());
//         atomStart += rhoAs.size();
//     }        
// }


// void DfCalcGridX_Parallel::addPreviousPhys(std::vector<double>* pRhoAs,
//                                            std::vector<double>* pGradRhoAXs,
//                                            std::vector<double>* pGradRhoAYs,
//                                            std::vector<double>* pGradRhoAZs,
//                                            std::vector<double>* pRhoBs,
//                                            std::vector<double>* pGradRhoBXs,
//                                            std::vector<double>* pGradRhoBYs,
//                                            std::vector<double>* pGradRhoBZs)
// {
//     assert(pRhoAs != NULL);
//     assert(pGradRhoAXs != NULL);
//     assert(pGradRhoAYs != NULL);
//     assert(pGradRhoAZs != NULL);
//     assert(pRhoBs != NULL);
//     assert(pGradRhoBXs != NULL);
//     assert(pGradRhoBYs != NULL);
//     assert(pGradRhoBZs != NULL);

//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const std::vector<double> rhoAs = gdm.getData(atom, GridDataManager::DENSITY_ALPHA);
//         const std::vector<double> gradRhoAXs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_X_ALPHA);
//         const std::vector<double> gradRhoAYs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Y_ALPHA);
//         const std::vector<double> gradRhoAZs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Z_ALPHA);
//         const std::vector<double> rhoBs = gdm.getData(atom, GridDataManager::DENSITY_BETA);
//         const std::vector<double> gradRhoBXs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_X_BETA);
//         const std::vector<double> gradRhoBYs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Y_BETA);
//         const std::vector<double> gradRhoBZs =
//             gdm.getData(atom, GridDataManager::GRAD_DENSITY_Z_BETA);
        
//         std::transform(rhoAs.begin(), rhoAs.end(),
//                        pRhoAs->begin() + atomStart,
//                        pRhoAs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAXs.begin(), gradRhoAXs.end(),
//                        pGradRhoAXs->begin() + atomStart,
//                        pGradRhoAXs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAYs.begin(), gradRhoAYs.end(),
//                        pGradRhoAYs->begin() + atomStart,
//                        pGradRhoAYs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoAZs.begin(), gradRhoAZs.end(),
//                        pGradRhoAZs->begin() + atomStart,
//                        pGradRhoAZs->begin() + atomStart, std::plus<double>());
//         std::transform(rhoBs.begin(), rhoBs.end(),
//                        pRhoBs->begin() + atomStart,
//                        pRhoBs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoBXs.begin(), gradRhoAXs.end(),
//                        pGradRhoBXs->begin() + atomStart,
//                        pGradRhoBXs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoBYs.begin(), gradRhoAYs.end(),
//                        pGradRhoBYs->begin() + atomStart,
//                        pGradRhoBYs->begin() + atomStart, std::plus<double>());
//         std::transform(gradRhoBZs.begin(), gradRhoAZs.end(),
//                        pGradRhoBZs->begin() + atomStart,
//                        pGradRhoBZs->begin() + atomStart, std::plus<double>());
//         atomStart += rhoAs.size();
//     }        
// }


// void DfCalcGridX_Parallel::saveCurrentPhys(const std::vector<int>& gridCounts,
//                                            const std::vector<double>& rhoAs)
// {
//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     std::vector<double>::const_iterator it_rhoA = rhoAs.begin();
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const int numOfGrids = gridCounts[atom];
//         std::vector<double> buf(numOfGrids);

//         std::copy(it_rhoA, it_rhoA + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY, buf);
//         it_rhoA += numOfGrids;
//     }
// }


// void DfCalcGridX_Parallel::saveCurrentPhys(const std::vector<int>& gridCounts,
//                                            const std::vector<double>& rhoAs,
//                                            const std::vector<double>& rhoBs)
// {
//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const int numOfGrids = gridCounts[atom];
//         std::vector<double> buf(numOfGrids);

//         std::copy(rhoAs.begin() + atomStart,
//                   rhoAs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY_ALPHA, buf);

//         std::copy(rhoBs.begin() + atomStart,
//                   rhoBs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY_BETA, buf);

//         atomStart += numOfGrids;
//     }
// }


// void DfCalcGridX_Parallel::saveCurrentPhys(const std::vector<int>& gridCounts,
//                                            const std::vector<double>& rhoAs,
//                                            const std::vector<double>& gradRhoAXs,
//                                            const std::vector<double>& gradRhoAYs,
//                                            const std::vector<double>& gradRhoAZs)
// {
//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const int numOfGrids = gridCounts[atom];
//         std::vector<double> buf(numOfGrids);

//         std::copy(rhoAs.begin() + atomStart,
//                   rhoAs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY, buf);

//         std::copy(gradRhoAXs.begin() + atomStart,
//                   gradRhoAXs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_X, buf);
        
//         std::copy(gradRhoAYs.begin() + atomStart,
//                   gradRhoAYs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Y, buf);

//         std::copy(gradRhoAZs.begin() + atomStart,
//                   gradRhoAZs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Z, buf);

//         atomStart += numOfGrids;
//     }
// }


// void DfCalcGridX_Parallel::saveCurrentPhys(const std::vector<int>& gridCounts,
//                                            const std::vector<double>& rhoAs,
//                                            const std::vector<double>& gradRhoAXs,
//                                            const std::vector<double>& gradRhoAYs,
//                                            const std::vector<double>& gradRhoAZs,
//                                            const std::vector<double>& rhoBs,
//                                            const std::vector<double>& gradRhoBXs,
//                                            const std::vector<double>& gradRhoBYs,
//                                            const std::vector<double>& gradRhoBZs)
// {
//     const int numOfRealAtoms = this->numOfRealAtoms_;
//     GridDataManager gdm(this->getGridDataFilePath());
//     int atomStart = 0;
//     for (int atom = 0; atom < numOfRealAtoms; ++atom) {
//         const int numOfGrids = gridCounts[atom];
//         std::vector<double> buf(numOfGrids);

//         std::copy(rhoAs.begin() + atomStart,
//                   rhoAs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY_ALPHA, buf);

//         std::copy(gradRhoAXs.begin() + atomStart,
//                   gradRhoAXs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_X_ALPHA, buf);
        
//         std::copy(gradRhoAYs.begin() + atomStart,
//                   gradRhoAYs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Y_ALPHA, buf);

//         std::copy(gradRhoAZs.begin() + atomStart,
//                   gradRhoAZs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Z_ALPHA, buf);

//         std::copy(rhoBs.begin() + atomStart,
//                   rhoBs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::DENSITY_BETA, buf);

//         std::copy(gradRhoBXs.begin() + atomStart,
//                   gradRhoBXs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_X_BETA, buf);
        
//         std::copy(gradRhoBYs.begin() + atomStart,
//                   gradRhoBYs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Y_BETA, buf);

//         std::copy(gradRhoBZs.begin() + atomStart,
//                   gradRhoBZs.begin() + atomStart + numOfGrids,
//                   buf.begin());
//         gdm.setData(atom, GridDataManager::GRAD_DENSITY_Z_BETA, buf);

//         atomStart += numOfGrids;
//     }
// }


// double DfCalcGridX_Parallel::buildK(const std::vector<int>& gridCounts,
//                                     const std::vector<double>& gridInfo,
//                                     const std::vector<double>& rhoAs,
//                                     DfFunctional_LDA* pFunctional,
//                                     TlSparseSymmetricMatrix* pF)
// {
//     double energy = 0.0;

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
//     const std::size_t interval = (numOfAllGrids + rComm.getNumOfProc() -1) / rComm.getNumOfProc();
//     const std::size_t startGrid = interval * rComm.getRank();
//     const std::size_t endGrid = std::min(startGrid + interval, numOfAllGrids);
//     for (std::size_t grid = startGrid; grid < endGrid; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const double w = gridInfo[grid*4 +3];
//         const TlPosition position(x, y, z);
//         const double rhoA = rhoAs[grid];

//         // calc phi table
//         std::vector<WFGrid> phis;
//         this->getPhiTable(position, 0, this->m_nNumOfAOs, phis);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());

//         if (rhoA > densityCutOffValue) {
//             this->buildFock(rhoA,
//                             phis,
//                             pFunctional, w, pF);
//             energy += w * pFunctional->getFunctional(rhoA);
//         }
//     }

//     return energy;
// }


// double DfCalcGridX_Parallel::buildK(const std::vector<int>& gridCounts,
//                                     const std::vector<double>& gridInfo,
//                                     const std::vector<double>& rhoAs,
//                                     const std::vector<double>& rhoBs,
//                                     DfFunctional_LDA* pFunctional,
//                                     TlSparseSymmetricMatrix* pFA,
//                                     TlSparseSymmetricMatrix* pFB)
// {
//     double energy = 0.0;

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
//     const std::size_t interval = (numOfAllGrids + rComm.getNumOfProc() -1) / rComm.getNumOfProc();
//     const std::size_t startGrid = interval * rComm.getRank();
//     const std::size_t endGrid = std::min(startGrid + interval, numOfAllGrids);
//     for (std::size_t grid = startGrid; grid < endGrid; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const double w = gridInfo[grid*4 +3];
//         const TlPosition position(x, y, z);
//         const double rhoA = rhoAs[grid];
//         const double rhoB = rhoBs[grid];

//         // calc phi table
//         std::vector<WFGrid> phis;
//         this->getPhiTable(position, 0, this->m_nNumOfAOs, phis);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());

//         if ((rhoA > densityCutOffValue) || (rhoB > densityCutOffValue)) {
//             this->buildFock(rhoA, rhoB,
//                             phis,
//                             pFunctional, w, pFA, pFB);
//             energy += w * pFunctional->getFunctional(rhoA, rhoB);
//         }
//     }

//     return energy;
// }


// double DfCalcGridX_Parallel::buildK(const std::vector<int>& gridCounts,
//                                     const std::vector<double>& gridInfo,
//                                     const std::vector<double>& rhoAs,
//                                     const std::vector<double>& gradRhoAXs,
//                                     const std::vector<double>& gradRhoAYs,
//                                     const std::vector<double>& gradRhoAZs,
//                                     DfFunctional_GGA* pFunctional,
//                                     TlSparseSymmetricMatrix* pF)
// {
//     double energy = 0.0;

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
//     const std::size_t interval = (numOfAllGrids + rComm.getNumOfProc() -1) / rComm.getNumOfProc();
//     const std::size_t startGrid = interval * rComm.getRank();
//     const std::size_t endGrid = std::min(startGrid + interval, numOfAllGrids);
//     for (std::size_t grid = startGrid; grid < endGrid; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const double w = gridInfo[grid*4 +3];
//         const TlPosition position(x, y, z);
//         const double rhoA = rhoAs[grid];
//         const double gradRhoAX = gradRhoAXs[grid];
//         const double gradRhoAY = gradRhoAYs[grid];
//         const double gradRhoAZ = gradRhoAZs[grid];

//         // calc phi table
//         std::vector<WFGrid> phis;
//         std::vector<WFGrid> gradPhiAXs;
//         std::vector<WFGrid> gradPhiAYs;
//         std::vector<WFGrid> gradPhiAZs;
//         this->getPhiTable(position, 0, this->m_nNumOfAOs,
//                           phis, gradPhiAXs, gradPhiAYs, gradPhiAZs);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
//         std::sort(gradPhiAXs.begin(), gradPhiAXs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiAYs.begin(), gradPhiAYs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiAZs.begin(), gradPhiAZs.end(), WFGrid_sort_functional());

//         if (rhoA > densityCutOffValue) {
//             this->buildFock(rhoA, gradRhoAX, gradRhoAY, gradRhoAZ,
//                             phis, gradPhiAXs, gradPhiAYs, gradPhiAZs,
//                             pFunctional, w, pF);
//             const double gammaAA = gradRhoAX*gradRhoAX + gradRhoAY*gradRhoAY + gradRhoAZ*gradRhoAZ;
//             energy += w * pFunctional->getFunctional(rhoA, gammaAA);
//         }
//     }

//     return energy;
// }


// double DfCalcGridX_Parallel::buildK_rev2(const std::vector<int>& gridCounts,
//                                          const std::vector<double>& gridInfo,
//                                          const std::vector<double>& rhoAs,
//                                          const std::vector<double>& gradRhoAXs,
//                                          const std::vector<double>& gradRhoAYs,
//                                          const std::vector<double>& gradRhoAZs,
//                                          DfFunctional_GGA* pFunctional,
//                                          TlDistributeSymmetricMatrix* pF)
// {
//     pF->resize(this->m_nNumOfAOs);

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int numOfProcs = rComm.getNumOfProc();
//     const int rank = rComm.getRank();
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     this->loggerTime(" build the KS matrix of XC term");
//     this->logger(TlUtils::format("  assign AO range = %d\n", this->assignAoRange_));
//     this->logger(TlUtils::format("  assign jobs per proc = %d\n", this->assignJobsPerProc_));
//     //const static int JOB_PROTOCOL_SIZE = 4; // index=0: shell index
//     std::list<TlPartialSymmetricMatrix> partialFxcs;

//     if (rComm.isMaster() == true) {
//         // MASTER ==============================================================
//         this->buildK_Master(gridInfo,
//                             rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
//                             pFunctional);
//     } else {
//         // SLAVE ===============================================================
//         this->buildK_Slave(gridCounts, gridInfo,
//                            rhoAs, gradRhoAXs, gradRhoAYs, gradRhoAZs,
//                            pFunctional, &partialFxcs);
//     }
    
//     // merge
//     this->loggerTime(" reduce the KS matrix of XC term");
//     pF->mergePartialMatrix(partialFxcs);

//     // エネルギーはmasterのみで計算されている
//     // calc XC energy 
//     double energy = 0.0;
//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
//     const std::size_t interval = (numOfAllGrids / numOfProcs) +1;
//     const std::size_t startGrid = interval * rank;
//     const std::size_t endGrid = std::min(interval * (rank +1), numOfAllGrids);
//     for (std::size_t grid = startGrid; grid < endGrid; ++grid) {
//         const double w = gridInfo[grid*4 +3];
//         const double rhoA = rhoAs[grid];
//         const double gradRhoAX = gradRhoAXs[grid];
//         const double gradRhoAY = gradRhoAYs[grid];
//         const double gradRhoAZ = gradRhoAZs[grid];
//         if (rhoA > densityCutOffValue) {
//             const double gammaAA = gradRhoAX*gradRhoAX + gradRhoAY*gradRhoAY + gradRhoAZ*gradRhoAZ;
//             energy += w * pFunctional->getFunctional(rhoA, gammaAA);
//         }
//     }
//     rComm.allReduce_SUM(energy);
//     this->logger("\n");

//     return energy;
// }


// void DfCalcGridX_Parallel::buildK_Master(const std::vector<double>& gridInfo,
//                                          const std::vector<double>& rhoAs,
//                                          const std::vector<double>& gradRhoAXs,
//                                          const std::vector<double>& gradRhoAYs,
//                                          const std::vector<double>& gradRhoAZs,
//                                          DfFunctional_GGA* pFunctional)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int numOfProcs = rComm.getNumOfProc();
//     const int jobsPerProc = this->assignJobsPerProc_; // 1プロセスあたりの同時ジョブ数

//     this->getQueueX(NULL, NULL, NULL, NULL, NULL, true); // initialize only

//     // slaveからのjobのリクエスト
//     int cmdRequestJob = 0;
//     bool isWaitingCmdRequestJob = false;

//     // ジョブキュー
//     typedef std::list<int> AssignJobQueueType;
//     AssignJobQueueType assignJobQueue;

//     //
//     std::vector<bool> isWaitingJobList(numOfProcs);
//     std::vector<index_type*> jobList(numOfProcs);
//     for (int i = 0; i < numOfProcs; ++i) {
//         jobList[i] = new int[JOB_PROTOCOL_SIZE];
//     }

//     int terminateJobProc = 0;
//     bool isTerminate = false;
//     bool isBreakLoop = false;
//     int prevProgress = -1;

// #pragma omp parallel
//     {
// #pragma omp sections
//         {
// #pragma omp section
//             {
//                 while (isBreakLoop != true) {
//                     int src = 0;
//                     // phase1
//                     if (isWaitingCmdRequestJob != true) {
//                         rComm.iReceiveDataFromAnySource(cmdRequestJob, TAG_REQUEST_JOB);
//                         isWaitingCmdRequestJob = true;
//                     }
            
//                     // phase2
//                     if ((isWaitingCmdRequestJob == true) &&
//                         (rComm.test(cmdRequestJob, &src) == true)) {
//                         rComm.wait(cmdRequestJob);
//                         isWaitingCmdRequestJob = false;
//                         assignJobQueue.push_back(src);
//                     }

//                     // phase3
//                     AssignJobQueueType::iterator itEnd = assignJobQueue.end();
//                     for (AssignJobQueueType::iterator it = assignJobQueue.begin(); it != itEnd; ) {
//                         const int slave = *it;
//                         index_type startShell_I = 0;
//                         index_type startShell_J = 0;
//                         index_type endShell_I = 0;
//                         index_type endShell_J = 0;
//                         double progress = 0.0;
//                         const bool isContinued = this->getQueueX(&startShell_I, &startShell_J,
//                                                                  &endShell_I, &endShell_J,
//                                                                  &progress);
                
//                         // for progress
//                         int currentProgress = int(progress / 10.0);
//                         if (currentProgress > prevProgress) {
//                             this->loggerTime(TlUtils::format("  %3d%% done.", currentProgress * 10));
//                             prevProgress = currentProgress;
//                         }
                
//                         if (isContinued == true) {
//                             // assign job
//                             if (isWaitingJobList[slave] == true) {
//                                 rComm.wait(jobList[slave]);
//                                 isWaitingJobList[slave] = false;
//                             }
//                             jobList[slave][0] = startShell_I;
//                             jobList[slave][1] = startShell_J;
//                             jobList[slave][2] = endShell_I;
//                             jobList[slave][3] = endShell_J;
//                             rComm.iSendDataX(jobList[slave], JOB_PROTOCOL_SIZE, slave, TAG_ASSIGN_JOB);
//                             isWaitingJobList[slave] = true;
//                         } else {
//                             // finish job
//                             ++terminateJobProc;
//                             if (terminateJobProc >= (numOfProcs -1) * jobsPerProc) {
//                                 int cmdTerminateSlave = -1;
//                                 for (int i = 0; i < numOfProcs -1; ++i) {
//                                     rComm.sendData(cmdTerminateSlave, i +1, TAG_TERMINATE_SLAVE);
//                                 }
//                                 isTerminate = true;
//                             }
//                         }
//                         it = assignJobQueue.erase(it);
//                     }
            
//                     // terminate
//                     if ((isTerminate == true) &&
//                         (assignJobQueue.empty() == true)) {
//                         isBreakLoop = true;
// #pragma omp flush(isBreakLoop)
//                         {
//                         }
//                     }
//                 }
//             } // end section

// // #pragma omp section
// //             {
// //                 while (isBreakLoop != true) {
// //                 }
// //             }
//         } // end sections
//     } // end parallel

//     // 後始末
//     for (int i = 0; i < numOfProcs; ++i) {
//         if (isWaitingJobList[i] == true) {
//             rComm.cancel(jobList[i]);
//             isWaitingJobList[i] = false;
//         }
//         delete[] jobList[i];
//         jobList[i] = NULL;
//     }

//     if (isWaitingCmdRequestJob == true) {
//         rComm.cancel(cmdRequestJob);
//         isWaitingCmdRequestJob = false;
//     }
// }

// void DfCalcGridX_Parallel::buildK_Slave(const std::vector<int>& gridCounts,
//                                         const std::vector<double>& gridInfo,
//                                         const std::vector<double>& rhoAs,
//                                         const std::vector<double>& gradRhoAXs,
//                                         const std::vector<double>& gradRhoAYs,
//                                         const std::vector<double>& gradRhoAZs,
//                                         DfFunctional_GGA* pFunctional,
//                                         std::list<TlPartialSymmetricMatrix>* pFs)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     const int root = 0;
//     const int jobsPerProc = this->assignJobsPerProc_; // 1プロセスあたりの同時ジョブ数
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     // masterからのジョブの割り当て
//     std::vector<std::vector<int> > cmdAssignJobs(jobsPerProc);
//     std::vector<bool> isWaitingCmdAssignJobs(jobsPerProc);
//     for (int i = 0; i < jobsPerProc; ++i) {
//         cmdAssignJobs[i] = std::vector<int>(JOB_PROTOCOL_SIZE, 0);
//         isWaitingCmdAssignJobs[i] = false;
//     }
    
//     // masterからの終了メッセージ
//     int cmdTerminateSlave = 0;
//     bool isWaitingCmdTerminateSlave = false;
    
//     // ジョブのリスト
//     typedef std::list<DensityCalcJob> JobListType;
//     JobListType jobList;
    
//     // 初期メッセージ
//     std::vector<int> cmdRequestJob(jobsPerProc, 0);
//     std::vector<bool> isWaitingCmdRequestJob(jobsPerProc);
//     for (int i = 0; i < jobsPerProc; ++i) {
//         rComm.iSendData(cmdRequestJob[i], root, TAG_REQUEST_JOB);
//         isWaitingCmdRequestJob[i] = true;
//     }
    
//     bool isTerminate = false;
//     bool isBreakLoop = false;

//     // スレッド共通
//     JobListType::iterator it;
//     std::vector<index_type> AO_list;

// #pragma omp parallel
//     {
//         while (isBreakLoop != true) {

// #pragma omp single
//             {
//                 // phase1
//                 for (int i = 0; i < jobsPerProc; ++i) {
//                     if (isWaitingCmdAssignJobs[i] != true) {
//                         rComm.iReceiveDataX(&(cmdAssignJobs[i][0]), JOB_PROTOCOL_SIZE, root, TAG_ASSIGN_JOB);
//                         isWaitingCmdAssignJobs[i] = true;
//                     }
//                 }
//                 if (isWaitingCmdTerminateSlave != true) {
//                     rComm.iReceiveData(cmdTerminateSlave, root, TAG_TERMINATE_SLAVE);
//                     isWaitingCmdTerminateSlave = true;
//                 }

//                 // phase2
//                 for (int i = 0; i < jobsPerProc; ++i) {
//                     if ((isWaitingCmdAssignJobs[i] == true) &&
//                         (rComm.test(&(cmdAssignJobs[i][0])) == true)) {
//                         rComm.wait(&(cmdAssignJobs[i][0]));
//                         isWaitingCmdAssignJobs[i] = false;
                        
//                         DensityCalcJob job(this->m_nNumOfAOs,
//                                            cmdAssignJobs[i][0],
//                                            cmdAssignJobs[i][1],
//                                            cmdAssignJobs[i][2],
//                                            cmdAssignJobs[i][3]);
//                         jobList.push_back(job);
//                     }
//                 }
//                 if ((isWaitingCmdTerminateSlave == true) &&
//                     (rComm.test(cmdTerminateSlave) == true)) {
//                     this->log_.debug(TlUtils::format(" [%d] BUILD_XC: RECV[**] terminate from [0]\n",
//                                                      rComm.getRank()));
//                     rComm.wait(cmdTerminateSlave);
//                     isWaitingCmdTerminateSlave = false;
                    
//                     isTerminate = true;
//                 }
//             } // end single

//             {
//                 // phase3
//                 JobListType::iterator itEnd = jobList.end();
//                 for (it = jobList.begin(); it != itEnd; ) {
//                     // 計算条件設定
//                     const index_type startRow = it->startRow;
//                     const index_type startCol = it->startCol;
//                     const index_type endRow = it->endRow;
//                     const index_type endCol = it->endCol;
//                     const index_type rangeRow = endRow - startRow;
//                     const index_type rangeCol = endCol - startCol;
                    
//                     // AOリストの作成
// #pragma omp single
//                     {
//                         if (startRow != startCol) {
//                             AO_list.resize(rangeRow + rangeCol);
//                             for (index_type i = 0; i < rangeRow; ++i) {
//                                 AO_list[i] = startRow + i;
//                             }
//                             for (index_type i = 0; i < rangeCol; ++i) {
//                                 AO_list[rangeRow + i] = startCol + i;
//                             }
//                         } else {
//                             assert(rangeRow == rangeCol);
//                             AO_list.resize(rangeRow);
//                             for (index_type i = 0; i < rangeRow; ++i) {
//                                 AO_list[i] = startRow + i;
//                             }
//                         }
//                     }

//                     // XC行列要素を求め、部分行列へ
//                     const std::size_t numOfAllGrids = gridInfo.size() / 4;
// #pragma omp for schedule(runtime)
//                     for (std::size_t grid = 0; grid < numOfAllGrids; ++grid) {
//                         const double x = gridInfo[grid*4   ];
//                         const double y = gridInfo[grid*4 +1];
//                         const double z = gridInfo[grid*4 +2];
//                         const double w = gridInfo[grid*4 +3];
//                         const TlPosition position(x, y, z);
//                         const double rhoA = rhoAs[grid];
//                         const double gradRhoAX = gradRhoAXs[grid];
//                         const double gradRhoAY = gradRhoAYs[grid];
//                         const double gradRhoAZ = gradRhoAZs[grid];
                        
//                         if (rhoA > densityCutOffValue) {
//                             // calc phi table
//                             std::vector<WFGrid> phis;
//                             std::vector<WFGrid> gradPhiAXs;
//                             std::vector<WFGrid> gradPhiAYs;
//                             std::vector<WFGrid> gradPhiAZs;
//                             this->getPhiTable(position, &(AO_list[0]), AO_list.size(),
//                                               &phis, &gradPhiAXs, &gradPhiAYs, &gradPhiAZs);
//                             std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
//                             std::sort(gradPhiAXs.begin(), gradPhiAXs.end(), WFGrid_sort_functional());
//                             std::sort(gradPhiAYs.begin(), gradPhiAYs.end(), WFGrid_sort_functional());
//                             std::sort(gradPhiAZs.begin(), gradPhiAZs.end(), WFGrid_sort_functional());
                            
// #pragma omp critical (DfCalcGridX_Parallel__buildK_Slave)
//                             {
//                                 this->buildFock(rhoA, gradRhoAX, gradRhoAY, gradRhoAZ,
//                                                 phis, gradPhiAXs, gradPhiAYs, gradPhiAZs,
//                                                 pFunctional, w, &(it->partialMatrix));
//                             }
//                             // エネルギーは既にmasterで計算されている
//                             //const double gammaAA = gradRhoAX*gradRhoAX + gradRhoAY*gradRhoAY + gradRhoAZ*gradRhoAZ;
//                             //energy += w * pFunctional->getFunctional(rhoA, gammaAA);
//                         }
//                     }

// #pragma omp single
//                     {
//                         pFs->push_back(it->partialMatrix);
            
//                         // 計算終了処理
//                         for (int i = 0; i < jobsPerProc; ++i) {
//                             if ((isWaitingCmdRequestJob[i] == true) &&
//                                 (rComm.test(cmdRequestJob[i]) == true)) {
//                                 rComm.wait(cmdRequestJob[i]);
//                                 isWaitingCmdRequestJob[i] = false; 
                                
//                                 rComm.iSendData(cmdRequestJob[i], root, TAG_REQUEST_JOB);
//                                 isWaitingCmdRequestJob[i] = true;
//                             }
//                             break; // 1つ処理すればよい
//                         }
                    
//                         it = jobList.erase(it);
//                     }
//                 }
//             }

// #pragma omp single
//             {
//                 // phase4
//                 if ((isTerminate == true) &&(jobList.empty() == true)) {
//                     isBreakLoop = true; // break loop
//                 }
//             }
//         } // end while
//     }
    
//     // 後始末
//     for (int i = 0; i < jobsPerProc; ++i) {
//         if (isWaitingCmdAssignJobs[i] == true) {
//             rComm.cancel(&(cmdAssignJobs[i][0]));
//             isWaitingCmdAssignJobs[i] = false;
//         }
//     }
//     if (isWaitingCmdTerminateSlave == true) {
//         rComm.cancel(cmdTerminateSlave);
//         isWaitingCmdTerminateSlave = false;
//     }
//     for (int i = 0; i < jobsPerProc; ++i) {
//         if (isWaitingCmdRequestJob[i] == true) {
//             rComm.cancel(cmdRequestJob[i]); // 送信
//             isWaitingCmdRequestJob[i] = false;
//         }
//     }
// }


// double DfCalcGridX_Parallel::buildK(const std::vector<int>& gridCounts,
//                                     const std::vector<double>& gridInfo,
//                                     const std::vector<double>& rhoAs,
//                                     const std::vector<double>& gradRhoAXs,
//                                     const std::vector<double>& gradRhoAYs,
//                                     const std::vector<double>& gradRhoAZs,
//                                     const std::vector<double>& rhoBs,
//                                     const std::vector<double>& gradRhoBXs,
//                                     const std::vector<double>& gradRhoBYs,
//                                     const std::vector<double>& gradRhoBZs,
//                                     DfFunctional_GGA* pFunctional,
//                                     TlSparseSymmetricMatrix* pFA,
//                                     TlSparseSymmetricMatrix* pFB)
// {
//     double energy = 0.0;

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const double densityCutOffValue = this->m_densityCutOffValueA;

//     const std::size_t numOfAllGrids = gridInfo.size() / 4;
//     const std::size_t interval = (numOfAllGrids + rComm.getNumOfProc() -1) / rComm.getNumOfProc();
//     const std::size_t startGrid = interval * rComm.getRank();
//     const std::size_t endGrid = std::min(startGrid + interval, numOfAllGrids);
//     for (std::size_t grid = startGrid; grid < endGrid; ++grid) {
//         const double x = gridInfo[grid*4   ];
//         const double y = gridInfo[grid*4 +1];
//         const double z = gridInfo[grid*4 +2];
//         const double w = gridInfo[grid*4 +3];
//         const TlPosition position(x, y, z);
//         const double rhoA = rhoAs[grid];
//         const double gradRhoAX = gradRhoAXs[grid];
//         const double gradRhoAY = gradRhoAYs[grid];
//         const double gradRhoAZ = gradRhoAZs[grid];
//         const double rhoB = rhoBs[grid];
//         const double gradRhoBX = gradRhoBXs[grid];
//         const double gradRhoBY = gradRhoBYs[grid];
//         const double gradRhoBZ = gradRhoBZs[grid];
        
//         // calc phi table
//         std::vector<WFGrid> phis;
//         std::vector<WFGrid> gradPhiXs;
//         std::vector<WFGrid> gradPhiYs;
//         std::vector<WFGrid> gradPhiZs;
//         this->getPhiTable(position, 0, this->m_nNumOfAOs,
//                           phis, gradPhiXs, gradPhiYs, gradPhiZs);
//         std::sort(phis.begin(), phis.end(), WFGrid_sort_functional());
//         std::sort(gradPhiXs.begin(), gradPhiXs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiYs.begin(), gradPhiYs.end(), WFGrid_sort_functional());
//         std::sort(gradPhiZs.begin(), gradPhiZs.end(), WFGrid_sort_functional());
        
//         if ((rhoA > densityCutOffValue) || (rhoB > densityCutOffValue)){
//             this->buildFock(rhoA, gradRhoAX, gradRhoAY, gradRhoAZ,
//                             rhoB, gradRhoBX, gradRhoBY, gradRhoBZ,                                
//                             phis, gradPhiXs, gradPhiYs, gradPhiZs,
//                             pFunctional, w, pFA, pFB);
//             const double gammaAA = gradRhoAX*gradRhoAX + gradRhoAY*gradRhoAY + gradRhoAZ*gradRhoAZ;
//             const double gammaAB = gradRhoAX*gradRhoBX + gradRhoAY*gradRhoBY + gradRhoAZ*gradRhoBZ;
//             const double gammaBB = gradRhoBX*gradRhoBX + gradRhoBY*gradRhoBY + gradRhoBZ*gradRhoBZ;
//             energy += w * pFunctional->getFunctional(rhoA, rhoB, gammaAA, gammaAB, gammaBB);
//         }
//     }

//     return energy;
// }

// experimental code -----------------------------------------------------------
// double DfCalcGridX_Parallel::calcXC_DC(const TlDistributeSymmetricMatrix& P,
//                                        DfFunctional_GGA* pFunctional,
//                                        TlDistributeSymmetricMatrix* pF)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(pFunctional != NULL);
//     assert(pF != NULL);

//     // setup
//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     this->physicalValues_.clear();
//     this->defineCutOffValues(P);

//     const TlMatrix localP = 0.5 * TlDistributeMatrix(P).getLocalMatrix();
//     const std::vector<index_type> rowIndexes = P.getRowIndexTable();
//     const std::vector<index_type> colIndexes = P.getColIndexTable();
//     const index_type numOfLocalRows = rowIndexes.size();
//     const index_type numOfLocalCols = colIndexes.size();
//     assert(localP.getNumOfRows() == numOfLocalRows);
//     assert(localP.getNumOfCols() == numOfLocalCols);
    
    
//     TlMatrix gridMatrix = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
//     gridMatrix.resize(gridMatrix.getNumOfRows(), 8);

//     //
//     this->loggerTime(" calculate population of each grid.");
//     TlMatrix prevRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     this->calcRhoVals_GGA(rowIndexes, colIndexes, localP,
//                           &gridMatrix);
//     TlMatrix currRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     currRhoValsMtx -= prevRhoValsMtx;
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));
//     rComm.allReduce_SUM(currRhoValsMtx);
//     currRhoValsMtx += prevRhoValsMtx;
//     gridMatrix.setBlockMatrix(0, 4, currRhoValsMtx);

//     // 
//     this->loggerTime(" build the KS matrix of XC term");
//     TlMatrix localF(localP.getNumOfRows(), localP.getNumOfCols());
//     const double energy = this->buildK_2(gridMatrix,
//                                          rowIndexes, colIndexes,
//                                          pFunctional, &localF);
//     for (index_type r = 0; r < numOfLocalRows; ++r) {
//         const index_type globalRow = rowIndexes[r];
//         for (index_type c = 0; c < numOfLocalCols; ++c) {
//             const index_type globalCol = colIndexes[c];
//             if (globalRow >= globalCol) {
//                 const double value = localF.get(r, c);
//                 pF->add(globalRow, globalCol, value);
//             }
//         }
//     }
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));

//     if (rComm.isMaster() == true) {
//         if (this->m_bIsUpdateXC == true) {
//             this->loggerTime(" save grid matrix");
//             this->saveGridMatrix(gridMatrix);
//         }
//     }
//     this->loggerTime(" finished");

//     return energy;
// }


// double DfCalcGridX_Parallel::calcXC_DC(const TlDistributeSymmetricMatrix& P_A,
//                                        DfFunctional_GGA* pFunctional,
//                                        TlDistributeSymmetricMatrix* pF_A)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(pFunctional != NULL);
//     assert(pF_A != NULL);
    
//     //this->defineCutOffValues(P_A, P_B);
//     this->defineCutOffValues(P_A);

//     // 大域分散行列である密度行列において、ローカルな部分行列を取得する
//     const TlMatrix localP_A = TlDistributeMatrix(P_A).getLocalMatrix();
//     const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
//     const std::vector<index_type> colIndexes = P_A.getColIndexTable();
//     const index_type numOfLocalRows = rowIndexes.size();
//     const index_type numOfLocalCols = colIndexes.size();
//     assert(localP_A.getNumOfRows() == numOfLocalRows);
//     assert(localP_A.getNumOfCols() == numOfLocalCols);
    
//     TlMatrix gridMatrix = this->getGridMatrix();
//     gridMatrix.resize(gridMatrix.getNumOfRows(), 8);

//     this->log_.info(" calculate population of each grid.");
//     TlMatrix prevRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     this->calcRhoVals_GGA(rowIndexes, colIndexes, localP_A,
//                           &gridMatrix);
//     TlMatrix currRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     currRhoValsMtx -= prevRhoValsMtx;
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));
//     rComm.allReduce_SUM(currRhoValsMtx);
//     currRhoValsMtx += prevRhoValsMtx;
//     gridMatrix.setBlockMatrix(0, 4, currRhoValsMtx);

//     this->log_.info(" build the KS matrix of XC term");
//     TlMatrix localF_A(localP_A.getNumOfRows(), localP_A.getNumOfCols());
//     const double energy = this->buildK_2(gridMatrix,
//                                          rowIndexes, colIndexes,
//                                          pFunctional, &localF_A);
//     for (index_type r = 0; r < numOfLocalRows; ++r) {
//         const index_type globalRow = rowIndexes[r];
//         for (index_type c = 0; c < numOfLocalCols; ++c) {
//             const index_type globalCol = colIndexes[c];
//             if (globalRow >= globalCol) {
//                 const double value = localF_A.get(r, c);
//                 pF_A->add(globalRow, globalCol, value);
//             }
//         }
//     }
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));

//     if (rComm.isMaster() == true) {
//         if (this->m_bIsUpdateXC == true) {
//             this->loggerTime(" save grid matrix");
//             this->saveGridMatrix(gridMatrix);
//         }
//     }
//     this->loggerTime(" finished");

//     return energy;
// }


// void DfCalcGridX_Parallel::calcRhoVals_GGA(const std::vector<index_type>& P_rowIndexes,
//                                            const std::vector<index_type>& P_colIndexes,
//                                            const TlDistributeMatrix& P,
//                                            TlMatrix* pGridMatrix)
// {
//     assert(pGridMatrix != NULL);
    
//     const std::size_t numOfRows = P_rowIndexes.size();
//     const std::size_t numOfCols = P_colIndexes.size();
    
//     const std::size_t numOfGrids = pGridMatrix->getNumOfRows();
// #pragma omp parallel for schedule(runtime)
//     for (std::size_t grid = 0; grid < numOfGrids; ++grid) {
//         const double x = pGridMatrix->get(grid, 0);
//         const double y = pGridMatrix->get(grid, 1);
//         const double z = pGridMatrix->get(grid, 2);
//         const TlPosition r(x, y, z);
        
//         TlMatrix wf_row;
//         this->getWaveFunctionValues(P_rowIndexes, r, &wf_row);

//         TlMatrix wf_col, wf_dX, wf_dY, wf_dZ;
//         this->getWaveFunctionValues(P_colIndexes, r,
//                                     &wf_col, &wf_dX, &wf_dY, &wf_dZ);
//         wf_col.transpose();
//         wf_dX.transpose();
//         wf_dY.transpose();
//         wf_dZ.transpose();
//         {
//             TlMatrix wf_rc = wf_row * wf_col;
//             assert(wf_rc.getNumOfRows() == numOfRows);
//             assert(wf_rc.getNumOfCols() == numOfCols);
//             wf_rc.dot(P);
//             pGridMatrix->add(grid, 4, wf_rc.sum());
//         }
//         {
//             TlMatrix wf_rc = wf_row * wf_dX;
//             assert(wf_rc.getNumOfRows() == numOfRows);
//             assert(wf_rc.getNumOfCols() == numOfCols);
//             wf_rc.dot(P);
//             pGridMatrix->add(grid, 5, 2.0 * wf_rc.sum());
//         }
//         {
//             TlMatrix wf_rc = wf_row * wf_dY;
//             assert(wf_rc.getNumOfRows() == numOfRows);
//             assert(wf_rc.getNumOfCols() == numOfCols);
//             wf_rc.dot(P);
//             pGridMatrix->add(grid, 6, 2.0 * wf_rc.sum());
//         }
//         {
//             TlMatrix wf_rc = wf_row * wf_dZ;
//             assert(wf_rc.getNumOfRows() == numOfRows);
//             assert(wf_rc.getNumOfCols() == numOfCols);
//             wf_rc.dot(P);
//             pGridMatrix->add(grid, 7, 2.0 * wf_rc.sum());
//         }
//     }
// }

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy_DC(const TlSymmetricMatrix& P_A,
//                                                             DfFunctional_LDA* pFunctional,
//                                                             TlSymmetricMatrix* pF_A)
// {
//     this->calcRho_LDA(P_A);
//     double energy = this->buildVxc(pFunctional,
//                                    pF_A);
//     return energy;
// }

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlSymmetricMatrix* pF_A)
{
    this->calcRho_LDA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         const TlSymmetricMatrix& P_B,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlSymmetricMatrix* pF_A,
                                                         TlSymmetricMatrix* pF_B)
{
    this->calcRho_LDA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlSymmetricMatrix* pF_A)
{
    this->calcRho_GGA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                                         const TlSymmetricMatrix& P_B,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlSymmetricMatrix* pF_A,
                                                         TlSymmetricMatrix* pF_B)
{
    this->calcRho_GGA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_LDA(P_A, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_LDA(P_A, P_B, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_GGA(P_A, &gridMat);
    this->gatherGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration -1);
    DfCalcGridX::calcRho_GGA(P_A, P_B, &gridMat);
    this->gatherGridMatrix(gridMat);
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlSymmetricMatrix* pF_A,
                                      TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlSymmetricMatrix* pF_A)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlSymmetricMatrix* pF_A,
                                      TlSymmetricMatrix* pF_B)
{
    TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, pF_A, pF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

TlMatrix DfCalcGridX_Parallel::distributeGridMatrix(const int iteration)
{
    this->log_.info("distribute grid matrix: start");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    const int tag = TAG_CALC_GRID_DISTRIBUTE;
    TlMatrix gridMat;
    if (rComm.isMaster() == true) {
#ifdef USE_FILE_MATRIX
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(iteration));
#else
        TlMatrix globalGridMat;
        globalGridMat.load(DfObject::getGridMatrixPath(iteration));
#endif // USE_FILE_MATRIX
        this->numOfRows_gridMatrix_ = globalGridMat.getNumOfRows();
        this->numOfCols_gridMatrix_ = globalGridMat.getNumOfCols();
        
        const index_type range = (this->numOfRows_gridMatrix_ + numOfProcs -1) / numOfProcs;

        // for self(master)
        {
            const index_type startGrid = 0;
            const index_type endGrid = std::min<index_type>(startGrid + range,
                                                            this->numOfRows_gridMatrix_);
            gridMat = globalGridMat.getBlockMatrix(startGrid, 0,
                                                   endGrid - startGrid,
                                                   globalGridMat.getNumOfCols());
        }
        
        // for slave
        for (int i = 1; i < numOfProcs; ++i) {
            const index_type startGrid = range * i;
            const index_type endGrid = std::min<index_type>(startGrid + range,
                                                            this->numOfRows_gridMatrix_);
            TlMatrix tmpMat = globalGridMat.getBlockMatrix(startGrid, 0,
                                                           endGrid - startGrid,
                                                           globalGridMat.getNumOfCols());
            rComm.sendData(tmpMat, i, tag);
            this->log_.info(TlUtils::format("send grid data to %d", i));
        }
    } else {
      rComm.receiveData(gridMat, 0, tag);
      this->log_.info("recv grid data");
    }

    this->log_.info("distribute grid matrix: end");
    return gridMat;
}

void DfCalcGridX_Parallel::gatherGridMatrix(const TlMatrix& gridMat)
{
  this->log_.info("gather grid matrix: start");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    const int tag = TAG_CALCGRID_GATHER;

    if (rComm.isMaster() == true) {
#ifdef USE_FILE_MATRIX
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(this->m_nIteration),
                                   this->numOfRows_gridMatrix_,
                                   this->numOfCols_gridMatrix_);
#else
        TlMatrix globalGridMat(this->numOfRows_gridMatrix_,
                               this->numOfCols_gridMatrix_);
#endif // USE_FILE_MATRIX
        globalGridMat.setBlockMatrix(0, 0,
                                     gridMat);
        index_type currentGridIndex = gridMat.getNumOfRows();
        std::vector<bool> recvCheck(numOfProcs, false);
        for (int i = 1; i < numOfProcs; ++i) {
	    int proc = 0;
            TlMatrix tmpMat;
            rComm.receiveDataFromAnySource(tmpMat, &proc, tag);
	    if (recvCheck[proc] != false) {
	      this->log_.warn(TlUtils::format("already receive grid data from %d",
					      proc));
	    }
	    recvCheck[proc] = true;
	    this->log_.debug(TlUtils::format("recv grid data from %d", proc));

            assert(globalGridMat.getNumOfCols() == gridMat.getNumOfCols());
            globalGridMat.setBlockMatrix(currentGridIndex, 0,
                                         tmpMat);
            currentGridIndex += gridMat.getNumOfRows();
        }

#ifndef USE_FILE_MATRIX
        globalGridMat.save(DfObject::getGridMatrixPath(this->m_nIteration));
#endif // USE_FILE_MATRIX
    } else {
      rComm.sendData(gridMat, 0, tag);
    }
  this->log_.info("gather grid matrix: end");
}

// =============================================================================
// double DfCalcGridX_Parallel::calcXC_DC(const TlDistributeSymmetricMatrix& P,
//                                        DfFunctional_GGA* pFunctional,
//                                        TlDistributeSymmetricMatrix* pF)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(pFunctional != NULL);
//     assert(pF != NULL);

//     // setup
//     if (rComm.isMaster() == true) {
//         DfCalcGridX::backupGridData();
//     }
//     this->physicalValues_.clear();
//     this->defineCutOffValues(P);

//     const TlMatrix localP = 0.5 * TlDistributeMatrix(P).getLocalMatrix();
//     const std::vector<index_type> rowIndexes = P.getRowIndexTable();
//     const std::vector<index_type> colIndexes = P.getColIndexTable();
//     const index_type numOfLocalRows = rowIndexes.size();
//     const index_type numOfLocalCols = colIndexes.size();
//     assert(localP.getNumOfRows() == numOfLocalRows);
//     assert(localP.getNumOfCols() == numOfLocalCols);
    
    
//     TlMatrix gridMatrix = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration -1);
//     gridMatrix.resize(gridMatrix.getNumOfRows(), 8);

//     //
//     this->loggerTime(" calculate population of each grid.");
//     TlMatrix prevRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     this->calcRhoVals_GGA(rowIndexes, colIndexes, localP,
//                           &gridMatrix);
//     TlMatrix currRhoValsMtx =
//         gridMatrix.getBlockMatrix(0, 4,
//                                   gridMatrix.getNumOfRows(), 4);
//     currRhoValsMtx -= prevRhoValsMtx;
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));
//     rComm.allReduce_SUM(currRhoValsMtx);
//     currRhoValsMtx += prevRhoValsMtx;
//     gridMatrix.setBlockMatrix(0, 4, currRhoValsMtx);

//     // 
//     this->loggerTime(" build the KS matrix of XC term");
//     TlMatrix localF(localP.getNumOfRows(), localP.getNumOfCols());
//     const double energy = this->buildK_2(gridMatrix,
//                                          rowIndexes, colIndexes,
//                                          pFunctional, &localF);
//     for (index_type r = 0; r < numOfLocalRows; ++r) {
//         const index_type globalRow = rowIndexes[r];
//         for (index_type c = 0; c < numOfLocalCols; ++c) {
//             const index_type globalCol = colIndexes[c];
//             if (globalRow >= globalCol) {
//                 const double value = localF.get(r, c);
//                 pF->add(globalRow, globalCol, value);
//             }
//         }
//     }
//     this->log_.info(TlUtils::format("  waiting...#%6d", rComm.getRank()));

//     if (rComm.isMaster() == true) {
//         if (this->m_bIsUpdateXC == true) {
//             this->loggerTime(" save grid matrix");
//             this->saveGridMatrix(gridMatrix);
//         }
//     }
//     this->loggerTime(" finished");

//     return energy;
// }

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A)
{
    this->calcRho_LDA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         const TlDistributeSymmetricMatrix& P_B,
                                                         DfFunctional_LDA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A,
                                                         TlDistributeSymmetricMatrix* pF_B)
{
    this->calcRho_LDA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A)
{
    this->calcRho_GGA(P_A);
    double energy = this->buildVxc(pFunctional, pF_A);
    return energy;
}

double DfCalcGridX_Parallel::calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                                         const TlDistributeSymmetricMatrix& P_B,
                                                         DfFunctional_GGA* pFunctional,
                                                         TlDistributeSymmetricMatrix* pF_A,
                                                         TlDistributeSymmetricMatrix* pF_B)
{
    this->calcRho_GGA(P_A, P_B);
    double energy = this->buildVxc(pFunctional, pF_A, pF_B);
    return energy;
}

TlMatrix DfCalcGridX_Parallel::getGlobalGridMatrix(const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlMatrix gridMat;
    TlMatrix rhoMat;
    index_type numOfGrids = 0;
    index_type numOfCols = 0;
    if (rComm.isMaster() == true) {
        gridMat = DfObject::getGridMatrix<TlMatrix>(iteration);
        numOfGrids = gridMat.getNumOfRows();
        numOfCols = gridMat.getNumOfCols();

        TlMatrix crdMat = gridMat.getBlockMatrix(0, 0,
                                                 numOfGrids,
                                                 GM_ATOM_INDEX +1);
        rhoMat = gridMat.getBlockMatrix(0, GM_ATOM_INDEX +1,
                                        numOfGrids,
                                        numOfCols - (GM_ATOM_INDEX +1));
        gridMat = crdMat;
    }
    rComm.broadcast(numOfGrids);
    rComm.broadcast(numOfCols);
    rComm.broadcast(gridMat);
    gridMat.resize(numOfGrids, numOfCols);

    if (this->m_bIsUpdateXC == true) {
        if (rComm.isMaster() == true) {
            gridMat.setBlockMatrix(0, GM_ATOM_INDEX +1,
                                   rhoMat);
        }
    }
    
    return gridMat;
}

void DfCalcGridX_Parallel::allReduceGridMatrix(const TlMatrix& gridMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const index_type numOfAllGrids = gridMat.getNumOfRows();
    const index_type numOfCols = gridMat.getNumOfCols();
    
    TlMatrix rhoMat = gridMat.getBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                             gridMat.getNumOfRows(),
                                             gridMat.getNumOfCols() - GM_LDA_RHO_ALPHA);
    rComm.allReduce_SUM(rhoMat);

    if (rComm.isMaster() == true) {
        TlFileMatrix globalGridMat(DfObject::getGridMatrixPath(this->m_nIteration),
                                   gridMat.getNumOfRows(),
                                   gridMat.getNumOfCols());
        assert(globalGridMat.getNumOfRows() == numOfAllGrids);
        assert(globalGridMat.getNumOfCols() == numOfCols);
        
        const TlMatrix crdMat = gridMat.getBlockMatrix(0, 0,
                                                       gridMat.getNumOfRows(),
                                                       GM_LDA_RHO_ALPHA);
        globalGridMat.setBlockMatrix(0, 0,
                                     crdMat);
        globalGridMat.setBlockMatrix(0, GM_LDA_RHO_ALPHA,
                                     rhoMat);
    }
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_LDA(TlDistributeMatrix(P_A), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_LDA(TlDistributeMatrix(P_A),
                      TlDistributeMatrix(P_B), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeSymmetricMatrix& P_A)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_GGA(TlDistributeMatrix(P_A), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B)
{
    TlMatrix gridMat = this->getGlobalGridMatrix(this->m_nIteration -1);
    this->calcRho_GGA(TlDistributeMatrix(P_A),
                      TlDistributeMatrix(P_B), &gridMat);
    this->allReduceGridMatrix(gridMat);
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeMatrix& P_A,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi);

        // get rho at grid point
        double rhoA = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi, colPhi,
                                &rhoA);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
    }
}

void DfCalcGridX_Parallel::calcRho_LDA(const TlDistributeMatrix& P_A,
                                       const TlDistributeMatrix& P_B,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const TlMatrix localP_B = P_B.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();
    // rowIndexes, colIndexes はP_Bと共通

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi);

        // get rho at grid point
        double rhoA = 0.0;
        double rhoB = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi, colPhi,
                                &rhoA);
        this->getRhoAtGridPoint(localP_B,
                                rowPhi, colPhi,
                                &rhoB);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_RHO_BETA, rhoB);
    }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeMatrix& P_A,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        std::vector<WFGrid> colGradPhiX;
        std::vector<WFGrid> colGradPhiY;
        std::vector<WFGrid> colGradPhiZ;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi, &colGradPhiX, &colGradPhiY, &colGradPhiZ);

        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi,
                                colPhi, colGradPhiX, colGradPhiY, colGradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoXA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoYA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoZA);
    }
}

void DfCalcGridX_Parallel::calcRho_GGA(const TlDistributeMatrix& P_A,
                                       const TlDistributeMatrix& P_B,
                                       TlMatrix* pGridMat)
{
    const TlMatrix localP_A = P_A.getLocalMatrix();
    const TlMatrix localP_B = P_B.getLocalMatrix();
    const std::vector<index_type> rowIndexes = P_A.getRowIndexTable();
    const std::vector<index_type> colIndexes = P_A.getColIndexTable();
    // rowIndexes, colIndexes はP_Bと共通

    const index_type numOfGrids = pGridMat->getNumOfRows();
#pragma omp parallel for schedule(runtime)
    for (index_type grid = 0; grid < numOfGrids; ++grid) {
        const TlPosition gridPosition(pGridMat->get(grid, GM_X),
                                      pGridMat->get(grid, GM_Y),
                                      pGridMat->get(grid, GM_Z));

        // calc phi table
        std::vector<WFGrid> rowPhi;
        std::vector<WFGrid> colPhi;
        std::vector<WFGrid> colGradPhiX;
        std::vector<WFGrid> colGradPhiY;
        std::vector<WFGrid> colGradPhiZ;
        this->getPhiTable(gridPosition, rowIndexes,
                          &rowPhi);
        this->getPhiTable(gridPosition, colIndexes,
                          &colPhi, &colGradPhiX, &colGradPhiY, &colGradPhiZ);

        // get rho at grid point
        double rhoA = 0.0;
        double gradRhoXA = 0.0;
        double gradRhoYA = 0.0;
        double gradRhoZA = 0.0;
        double rhoB = 0.0;
        double gradRhoXB = 0.0;
        double gradRhoYB = 0.0;
        double gradRhoZB = 0.0;
        this->getRhoAtGridPoint(localP_A,
                                rowPhi,
                                colPhi, colGradPhiX, colGradPhiY, colGradPhiZ,
                                &rhoA, &gradRhoXA, &gradRhoYA, &gradRhoZA);
        this->getRhoAtGridPoint(localP_B,
                                rowPhi,
                                colPhi, colGradPhiY, colGradPhiY, colGradPhiZ,
                                &rhoB, &gradRhoXB, &gradRhoYB, &gradRhoZB);

        pGridMat->add(grid, GM_GGA_RHO_ALPHA, rhoA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_ALPHA, gradRhoXA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_ALPHA, gradRhoYA);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_ALPHA, gradRhoZA);
        pGridMat->add(grid, GM_GGA_RHO_BETA, rhoB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_X_BETA, gradRhoXB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Y_BETA, gradRhoYB);
        pGridMat->add(grid, GM_GGA_GRAD_RHO_Z_BETA, gradRhoZB);
    }
}

void DfCalcGridX_Parallel::getPhiTable(const TlPosition& gridPosition,
                                       const std::vector<index_type>& AO_List,
                                       std::vector<WFGrid>* pPhis)
{
    assert(pPhis != NULL);
    double densityCutOffValue = 1.0E-16;
    
    const std::size_t AO_ListSize = AO_List.size();
    pPhis->clear();
    pPhis->reserve(AO_ListSize);

    // orbital loop
    for (std::size_t AO_index = 0; AO_index < AO_ListSize; ++AO_index) {
        const index_type AO = AO_List[AO_index];

        double phi = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = DfCalcGridX::getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            double prefactorX = 0.0;
            double prefactorY = 0.0;
            double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = prefactor
                    * this->m_tlOrbInfo.getCoefficient(AO, PGTO)
                    * std::exp(-shoulder);
                phi += g;
            }
        }
        
        if (std::fabs(phi) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, phi);
            const WFGrid wfGrid(AO_index, phi);
            pPhis->push_back(wfGrid);
        }
    }

    std::sort(pPhis->begin(), pPhis->end(), WFGrid_sort_functional());
}

void DfCalcGridX_Parallel::getPhiTable(const TlPosition& gridPosition,
                                       const std::vector<index_type>& AO_List,
                                       std::vector<WFGrid>* pPhis,
                                       std::vector<WFGrid>* pGradPhiXs,
                                       std::vector<WFGrid>* pGradPhiYs,
                                       std::vector<WFGrid>* pGradPhiZs)
{
    assert(pPhis != NULL);
    assert(pGradPhiXs != NULL);
    assert(pGradPhiYs != NULL);
    assert(pGradPhiZs != NULL);
    const double densityCutOffValue = 1.0E-16;
    
    const std::size_t AO_ListSize = AO_List.size();
    pPhis->clear();
    pPhis->reserve(AO_ListSize);
    pGradPhiXs->clear();
    pGradPhiXs->reserve(AO_ListSize);
    pGradPhiYs->clear();
    pGradPhiYs->reserve(AO_ListSize);
    pGradPhiZs->clear();
    pGradPhiZs->reserve(AO_ListSize);

    // orbital loop
    for (std::size_t AO_index = 0; AO_index < AO_ListSize; ++AO_index) {
        const index_type AO = AO_List[AO_index];

        double phi = 0.0;
        double gradPhiX = 0.0;
        double gradPhiY = 0.0;
        double gradPhiZ = 0.0;

        const TlPosition pos = gridPosition - this->m_tlOrbInfo.getPosition(AO);
        const double distance2 = pos.squareDistanceFrom();
        const int basisType = this->m_tlOrbInfo.getBasisType(AO);
        const double prefactor = DfCalcGridX::getPrefactor(basisType, pos);
        const int contraction = this->m_tlOrbInfo.getCgtoContraction(AO);

        for (int PGTO = 0; PGTO < contraction; ++PGTO) {
            double prefactorX = 0.0;
            double prefactorY = 0.0;
            double prefactorZ = 0.0;
            const double alpha = this->m_tlOrbInfo.getExponent(AO, PGTO);
            const double shoulder = alpha * distance2;

            if (shoulder <= TOOBIG) {
                const double g = this->m_tlOrbInfo.getCoefficient(AO, PGTO) * std::exp(-shoulder);
                phi += prefactor * g;
                DfCalcGridX::getPrefactorForDerivative(basisType, alpha, pos,
                                                       &prefactorX, &prefactorY, &prefactorZ);
                gradPhiX += prefactorX * g;
                gradPhiY += prefactorY * g;
                gradPhiZ += prefactorZ * g;
            }
        }
        
        if (std::fabs(phi) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, phi);
            const WFGrid wfGrid(AO_index, phi);
            pPhis->push_back(wfGrid);
        }

        if (std::fabs(gradPhiX) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiX);
            const WFGrid wfGrid(AO_index, gradPhiX);
            pGradPhiXs->push_back(wfGrid);
        }
        if (std::fabs(gradPhiY) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiY);
            const WFGrid wfGrid(AO_index, gradPhiY);
            pGradPhiYs->push_back(wfGrid);
        }
        if (std::fabs(gradPhiZ) > densityCutOffValue) {
            //const WFGrid wfGrid(AO, gradPhiZ);
            const WFGrid wfGrid(AO_index, gradPhiZ);
            pGradPhiZs->push_back(wfGrid);
        }
    }

    std::sort(pPhis->begin(), pPhis->end(), WFGrid_sort_functional());
    std::sort(pGradPhiXs->begin(), pGradPhiXs->end(), WFGrid_sort_functional());
    std::sort(pGradPhiYs->begin(), pGradPhiYs->end(), WFGrid_sort_functional());
    std::sort(pGradPhiZs->begin(), pGradPhiZs->end(), WFGrid_sort_functional());
}

// 分散行列専用
void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlMatrix& P,
                                             const std::vector<WFGrid>& rowPhi,
                                             const std::vector<WFGrid>& colPhi,
                                             double* pRhoA)
{
    const double densityCutOffValue = this->m_densityCutOffValueA;
    double dRho = 0.0;
    double dGradRhoX = 0.0;
    double dGradRhoY = 0.0;
    double dGradRhoZ = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    double value = densityCutOffValue;

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(rowPhi.begin(),
                                                                rowPhi.end(),
                                                                WFGrid(0, value),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(rowPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const std::size_t nOrb_p = rowPhi[p].index;
        const double phi_p = rowPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(colPhi.begin(),
                                                                    colPhi.end(),
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        const std::size_t max_q = std::distance(colPhi.begin(), qEnd);
        for (std::size_t q = 0; q < max_q; ++q) {
            const std::size_t nOrb_q = colPhi[q].index;
            const double phi_q = colPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
}


void DfCalcGridX_Parallel::getRhoAtGridPoint(const TlMatrix& P,
                                             const std::vector<WFGrid>& rowPhi,
                                             const std::vector<WFGrid>& colPhi,
                                             const std::vector<WFGrid>& colGradPhiX,
                                             const std::vector<WFGrid>& colGradPhiY,
                                             const std::vector<WFGrid>& colGradPhiZ,
                                             double* pRhoA,
                                             double* pGradRhoAX,
                                             double* pGradRhoAY,
                                             double* pGradRhoAZ)
{
    const double densityCutOffValue = 1.0E-16;
    double dRho = 0.0;
    double dGradRhoX = 0.0;
    double dGradRhoY = 0.0;
    double dGradRhoZ = 0.0;

    // 密度行列のループ(index = p)の最大を求める
    //double value = 1.0E-16;
    // if (aGradPhiX.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiX[0].value));
    // }
    // if (aGradPhiY.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiY[0].value));
    // }
    // if (aGradPhiZ.size() > 0) {
    //     value = std::min(value, std::fabs(aGradPhiZ[0].value));
    // }

    // std::size_t max_q = aPhi.size();
    // std::size_t max_qx = aGradPhiX.size();
    // std::size_t max_qy = aGradPhiY.size();
    // std::size_t max_qz = aGradPhiZ.size();

    std::vector<WFGrid>::const_iterator pEnd = std::upper_bound(rowPhi.begin(),
                                                                rowPhi.end(),
                                                                WFGrid(0, densityCutOffValue),
                                                                WFGrid_sort_functional());
    const std::size_t max_p = std::distance(rowPhi.begin(), pEnd);
    for (std::size_t p = 0; p < max_p; ++p) {
        const index_type nOrb_p = rowPhi[p].index;
        const double phi_p = rowPhi[p].value;
        const double cutValue = std::fabs(densityCutOffValue / phi_p);

        // 高速化
        const TlVector P_row = P.getRowVector(nOrb_p);

        std::vector<WFGrid>::const_iterator qEnd = std::upper_bound(colPhi.begin(),
                                                                    colPhi.end(),
                                                                    WFGrid(0, cutValue),
                                                                    WFGrid_sort_functional());
        const std::size_t max_q = std::distance(colPhi.begin(), qEnd);
        for (std::size_t q = 0; q < max_q; ++q) {
            const std::size_t nOrb_q = colPhi[q].index;
            const double phi_q = colPhi[q].value;
            dRho += P_row[nOrb_q] * phi_p * phi_q;
        }
        
        std::vector<WFGrid>::const_iterator qxEnd = std::upper_bound(colGradPhiX.begin(),
                                                                     colGradPhiX.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qx = std::distance(colGradPhiX.begin(), qxEnd);
        for (std::size_t qx = 0; qx < max_qx; ++qx) {
            const index_type nOrb_q = colGradPhiX[qx].index;
            const double phi_q = colGradPhiX[qx].value;
            dGradRhoX += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qyEnd = std::upper_bound(colGradPhiY.begin(),
                                                                     colGradPhiY.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qy = std::distance(colGradPhiY.begin(), qyEnd);
        for (std::size_t qy = 0; qy < max_qy; ++qy) {
            const std::size_t nOrb_q = colGradPhiY[qy].index;
            const double phi_q = colGradPhiY[qy].value;
            dGradRhoY += P_row[nOrb_q] * phi_p * phi_q;
        }

        std::vector<WFGrid>::const_iterator qzEnd = std::upper_bound(colGradPhiZ.begin(),
                                                                     colGradPhiZ.end(),
                                                                     WFGrid(0, cutValue),
                                                                     WFGrid_sort_functional());
        const std::size_t max_qz = std::distance(colGradPhiZ.begin(), qzEnd);
        for (std::size_t qz = 0; qz < max_qz; ++qz) {
            const std::size_t nOrb_q = colGradPhiZ[qz].index;
            const double phi_q = colGradPhiZ[qz].value;
            dGradRhoZ += P_row[nOrb_q] * phi_p * phi_q;
        }
    }

    *pRhoA = dRho;
    *pGradRhoAX = dGradRhoX;
    *pGradRhoAY = dGradRhoY;
    *pGradRhoAZ = dGradRhoZ;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_LDA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A,
                                      TlDistributeSymmetricMatrix* pF_B)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    TlSparseSymmetricMatrix tmpF_B(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A, &tmpF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    pF_B->mergeSparseMatrix(tmpF_B);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    rComm.allReduce_SUM(energy);
    return energy;
}

double DfCalcGridX_Parallel::buildVxc(DfFunctional_GGA* pFunctional,
                                      TlDistributeSymmetricMatrix* pF_A,
                                      TlDistributeSymmetricMatrix* pF_B)
{
    const TlMatrix gridMat = this->distributeGridMatrix(this->m_nIteration);
    TlSparseSymmetricMatrix tmpF_A(this->m_nNumOfAOs);
    TlSparseSymmetricMatrix tmpF_B(this->m_nNumOfAOs);
    double energy = DfCalcGridX::buildVxc(gridMat, pFunctional, &tmpF_A, &tmpF_B);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    pF_A->mergeSparseMatrix(tmpF_A);
    pF_B->mergeSparseMatrix(tmpF_B);
    rComm.allReduce_SUM(energy);
    return energy;
}

void DfCalcGridX_Parallel::getWholeDensity(double* pRhoA, double* pRhoB) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCalcGridX::getWholeDensity(pRhoA, pRhoB);
    }
}

