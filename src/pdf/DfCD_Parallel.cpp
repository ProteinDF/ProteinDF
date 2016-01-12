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

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "CnError.h"
#include "DfCD_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "DfOverlapEngine.h"
#include "DfEriEngine.h"
#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlSystem.h"
#include "TlUtils.h"

#define TRANS_MEM_SIZE (1 * 1024 * 1024 * 1024) // 1GB
// #define CD_DEBUG

DfCD_Parallel::DfCD_Parallel(TlSerializeData* pPdfParam) 
    : DfCD(pPdfParam) {

    this->isDebugSaveL_ = false;
    if (! (*this->pPdfParam_)["debug/saveL"].getStr().empty()) {
        this->isDebugSaveL_ = (*this->pPdfParam_)["debug/saveL"].getBoolean();
    }
}


DfCD_Parallel::~DfCD_Parallel()
{
}

void DfCD_Parallel::calcCholeskyVectorsForJK()
{
    this->log_.info("calc CholeskyVectors (parallel)");

    // for J & K
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_set"]);
    // productive code for J
    {
        this->createEngines<DfEriEngine>();
        
        const TlRowVectorMatrix Ljk 
            = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                     this->getI2pqVtrPath(),
                                                     this->epsilon_,
                                                     &DfCD::calcDiagonals,
                                                     &DfCD_Parallel::getSuperMatrixElements);
        this->saveL(Ljk, DfObject::getLjkMatrixPath());
        
        this->destroyEngines();
        this->log_.info("");
    }

    // for K
    switch (this->fastCDK_mode_) {
    case FASTCDK_PRODUCTIVE_FULL:
        {
            this->log_.info("fast CDK routine(parallel; full).");
            this->createEngines<DfEriEngine>();
            
            const TlRowVectorMatrix Lk 
                = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                         this->getI2prVtrPath(),
                                                         this->epsilon_K_,
                                                         &DfCD::calcDiagonals_K_full,
                                                         &DfCD_Parallel::getSuperMatrixElements_K_full);
            this->saveL(Lk, DfObject::getLkMatrixPath());
            
            this->destroyEngines();
            this->log_.info("");
        }
        break;

    case FASTCDK_PRODUCTIVE:
        {
            this->log_.info("fast CDK routine(parallel).");
            this->createEngines<DfEriEngine>();
            
            const TlRowVectorMatrix Lk 
                = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                         this->getI2prVtrPath(),
                                                         this->epsilon_K_,
                                                         &DfCD::calcDiagonals_K_half,
                                                         &DfCD_Parallel::getSuperMatrixElements_K_half);
            this->saveL(Lk, DfObject::getLkMatrixPath());
            
            this->destroyEngines();
            this->log_.info("");
        }
        break;

    case FASTCDK_NONE:
        this->log_.info("fast CDK routine is invalided.");
        break;

    default:
        this->log_.critical("program error");
        CnErr.abort();
        break;
    }
}


void DfCD_Parallel::calcCholeskyVectorsForGridFree()
{
    // for XC(gridfree)
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);

    if (this->isDedicatedBasisForGridFree_) {
        const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                      (*this->pPdfParam_)["basis_set_gridfree"]);
        // productive code
        const TlRowVectorMatrix Lxc 
            = DfCD::calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                 orbInfo_q,
                                                                 this->getI2pqVtrXCPath());
        this->saveL(Lxc, DfObject::getLxcMatrixPath());
    } else {
        // productive code
        this->log_.info("build Lxc matrix by on-the-fly method.");
        this->createEngines<DfOverlapEngine>();

        const TlRowVectorMatrix Lxc
            = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo_p,
                                                     this->getI2pqVtrXCPath(),
                                                     this->epsilon_,
                                                     &DfCD::calcDiagonals,
                                                     &DfCD_Parallel::getSuperMatrixElements);

        this->saveL(Lxc, DfObject::getLxcMatrixPath());
        this->destroyEngines();
    }
}

DfTaskCtrl* DfCD_Parallel::getDfTaskCtrlObject() const
{
    DfTaskCtrl *pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}

void DfCD_Parallel::finalize(TlSymmetricMatrix *pMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMat);
}

void DfCD_Parallel::finalize(TlSparseSymmetricMatrix *pMat)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    // rComm.gatherToMaster(*pMat);
    // rComm.broadcast(*pMat);
    rComm.allReduce_SUM(*pMat);
}


void DfCD_Parallel::finalize_I2PQ(PQ_PairArray* pI2PQ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    
    // gather to master
    if (rComm.isMaster() == true) {
        const int slaves = numOfProcs - 1;
        int proc = 0;
        std::vector<index_type> shellArray;
        for (int i = 0; i < slaves; ++i) {
            rComm.receiveDataFromAnySource(shellArray, &proc);

            const std::size_t I2PQ_size = shellArray.size() / 2;
            PQ_PairArray i2pq_tmp(I2PQ_size);
            for (std::size_t i = 0; i < I2PQ_size; ++i) {
                i2pq_tmp[i] = Index2(shellArray[i*2   ],
                                     shellArray[i*2 +1]);
            }
            
            pI2PQ->insert(pI2PQ->end(),
                          i2pq_tmp.begin(), i2pq_tmp.end());
        }
        //std::sort(pI2PQ->begin(), pI2PQ->end(), PQ_Pair_less());
        std::sort(pI2PQ->begin(), pI2PQ->end());
    } else {
        const std::size_t I2PQ_size = pI2PQ->size();
        std::vector<index_type> shellArray(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i*2   ] = (*pI2PQ)[i].index1();
            shellArray[i*2 +1] = (*pI2PQ)[i].index2();
        }
        rComm.sendData(shellArray);
    }
    
    // broadcast
    {
        std::vector<index_type> tmp;
        if (rComm.isMaster() == true) {
            const std::size_t size = pI2PQ->size();
            tmp.resize(size * 2);
            for (std::size_t i = 0; i < size; ++i) {
                tmp[i*2   ] = (*pI2PQ)[i].index1();
                tmp[i*2 +1] = (*pI2PQ)[i].index2();
            }
        }
        rComm.broadcast(tmp);
        if (rComm.isMaster() != true) {
            const std::size_t size = tmp.size() / 2;
            pI2PQ->resize(size);
            for (std::size_t i = 0; i < size; ++i) {
                (*pI2PQ)[i] = Index2(tmp[i*2   ],
                                     tmp[i*2 +1]);
            }
        }
    }
}


void DfCD_Parallel::saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCD::saveI2PQ(I2PQ, filepath);
    }
}


DfCD::PQ_PairArray DfCD_Parallel::getI2PQ(const std::string& filepath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("distribute I2PQ table.");

    PQ_PairArray I2PQ;
    std::vector<index_type> shellArray;
    if (rComm.isMaster() == true) {
        I2PQ = DfCD::getI2PQ(filepath);

        const std::size_t I2PQ_size = I2PQ.size(); 
        shellArray.resize(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i*2   ] = I2PQ[i].index1();
            shellArray[i*2 +1] = I2PQ[i].index2();
        }
    }
    
    rComm.broadcast(shellArray);
    if (rComm.isMaster() != true) {
        const std::size_t I2PQ_size = shellArray.size() / 2;
        I2PQ.resize(I2PQ_size);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            I2PQ[i] = Index2(shellArray[i*2   ],
                             shellArray[i*2 +1]);
        }
    }

    return I2PQ;
}


TlSymmetricMatrix DfCD_Parallel::getPMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlMatrix P;
    if (rComm.isMaster() == true) {
        P = DfCD::getPMatrix();
    }
    rComm.broadcast(P);
    return P;
}


void DfCD_Parallel::divideCholeskyBasis(const index_type numOfCVs,
                                        index_type *pStart, index_type *pEnd)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();
    
    const index_type range = (numOfCVs + numOfProcs -1) / numOfProcs;
    *pStart = range * rank;
    *pEnd = std::min(range * (rank +1), numOfCVs);
}


TlDistributeSymmetricMatrix 
DfCD_Parallel::getCholeskyVector_distribute(const TlVector& L_col,
                                            const PQ_PairArray& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlDistributeSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

TlDistributeSymmetricMatrix 
DfCD_Parallel::getCholeskyVectorA_distribute(const TlOrbitalInfoObject& orbInfo_p,
                                             const TlOrbitalInfoObject& orbInfo_q,
                                             const TlVector& L_col,
                                             const PQ_PairArray& I2PQ)
{
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();

    const index_type numOfItilde = L_col.getSize();
    TlDistributeMatrix answer(numOfOrbs_p, numOfOrbs_q);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}


TlRowVectorMatrix 
DfCD_Parallel::calcCholeskyVectorsOnTheFlyS_new(const TlOrbitalInfoObject& orbInfo,
                                                const std::string& I2PQ_path,
                                                const double threshold,
                                                void (DfCD_Parallel::*calcDiagonalsFunc)(
                                                    const TlOrbitalInfoObject&,
                                                    PQ_PairArray*,
                                                    TlVector*),
                                                void (DfCD_Parallel::*getSuperMatrixElementsFunc)(
                                                    const TlOrbitalInfoObject&,
                                                    const index_type,
                                                    const std::vector<index_type>&,
                                                    const PQ_PairArray&,
                                                    std::vector<double>*))
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int myRank = rComm.getRank();

    this->log_.info("call on-the-fly Cholesky Decomposition routine (parallel; symmetric)");
    assert(this->pEngines_ != NULL);

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs +1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));

    // CDAM
    assert(this->pEngines_ != NULL);
    PQ_PairArray I2PQ;
    TlVector diagonals; // 対角成分
    (this->*calcDiagonalsFunc)(orbInfo, &I2PQ, &diagonals);
    assert((std::size_t)diagonals.getSize() == I2PQ.size());

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);

    // debug
    // this->debug_I2PQ_ = I2PQ;

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", threshold));
    const index_type numOfPQtilde = I2PQ.size();
    TlRowVectorMatrix L(numOfPQtilde, 1,
                        rComm.getNumOfProcs(), myRank,
                        this->isEnableMmap_);

    double error = diagonals.getMaxAbsoluteElement();
    std::vector<TlVector::size_type> pivot(numOfPQtilde);
    for (index_type i = 0; i < numOfPQtilde; ++i) {
        pivot[i] = i;
    }

    int progress = 0;
    index_type division =  std::max<index_type>(numOfPQtilde * 0.01, 100);
    L.reserveColSize(division);

    index_type numOfCDVcts = 0;
    while ((error > threshold) && (numOfCDVcts < numOfPQtilde)) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", numOfCDVcts, numOfPQtilde, error));
#endif //DEBUG_CD

        // progress 
        if (numOfCDVcts >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e, local mem:%8.1f MB",
                                            numOfCDVcts, error, TlSystem::getMaxRSS()));
            ++progress;

            // メモリの確保
            L.reserveColSize(division * progress);
        }
        L.resize(numOfPQtilde, numOfCDVcts +1);

        // pivot
        {
            std::vector<TlVector::size_type>::const_iterator it = diagonals.argmax(pivot.begin() + numOfCDVcts,
                                                                                   pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[numOfCDVcts], pivot[i]);
        }

        const index_type pivot_m = pivot[numOfCDVcts];
        error = diagonals[pivot_m];
        if (error < threshold) {
            break;
        }

        const double l_m_pm = std::sqrt(diagonals[pivot_m]);
        const double inv_l_m_pm = 1.0 / l_m_pm;
        L.set(pivot_m, numOfCDVcts, l_m_pm);

        // get supermatrix elements
        const index_type numOf_G_cols = numOfPQtilde -(numOfCDVcts +1);
        std::vector<double> G_pm(numOf_G_cols);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[(numOfCDVcts +1) +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            (this->*getSuperMatrixElementsFunc)(orbInfo,
                                                pivot_m, G_col_list, I2PQ, &G_pm);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);
        std::vector<double> output_G_pm(numOf_G_cols, 0.0);
        rComm.iAllReduce_SUM(&(G_pm[0]), &(output_G_pm[0]), numOf_G_cols);

        // CD calc
        TlVector L_pm(numOfCDVcts +1);
        {
            const int PE_in_charge = L.getSubunitID(pivot_m);
            if (PE_in_charge == myRank) {
                L_pm = L.getVector(pivot_m);
                assert(L_pm.getSize() == numOfCDVcts +1);
            }
            rComm.broadcast(L_pm, PE_in_charge);
        }
        assert(L_pm.getSize() == numOfCDVcts +1);

        std::vector<double> L_xm(numOf_G_cols);
        TlVector update_diagonals(numOfPQtilde);
        rComm.wait(&(output_G_pm[0]));
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[(numOfCDVcts +1) +i]; // from (m+1) to N
            if (L.getSubunitID(pivot_i) == myRank) {
                TlVector L_pi = L.getVector(pivot_i);
                const double sum_ll = (L_pi.dot(L_pm)).sum();
                // const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;
                const double l_m_pi = (output_G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp atomic
                L_xm[i] += l_m_pi;
                
#pragma omp atomic
                update_diagonals[pivot_i] -= l_m_pi * l_m_pi;
            }
        }

        // allReduce_SUM (L_xm)
        std::vector<double> output_L_xm(numOf_G_cols, 0.0);
        rComm.iAllReduce_SUM(&(L_xm[0]), &(output_L_xm[0]), numOf_G_cols);

        // allReduce_SUM (update_diagonals)
        std::vector<double> output_update_diagonals(numOfPQtilde, 0.0);
        rComm.iAllReduce_SUM(&(update_diagonals[0]), &(output_update_diagonals[0]), numOfPQtilde);


        rComm.wait(&(output_L_xm[0]));
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[(numOfCDVcts +1) +i]; // from (m+1) to N
            L.set(pivot_i, numOfCDVcts, output_L_xm[i]);
        }

        rComm.wait(&(output_update_diagonals[0]));
        diagonals += output_update_diagonals;


        error = diagonals[pivot[numOfCDVcts]];
        ++numOfCDVcts;
    }
    L.resize(numOfPQtilde, numOfCDVcts);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", numOfCDVcts));

    return L;
}


TlRowVectorMatrix DfCD_Parallel::calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
                                                              const TlOrbitalInfoObject& orbInfo_q,
                                                              const std::string& I2PQ_path)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("call on-the-fly Cholesky Decomposition routine (MPI parallel; asymmetric)");
    assert(this->pEngines_ != NULL);
    this->initializeCutoffStats(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfOrbs_p * numOfOrbs_q;
    this->log_.info(TlUtils::format("number of orbitals1: %d", numOfOrbs_p));
    this->log_.info(TlUtils::format("number of orbitals2: %d", numOfOrbs_q));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));
    PQ_PairArray I2PQ;
    TlSparseMatrix schwartzTable(numOfOrbs_p, numOfOrbs_q);
    TlVector global_diagonals; // 対角成分

    assert(this->pEngines_ != NULL);
    this->calcDiagonalsA(orbInfo_p, orbInfo_q,
                         &I2PQ, &schwartzTable, &global_diagonals);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();

    // prepare variables (parallel)
    TlRowVectorMatrix L(N, 1,
                        rComm.getNumOfProcs(),
                        rComm.getRank(),
                        this->isEnableMmap_); // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
    const index_type local_N = L.getNumOfLocalVectors();
    TlVector L_pm(N);
    std::vector<int> global_pivot(N); // ERI計算リストを作成するために必要
    std::vector<int> reverse_pivot(N); // global_pivotの逆引き
    std::vector<int> local_pivot(local_N);
    std::vector<double> local_diagonals(local_N);
    const int myRank = rComm.getRank();

    double error = 0.0;
    index_type error_global_loc = 0;
    index_type error_local_loc = 0;
    {
        int local_i = 0;
        for (int global_i = 0; global_i < N; ++global_i) {
            global_pivot[global_i] = global_i;
            reverse_pivot[global_i] = global_i;
            if (L.getSubunitID(global_i) == myRank) {
                local_pivot[local_i] = global_i;
                local_diagonals[local_i] = global_diagonals[global_i];
                if (error < local_diagonals[local_i]) {
                    error_local_loc = local_i;
                    error_global_loc = local_pivot[local_i];
                    error = local_diagonals[local_i];
                }
                ++local_i;
            }
        }
        assert(local_i == local_N);
    }
    rComm.allReduce_MAXLOC(&error, &error_global_loc);

    index_type m = 0;
    index_type local_m = 0;
    if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
        std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
        ++local_m;
    }

    int progress = 0;
    index_type division = std::max<index_type>(N * 0.01, 100);
    while (error > threshold) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif // DEBUG_CD

        // progress 
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e",
                                            m, error));
            ++progress;

            // メモリの確保
            L.reserveColSize(progress * division);
        }
        L.resize(N, m+1);

        // pivot
        std::swap(global_pivot[m], global_pivot[error_global_loc]);
        reverse_pivot[global_pivot[m]] = m;
        reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

        const double l_m_pm = std::sqrt(error);
        const index_type pivot_m = global_pivot[m];
        L.set(pivot_m, m, l_m_pm); // 通信発生せず。関係無いPEは値を捨てる。
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // calc
        std::vector<double> G_pm;
        // const index_type numOf_G_cols = N -(m+1);
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = global_pivot[m+1 +i]; // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElementsA(orbInfo_p, orbInfo_q,
                                                 pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        {
            // 全PEに分配
            const int PEinCharge = L.getSubunitID(pivot_m);
            if (PEinCharge == rComm.getRank()) {
                //const index_type copySize = L.getVector(pivot_m, &(L_pm[0]), m +1);
                L_pm = L.getVector(pivot_m);
                //assert(L_pm[0].size() == m +1);
            }
            rComm.broadcast(&(L_pm[0]), m +1, PEinCharge);
        }

        error = 0.0;
#pragma omp parallel
        {
            TlVector L_pi(m +1);
            double my_error = 0.0;
            int my_error_global_loc = 0;
            int my_error_local_loc = 0;

#pragma omp for schedule(runtime)
            for (int i = local_m; i < local_N; ++i) {
                const int pivot_i = local_pivot[i];
                //const index_type copySize = L.getVector(pivot_i, &(L_pi[0]), m +1);
                L_pi = L.getVector(pivot_i);
                // assert(L_pi == m +1);
                double sum_ll = 0.0;
                for (index_type j = 0; j < m; ++j) {
                    sum_ll += L_pm[j] * L_pi[j];
                }

                const int G_pm_index = reverse_pivot[pivot_i] - (m+1);
                const double l_m_pi = (G_pm[G_pm_index] - sum_ll) * inv_l_m_pm;
                const double ll = l_m_pi * l_m_pi;

#pragma omp critical(DfCD_Parallel__calcCholeskyVectorsOnTheFlyA_updateL)
                {
                    L.set(pivot_i, m, l_m_pi);
                    global_diagonals[pivot_i] -= ll;
                }

                if (global_diagonals[pivot_i] > my_error) {
                    my_error = global_diagonals[pivot_i];
                    my_error_global_loc = reverse_pivot[pivot_i]; // == m +1 + i
                    my_error_local_loc = i;
                }
            }

#pragma omp critical(DfCD_Parallel__calcCholeskyVectorsOnTheFlyA_update_error)
            {
                if (error < my_error) {
                    error = my_error;
                    error_global_loc = my_error_global_loc;
                    error_local_loc = my_error_local_loc;
                }
            }
        }
        rComm.allReduce_MAXLOC(&error, &error_global_loc);
        global_diagonals[global_pivot[error_global_loc]] = error;

        ++m;
        if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
            std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
            ++local_m;
        }
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    return L;
}


void DfCD_Parallel::getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                                           const index_type G_row,
                                           const std::vector<index_type>& G_col_list,
                                           const PQ_PairArray& I2PQ,
                                           std::vector<double>* pSuperMatrixElements)
{
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs -1) / numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank +1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = this->getCalcList(orbInfo, G_row,
                                                                G_col_list, start, end, I2PQ);
    this->calcERIs(orbInfo, calcList);
    *pSuperMatrixElements = this->setERIs(orbInfo, G_row,
                                          G_col_list, start, end, I2PQ);

    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}


void DfCD_Parallel::getSuperMatrixElements_K_full(const TlOrbitalInfoObject& orbInfo,
                                                  const index_type G_row,
                                                  const std::vector<index_type>& G_col_list,
                                                  const PQ_PairArray& I2PR,
                                                  std::vector<double>* pSuperMatrixElements)
{
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs -1) / numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank +1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = DfCD::getCalcList_K_full(orbInfo, G_row,
                                                                       G_col_list, start, end, I2PR);
    DfCD::calcERIs_K(orbInfo, calcList);
    *pSuperMatrixElements = DfCD::setERIs_K_full(orbInfo, G_row,
                                                 G_col_list, start, end, I2PR);

    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}


void DfCD_Parallel::getSuperMatrixElements_K_half(const TlOrbitalInfoObject& orbInfo,
                                                  const index_type G_row,
                                                  const std::vector<index_type>& G_col_list,
                                                  const PQ_PairArray& I2PR,
                                                  std::vector<double>* pSuperMatrixElements)
{
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs -1)/ numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank +1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = DfCD::getCalcList_K_half(orbInfo, G_row,
                                                                       G_col_list, start, end, I2PR);
    DfCD::calcERIs_K(orbInfo, calcList);
    *pSuperMatrixElements = DfCD::setERIs_K_half(orbInfo, G_row,
                                                 G_col_list, start, end, I2PR);
    
    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}


void DfCD_Parallel::saveL(const TlRowVectorMatrix& L,
                          const std::string& path)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    this->log_.info("transform the Cholesky vectors.");
    TlColVectorMatrix colVecL = this->transLMatrix(L);

    this->log_.info("save the Cholesky vectors.");
    colVecL.save(path);
}


TlMatrix DfCD_Parallel::mergeL(const TlRowVectorMatrix& L)
{
    this->log_.info("merge L(row): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfRows, numOfProcs);
    const index_type localRows = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type r = 0; r < localRows; ++r) {
        const index_type row = r * numOfProcs + rank;
        TlVector rowVec = L.getVector(row);
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, rowVec[col]);
        }
    }
    rComm.allReduce_SUM(answer);
    
    this->log_.info("merge L: end");
    return answer;
}


TlMatrix DfCD_Parallel::mergeL(const TlColVectorMatrix& L)
{
    this->log_.info("merge L(col): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfCols, numOfProcs);
    const index_type localCols = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type c = 0; c < localCols; ++c) {
        const index_type col = c * numOfProcs + rank;
        TlVector colVec = L.getVector(col);
        for (index_type row = 0; row < numOfRows; ++row) {
            answer.set(row, col, colVec[row]);
        }
    }
    rComm.allReduce_SUM(answer);
    
    this->log_.info("merge L: end");
    return answer;
}


void DfCD_Parallel::getJ(TlSymmetricMatrix* pJ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc J by CD method (parallel).");

    // cholesky vector
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                        isUsingMemManager);
    // TlColVectorMatrix L;
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlVector vP = this->getScreenedDensityMatrix(I2PQ);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PQ.size());
    const index_type numOfCVs = L.getNumOfCols();

    TlVector vJ(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            const TlVector LI = L.getVector(I);
            
            TlVector tmpLI = LI;
            const double qi = tmpLI.dot(vP).sum();

            vJ += qi*LI;
        }
    }

    rComm.allReduce_SUM(vJ);
    this->expandJMatrix(vJ, I2PQ, pJ);
}


TlSymmetricMatrix DfCD_Parallel::getPMatrix(const RUN_TYPE runType, const int itr)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlSymmetricMatrix P;
    if (rComm.isMaster()) {
        P = DfCD::getPMatrix(runType, itr);
    }
    rComm.broadcast(P);
    return P;
}


void DfCD_Parallel::getK_S_woCD(const RUN_TYPE runType,
                                TlSymmetricMatrix* pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K by CD method (parallel).");

    // cholesky vector
    // TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    TlColVectorMatrix L;
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    
    TlSymmetricMatrix P = this->getPMatrix(runType, this->m_nIteration -1);
    
    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            TlSymmetricMatrix l = 
                this->getCholeskyVector(cv, I2PQ);
            
            TlMatrix X = l * P;
            l.transpose();
            X *= l;

            *pK += X;
        }
    }
    rComm.allReduce_SUM(*pK);
    
    *pK *= -1.0;
}


void DfCD_Parallel::getK_S_fast(const RUN_TYPE runType,
                                TlSymmetricMatrix* pK)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K(fast) by CD method (parallel).");

    // cholesky vector
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                        isUsingMemManager);
    // TlColVectorMatrix L;
    L.load(DfObject::getLkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PR = this->getI2PQ(this->getI2prVtrPath());
    const TlVector vP = this->getScreenedDensityMatrix(runType, I2PR);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PR.size());
    const index_type numOfCVs = L.getNumOfCols();

    TlVector vK(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            const TlVector LI = L.getVector(I);
            
            TlVector tmpLI = LI;
            const double qi = tmpLI.dot(vP).sum();

            vK += qi*LI;
        }
    }
    vK *= -1.0;

    rComm.allReduce_SUM(vK);
    this->expandKMatrix(vK, I2PR, pK);
}


void DfCD_Parallel::getJ_D(TlDistributeSymmetricMatrix* pJ)
{
    this->log_.info("calc J by CD method (parallel; distributed).");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // cholesky vector
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(), 
                        isUsingMemManager);
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlVector vP = this->getScreenedDensityMatrixD(I2PQ);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PQ.size());
    const index_type numOfCVs = L.getNumOfCols();

    TlVector vJ(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            const TlVector LI = L.getVector(I);
            
            TlVector tmpLI = LI;
            const double qi = tmpLI.dot(vP).sum();

            vJ += qi*LI;
        }
    }

    rComm.allReduce_SUM(vJ);
    this->expandJMatrixD(vJ, I2PQ, pJ);
}

TlVector DfCD_Parallel::getScreenedDensityMatrixD(const PQ_PairArray& I2PQ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlDistributeSymmetricMatrix P = DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    const std::size_t numOfI = I2PQ.size();
    TlVector answer(numOfI);
    
    for (std::size_t i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        const double coef = (r != c) ? 2.0 : 1.0;
        answer.set(i, coef * P.getLocal(r, c));
    }
    rComm.allReduce_SUM(answer);

    return answer;
}

void DfCD_Parallel::expandJMatrixD(const TlVector& vJ, const PQ_PairArray& I2PQ, TlDistributeSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    const index_type numOfI = I2PQ.size();
    for (index_type i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        pJ->set(pair.index1(), pair.index2(), vJ.get(i));
    }
}

void DfCD_Parallel::getK_D(const RUN_TYPE runType,
                           TlDistributeSymmetricMatrix* pK)
{
    this->log_.info("calc K by CD method (parallel; distributed).");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // cholesky vector
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                         isUsingMemManager);
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();

    this->log_.info("calc CD of density matrix");
    this->log_.info(TlUtils::format("epsilon = %8.3e", this->epsilon_));
    TlDistributeSymmetricMatrix P = 
        0.5 * DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1); // RKS
    const TlDistributeMatrix C = P.choleskyFactorization(this->epsilon_);
    //const TlDistributeMatrix C = P.choleskyFactorization_mod2(this->epsilon_);

    TlTime time_all;
    TlTime time_bcast;
    TlTime time_translate;
    TlTime time_multimat1;
    TlTime time_multimat2;
    TlTime time_sym2gen;
    TlTime time_transpose;
    TlTime time_add;
    
    time_all.start();
    this->log_.info("start loop");
    int progress = 0;
    const int division = numOfCVs * 0.1;
    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        // progress
        if (I >= progress * division) {
            const double rate = double(I) / double(numOfCVs) * 100.0;
            this->log_.info(TlUtils::format("K loop progress: %5.2f%%", rate));
            ++progress;
        }

        time_bcast.start();
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            //const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);
        time_bcast.stop();

        time_translate.start();
        TlDistributeSymmetricMatrix l = 
            this->getCholeskyVector_distribute(cv, I2PQ);
        time_translate.stop();

        time_multimat1.start();
        TlDistributeMatrix X = l * C;
        time_multimat1.stop();

        time_sym2gen.start();
        TlDistributeMatrix Xt = X;
        time_sym2gen.stop();

        time_transpose.start();
        Xt.transpose();
        time_transpose.stop();
        
        time_multimat2.start();
        TlDistributeSymmetricMatrix XX = X * Xt;
        time_multimat2.stop();

        time_add.start();
        *pK += XX;
        time_add.stop();
    }
    
    *pK *= -1.0;
    time_all.stop();

    // timing data
    this->log_.info(TlUtils::format("K all:       %10.1f sec.", time_all.getElapseTime()));
    this->log_.info(TlUtils::format("K bcast:     %10.1f sec.", time_bcast.getElapseTime()));
    this->log_.info(TlUtils::format("K translate: %10.1f sec.", time_translate.getElapseTime()));
    this->log_.info(TlUtils::format("K multi1:    %10.1f sec.", time_multimat1.getElapseTime()));
    this->log_.info(TlUtils::format("K multi2:    %10.1f sec.", time_multimat2.getElapseTime()));
    this->log_.info(TlUtils::format("K sym2gen:   %10.1f sec.", time_sym2gen.getElapseTime()));
    this->log_.info(TlUtils::format("K transpose: %10.1f sec.", time_transpose.getElapseTime()));
    this->log_.info(TlUtils::format("K add:       %10.1f sec.", time_add.getElapseTime()));
}

void DfCD_Parallel::getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    this->log_.info("DfCD_Parallel::getM()");
    if (this->isDedicatedBasisForGridFree_) {
        this->getM_A(P, pM);
    } else {
        this->getM_S(P, pM);
    }
}

void DfCD_Parallel::getM_S(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (symmetric routine; parallel)");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_set"]);
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pM->resize(numOfAOs);

    // cholesky vector
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            //const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            TlSymmetricMatrix LI = 
                this->getCholeskyVector(cv, I2PQ);
            assert(LI.getNumOfRows() == this->m_nNumOfAOs);
            assert(LI.getNumOfCols() == this->m_nNumOfAOs);

            TlSymmetricMatrix QI = LI;
            QI.dot(P);
            const double qi = QI.sum();
            
            *pM += qi*LI;
        }
    }

    rComm.allReduce_SUM(*pM);
}

void DfCD_Parallel::getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (asymmetric routine; parallel)");

    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set_gridfree"]);
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            //const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            TlMatrix l = 
                this->getCholeskyVectorA(orbInfo_p, orbInfo_q,
                                         cv, I2PQ);
            l.transpose();
            
            TlMatrix X = l * C;
            TlMatrix Xt = X;
            Xt.transpose();
            
            TlSymmetricMatrix XX = X * Xt;
            *pM += XX;
        }
    }
    rComm.allReduce_SUM(*pM);
}

void DfCD_Parallel::getM(const TlDistributeSymmetricMatrix& P,
                         TlDistributeSymmetricMatrix* pM)
{
    this->log_.critical("sorry, NO IMPLEMENTED!");
    this->log_.critical("DfCD_Parallel::getM(const TlDistributeSymmetricMatrix&, TlDistributeSymmetricMatrix*)");
    abort();
    // if (this->isDedicatedBasisForGridFree_) {
    //     this->getM_A(P, pM);
    // } else {
    //     this->getM_S(P, pM);
    // }
}

void DfCD_Parallel::getM_S(const TlDistributeSymmetricMatrix& P,
                           TlDistributeSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (symmetric routine; parallel; distributed)");

    // const TlDistributeSymmetricMatrix P = 
    //     DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

    // cholesky vector
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(), 
                         this->isEnableMmap_);
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);

        TlDistributeSymmetricMatrix LI = 
            this->getCholeskyVector_distribute(cv, I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlDistributeSymmetricMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();
        
        *pM += qi*LI;
    }
}

void DfCD_Parallel::getM_A(const TlDistributeSymmetricMatrix& P,
                           TlDistributeSymmetricMatrix* pM)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (asymmetric routine; parallel; distributed)");

    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set_gridfree"]);
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlColVectorMatrix L(1, 1, rComm.getNumOfProcs(), rComm.getRank(),
                         this->isEnableMmap_);
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    this->log_.info("calc CD of density matrix");
    this->log_.info(TlUtils::format("epsilon = %8.3e", this->epsilon_));
    // TlDistributeSymmetricMatrix P = 
    //     0.5 * DfObject::getPpqMatrix<TlDistributeSymmetricMatrix>(RUN_RKS, this->m_nIteration -1); // RKS
    const TlDistributeMatrix C = P.choleskyFactorization_mod2(this->epsilon_);

    this->log_.info("start loop");
    int progress = 0;
    const int division = numOfCVs * 0.1;
    TlVector cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        // progress
        if (I >= progress * division) {
            const double rate = double(I) / double(numOfCVs) * 100.0;
            this->log_.info(TlUtils::format("K loop progress: %5.2f%%", rate));
            ++progress;
        }

        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);

        TlDistributeMatrix l = 
            this->getCholeskyVectorA_distribute(orbInfo_p, orbInfo_q,
                                                cv, I2PQ);
        l.transpose();

        TlDistributeMatrix X = l * C;
        TlDistributeMatrix Xt = X;
        Xt.transpose();
        TlDistributeSymmetricMatrix XX = X * Xt;
        *pM += XX;
    }
}


TlColVectorMatrix DfCD_Parallel::transLMatrix(const TlRowVectorMatrix& rowVectorMatrix)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();
    const int tag = 12345; // 適当

    const index_type numOfRows = rowVectorMatrix.getNumOfRows();
    const index_type numOfCols = rowVectorMatrix.getNumOfCols();
    
    const bool isUsingMemManager = this->isEnableMmap_;
    TlColVectorMatrix colVectorMatrix(numOfRows, numOfCols,
                                      numOfProcs, myRank, isUsingMemManager);

    // prepare partial matrix
    std::vector<unsigned long> matrixElementsSizes(numOfProcs);
    std::vector<std::vector<TlMatrixElement> > matrixElements(numOfProcs);
    {
        std::vector<TlSparseMatrix> partMat(numOfProcs, TlSparseMatrix(numOfRows, numOfCols));
        for (index_type row = 0; row < numOfRows; ++row) {
            if (rowVectorMatrix.getSubunitID(row) == myRank) {
                for (index_type col = 0; col < numOfCols; ++col) {
                    const double v = rowVectorMatrix.get(row, col);
                    
                    const int unit = colVectorMatrix.getSubunitID(col);
                    partMat[unit].set(row, col, v); // TODO: ベクトルのアクセスにするべき
                }
            }
        }

        for (int proc = 0; proc < numOfProcs; ++proc) {
            matrixElementsSizes[proc] = partMat[proc].getSize();
            matrixElements[proc] = partMat[proc].getMatrixElements();
        }
    }

    // send 
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc != myRank) {
            rComm.iSendData(matrixElementsSizes[proc], proc, tag);

            if (matrixElementsSizes[proc] > 0) {
                rComm.iSendDataX(&(matrixElements[proc][0]),
                                 matrixElementsSizes[proc],
                                 proc, tag +1);
            }
        }
    }

    // my partMat
    {
        std::vector<TlMatrixElement>::const_iterator itEnd = matrixElements[myRank].end();
        for (std::vector<TlMatrixElement>::const_iterator it = matrixElements[myRank].begin(); it != itEnd; ++it) {
            const index_type r = it->row;
            const index_type c = it->col;
            colVectorMatrix.set(r, c, it->value);
        }
    }

    // recv
    {
        for (int i = 0; i < (numOfProcs -1); ++i) {
            int proc = 0;
            unsigned long tmp_size = 0;
            rComm.receiveDataFromAnySource(tmp_size, &proc, tag);

            if (tmp_size > 0) {
                std::vector<TlMatrixElement> tmp_elements(tmp_size);
                rComm.receiveDataX(&(tmp_elements[0]), tmp_size, proc, tag +1);

                std::vector<TlMatrixElement>::const_iterator itEnd = tmp_elements.end();
                for (std::vector<TlMatrixElement>::const_iterator it = tmp_elements.begin(); it != itEnd; ++it) {
                    colVectorMatrix.set(it->row, it->col, it->value);
                }
            }
        }
    }

    // wait
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc != myRank) {
            rComm.wait(&(matrixElementsSizes[proc]));
            rComm.wait(&(matrixElements[proc][0]));
        }
    }

    return colVectorMatrix;
}


