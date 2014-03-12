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

#include <set>
#include "DfCD.h"
#include "DfEngineObject.h"
#include "DfEriEngine.h"
#include "DfOverlapEngine.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "TlSystem.h"

DfCD::DfCD(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEngines_(NULL)
{
    // this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    
    this->cutoffEpsilon3_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon3"].getStr().empty() != true) {
        this->cutoffEpsilon3_ = (*pPdfParam)["cutoff_epsilon3"].getDouble();
    }    

    this->CDAM_tau_ = 1.0E-10;
    if ((*pPdfParam)["CDAM_tau"].getStr().empty() != true) {
        this->CDAM_tau_ = (*pPdfParam)["CDAM_tau"].getDouble();
    }    

    this->epsilon_ = 1.0E-4;
    if ((*pPdfParam)["CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["CD_epsilon"].getDouble();
    }    

    this->isStoreERIs_ = true;
    if ((*pPdfParam)["CD_store_ERIs"].getStr().empty() != true) {
        this->isStoreERIs_ = (*pPdfParam)["CD_store_ERIs"].getBoolean();
    }    

    this->debugBuildSuperMatrix_ = false;
    if ((*pPdfParam)["debug/DfGridFreeXC/build_supermatrix"].getStr().empty() != true) {
        this->isStoreERIs_ = (*pPdfParam)["debug/DfGridFreeXC/build_supermatrix"].getBoolean();
    }    
}

DfCD::~DfCD()
{
}

void DfCD::destroyEngines()
{
    this->log_.info("delete engine");
    const int numOfThreads = this->numOfThreads_;
    if (this->pEngines_ != NULL) {
        for (int i = 0; i < numOfThreads; ++i) {
            delete this->pEngines_[i];
            this->pEngines_[i] = NULL;
        }
        delete[] this->pEngines_;
    }
    this->pEngines_ = NULL;
}

void DfCD::calcCholeskyVectorsForJK()
{
    // for J & K
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets"]);
    if (this->debugBuildSuperMatrix_) {
        // DEBUG code
        this->log_.info("call DEBUG routine:");
        this->log_.info("build L matrix by supermatrix.");
        this->createEngines<DfEriEngine>();
        TlSymmetricMatrix V = this->getSuperMatrix(orbInfo);
        this->destroyEngines();
        V.save("fl_Work/debug_Vjk.mat");
            
        TlMatrix L = this->calcCholeskyVectors(V);
        this->saveLjk(L);
    } else {
        // productive code
        const TlRowVectorMatrix2 Ljk 
            = this->calcCholeskyVectorsOnTheFly<DfEriEngine>(orbInfo,
                                                             this->getI2pqVtrPath());
        this->saveLjk(Ljk.getTlMatrix());
    }

    // check
    // {
    //     this->log_.info("check: LL = L * L^t");
    //     TlMatrix L = this->getLxc();
    //     TlMatrix tL = L;
    //     tL.transpose();
    //     TlMatrix LL = L * tL;
    //     LL.save("fl_Work/debug_LL.mat");
    // }
}

void DfCD::calcCholeskyVectorsForGridFree()
{
    // for XC(gridfree)
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);

    if (this->isDedicatedBasisForGridFree_) {
        const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                      (*this->pPdfParam_)["basis_sets_GF"]);
        if (this->debugBuildSuperMatrix_) {
            // DEBUG code
            this->log_.info("call DEBUG routine:");
            this->log_.info("build Lxc matrix by supermatrix.");
            this->createEngines<DfOverlapEngine>();
            TlSymmetricMatrix V = this->getSuperMatrix(orbInfo_p, orbInfo_q);
            this->destroyEngines();
            V.save("fl_Work/debug_Vxc.mat");
            
            TlMatrix L = this->calcCholeskyVectors(V);
            this->saveLxc(L);
        } else {
            // productive code
            const TlRowVectorMatrix2 Lxc 
                = this->calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                     orbInfo_q,
                                                                     this->getI2pqVtrXCPath());
            this->saveLxc(Lxc.getTlMatrix());
        }
    } else {
        if (this->debugBuildSuperMatrix_) {
            // DEBUG code
            this->log_.info("call DEBUG routine:");
            this->log_.info("build Lxc matrix by supermatrix.");
            this->createEngines<DfOverlapEngine>();
            TlSymmetricMatrix V = this->getSuperMatrix(orbInfo_p);
            this->destroyEngines();
            V.save("fl_Work/debug_V.mat");
            
            TlMatrix L = this->calcCholeskyVectors(V);
            this->saveLxc(L);
        } else {
            // productive code
            this->log_.info("build Lxc matrix by on-the-fly method.");
            const TlRowVectorMatrix2 Lxc 
                = this->calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                     this->getI2pqVtrXCPath());
            this->saveLxc(Lxc.getTlMatrix());
        }

        // check
        // {
        //     this->log_.info("check: LL = L * L^t");
        //     TlMatrix L = this->getLxc();
        //     TlMatrix tL = L;
        //     tL.transpose();
        //     TlMatrix LL = L * tL;
        //     LL.save("fl_Work/debug_LL.mat");
        // }
    }
}



void DfCD::saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath) 
{
    //std::string filepath = this->getI2pqVtrPath();
    std::ofstream ofs;
    ofs.open(filepath.c_str(), std::ofstream::out | std::ofstream::binary);

    const std::size_t size = I2PQ.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
    for (std::size_t i = 0; i < size; ++i) {
        const index_type index1 = I2PQ[i].index1();
        const index_type index2 = I2PQ[i].index2();
        ofs.write(reinterpret_cast<const char*>(&index1), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&index2), sizeof(index_type));
    }

    ofs.close();
}

DfCD::PQ_PairArray DfCD::getI2PQ(const std::string& filepath)
{
    // std::string filepath = this->getI2pqVtrPath();
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        this->log_.critical(TlUtils::format("could not found: %s", filepath.c_str()));
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = Index2(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}

void DfCD::saveLjk(const TlMatrix& Ljk)
{
    DfObject::saveLjkMatrix(Ljk);
}

void DfCD::saveLxc(const TlMatrix& Lxc)
{
    DfObject::saveLxcMatrix(Lxc);
}

TlMatrix DfCD::getLjk()
{
    TlMatrix Ljk = DfObject::getLjkMatrix<TlMatrix>();
    return Ljk;
}

TlMatrix DfCD::getLxc()
{
    TlMatrix Lxc = DfObject::getLxcMatrix<TlMatrix>();
    return Lxc;
}

DfTaskCtrl* DfCD::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    // pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    // pDfTaskCtrl->setCutoffEpsilon_density(0.0);
    // pDfTaskCtrl->setCutoffEpsilon_distribution(this->CDAM_tau_);

    return pDfTaskCtrl;
}

void DfCD::finalize(TlSymmetricMatrix* pMat)
{
    // do nothing
}

void DfCD::finalize(TlSparseSymmetricMatrix *pMat) 
{
    // do nothing
}

void DfCD::finalize_I2PQ(PQ_PairArray *pI2PQ)
{
    std::sort(pI2PQ->begin(), pI2PQ->end());
}

void DfCD::finalize(TlMatrix* pMat)
{
    // do nothing
}

void DfCD::finalize(TlSparseMatrix *pMat) 
{
    // do nothing
}

// void DfCD::finalizeI2PQ_A(PQ_PairArray_A *pI2PQ)
// {
//     std::sort(pI2PQ->begin(), pI2PQ->end());
// }

TlSymmetricMatrix DfCD::getCholeskyVector(const TlVector& L_col,
                                          const PQ_PairArray& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    assert(static_cast<std::size_t>(numOfItilde) == I2PQ.size());

    TlSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

TlMatrix DfCD::getCholeskyVectorA(const TlOrbitalInfoObject& orbInfo_p,
                                  const TlOrbitalInfoObject& orbInfo_q,
                                  const TlVector& L_col,
                                  const PQ_PairArray& I2PQ)
{
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();

    const index_type numOfItilde = L_col.getSize();
    TlMatrix answer(numOfOrbs_p, numOfOrbs_q);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

void DfCD::getJ(TlSymmetricMatrix* pJ)
{
    this->getJ_S(pJ);
}

void DfCD::getJ_S(TlSymmetricMatrix* pJ)
{
    this->log_.info("calc J by CD method (parallel).");
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    TlMatrix L = this->getLjk();
    this->log_.info(TlUtils::format("L(J): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();
    
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        // LI.save(TlUtils::format("fl_Work/debug_LI_J.%d.mat", I));
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();

        *pJ += qi*LI;
    }

    this->finalize(pJ);
}

void DfCD::getJ_A(TlSymmetricMatrix* pJ)
{
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    TlMatrix L = this->getLjk();
    const index_type numOfCBs = L.getNumOfCols();

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlMatrix LI = this->getCholeskyVectorA(orbInfo_p, orbInfo_q,
                                               L.getColVector(I), I2PQ);
        //LI.save(TlUtils::format("fl_Work/debug_LI.%d.mat", I));
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();
        this->log_.info(TlUtils::format("qi [%d] = % f", I, qi));

        *pJ += qi*LI;
    }

    this->finalize(pJ);
}

void DfCD::divideCholeskyBasis(const index_type numOfCBs,
                               index_type *pStart, index_type *pEnd)
{
    *pStart = 0;
    *pEnd = numOfCBs;
}


void DfCD::getK(const RUN_TYPE runType,
                TlSymmetricMatrix *pK)
{
    this->getK_S(runType, pK);
}

void DfCD::getK_S(const RUN_TYPE runType,
                  TlSymmetricMatrix *pK)
{
    TlMatrix L = this->getLjk();
    this->log_.info(TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P;
    if (runType == RUN_RKS) {
        P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    } else {
        P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
    }
    this->log_.info("CD: density matrix");
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    this->log_.info("start loop");
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix l = this->getCholeskyVector(L.getColVector(I), I2PQ);
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
    this->log_.info("finalize");
    this->finalize(pK);
}

void DfCD::getK_A(const RUN_TYPE runType,
                  TlSymmetricMatrix *pK)
{
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);

    TlMatrix L = this->getLjk();
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P = 0.5 * this->getPMatrix(); // RKS
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlMatrix l = this->getCholeskyVectorA(orbInfo_p, 
                                              orbInfo_q,
                                              L.getColVector(I), I2PQ);
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
    this->finalize(pK);
}

void DfCD::getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    if (this->isDedicatedBasisForGridFree_) {
        this->getM_A(P, pM);
    } else {
        this->getM_S(P, pM);
    }
}

void DfCD::getM_S(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    this->log_.info("calc M by CD method. (symmetric routine)");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets"]);
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();

    pM->resize(numOfAOs);

    // cholesky vector
    TlMatrix L = this->getLxc();
    this->log_.info(TlUtils::format("L(xc): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == numOfAOs);
        assert(LI.getNumOfCols() == numOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();

        *pM += qi*LI;
    }

    this->finalize(pM);
}

void DfCD::getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    this->log_.info("calc M by CD method. (asymmetric routine)");
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_sets_GF"]);
    const index_type numOfAOs = orbInfo_p.getNumOfOrbitals();
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlMatrix L = this->getLxc();
    this->log_.info(TlUtils::format("L(xc; A): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();
    
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());

    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlMatrix l = this->getCholeskyVectorA(orbInfo_p,
                                                    orbInfo_q,
                                                    L.getColVector(I), I2PQ);
        // l.save(TlUtils::format("fl_Work/debug_LI_xc_%d.mat", I));
        assert(l.getNumOfRows() == numOfAOs);
        assert(l.getNumOfCols() == dim_M);
        l.transpose();
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        assert(XX.getNumOfRows() == dim_M);
        *pM += XX;
    }

    this->finalize(pM);
}


TlSymmetricMatrix DfCD::getPMatrix()
{
    TlSymmetricMatrix P;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
        break;

    case METHOD_UKS:
        P  = this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_ALPHA, this->m_nIteration -1);
        P += this->getPpqMatrix<TlSymmetricMatrix>(RUN_UKS_BETA, this->m_nIteration -1);
        break;

    case METHOD_ROKS:
        P  = this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSE, this->m_nIteration -1);
        P += this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_OPEN, this->m_nIteration -1);
        break;
        
    default:
        this->log_.critical("program error");
        break;
    }
    return P;
}

TlRowVectorMatrix2 DfCD::calcCholeskyVectorsOnTheFlyS(const TlOrbitalInfoObject& orbInfo,
                                                      const std::string& I2PQ_path)
{
    this->log_.info("call on-the-fly Cholesky Decomposition routine (symmetric)");
    assert(this->pEngines_ != NULL);
    this->initializeCutoffStats(orbInfo.getMaxShellType());

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs +1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));
    TlSparseSymmetricMatrix schwartzTable(numOfAOs);
    PQ_PairArray I2PQ;
    TlVector d; // 対角成分

    assert(this->pEngines_ != NULL);
    this->calcDiagonals(orbInfo, &I2PQ, &schwartzTable, &d);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();

    TlRowVectorMatrix2 L(N, 1, 1, 0, this->isEnableMmap_);

    double error = d.getMaxAbsoluteElement();
    // double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    int progress = 0;
    index_type division =  std::max<index_type>(N * 0.01, 100);
    L.reserve_cols(division);
    index_type m = 0;
    while ((error > threshold) && (m < N)) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif //DEBUG_CD

        // progress 
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e, local mem:%8.1f MB",
                                            m, error, TlSystem::getMaxRSS()));
            ++progress;

            // メモリの確保
            L.reserve_cols(division * progress);
        }
        L.resize(N, m+1);

        // pivot
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }

        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        const index_type pivot_m = pivot[m];
        std::vector<double> G_pm;
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[m+1 +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(orbInfo,
                                                pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        const TlVector L_pm = L.getRowVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getRowVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();
            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp atomic
            L_xm[i] += l_m_pi; // for OpenMP
            
#pragma omp atomic
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            L.set(pivot_i, m, L_xm[i]);
        }

        error = d[pivot[m]];
        // {
        //     error = 0.0;
        //     for (index_type i = m+1; i < N; ++i) {
        //         error += d[pivot[i]];
        //     }
        // }
        ++m;
    }
    L.resize(N, m);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(orbInfo.getMaxShellType());

    return L;
}


TlRowVectorMatrix2 DfCD::calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
                                                      const TlOrbitalInfoObject& orbInfo_q,
                                                      const std::string& I2PQ_path)
{
    this->log_.info("call on-the-fly Cholesky Decomposition routine");
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
    TlVector d; // 対角成分
    this->calcDiagonalsA(orbInfo_p, orbInfo_q,
                         &I2PQ, &schwartzTable, &d);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    const TlVector::size_type N = I2PQ.size();

    TlRowVectorMatrix2 L(N, 1, 1, 0, this->isEnableMmap_);
    //TlMatrix tmpL(N, N);
    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));

    double error = d.getMaxAbsoluteElement();
    // double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (TlVector::size_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    int progress = 0;
    index_type division = std::max<index_type>(N * 0.01, 100);
    L.reserve_cols(division);
    index_type m = 0;
    while (error > threshold) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif //DEBUG_CD

        // progress 
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e, local mem:%8.1f MB",
                                            m, error, TlSystem::getMaxRSS()));
            ++progress;

            // メモリの確保
            L.reserve_cols(division * progress);
        }
        L.resize(N, m+1);

        // pivot
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }

        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        //tmpL.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // calc
        const index_type pivot_m = pivot[m];
        std::vector<double> G_pm;
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[m+1 +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            G_pm = this->getSuperMatrixElementsA(orbInfo_p, orbInfo_q,
                                                 pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        const TlVector L_pm = L.getRowVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getRowVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();

            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp atomic
            L_xm[i] += l_m_pi; // for OpenMP
            
#pragma omp atomic
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            L.set(pivot_i, m, L_xm[i]);
            //tmpL.set(pivot_i, m, L_xm[i]);
        }

        // error計算
        error = d[pivot[m]];
        // {
        //     error = 0.0;
        //     for (std::size_t i = m +1; i < N; ++i) {
        //         error += d[pivot[i]];
        //     }
        // }
        ++m;
    }
    L.resize(N, m);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    return L;
}


void DfCD::calcDiagonals(const TlOrbitalInfoObject& orbInfo,
                         PQ_PairArray *pI2PQ,
                         TlSparseSymmetricMatrix *pSchwartzTable,
                         TlVector *pDiagonals)
{
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs + 1) / 2;

    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(numOfPQs);
    pSchwartzTable->clear();
    pSchwartzTable->resize(numOfAOs);
    TlSparseSymmetricMatrix diagonalMat(numOfAOs);

    // task
    this->log_.info("diagonal calculation: start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_kernel(orbInfo,
                                   taskList,
                                   pSchwartzTable,
                                   &diagonalMat, pI2PQ);
        hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    this->finalize_I2PQ(pI2PQ);
    this->finalize(&diagonalMat);
    this->finalize(pSchwartzTable);

    // set diagonals
    const index_type numOfI = pI2PQ->size();
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PQ)[i].index1();
        const index_type col = (*pI2PQ)[i].index2();
        const double value = diagonalMat.get(row, col);
        (*pDiagonals)[i] = value;
    }
}


void DfCD::calcDiagonalsA(const TlOrbitalInfoObject& orbInfo_p,
                          const TlOrbitalInfoObject& orbInfo_q,
                          PQ_PairArray *pI2PQ,
                          TlSparseMatrix *pSchwartzTable,
                          TlVector *pDiagonals)
{
    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();
    const index_type numOfPQs = numOfOrbs_p * numOfOrbs_q;
    pI2PQ->clear();
    pI2PQ->reserve(numOfPQs);
    pSchwartzTable->clear();
    pSchwartzTable->resize(numOfOrbs_p, numOfOrbs_q);
    TlSparseMatrix diagonalMat(numOfOrbs_p, numOfOrbs_q);

    // TODO: pq のスクリーニング

    // task
    this->log_.info("diagonal calculation: start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbInfo_p,
                                          orbInfo_q,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonalsA_kernel(orbInfo_p,
                                    orbInfo_q,
                                    taskList,
                                    pI2PQ,
                                    pSchwartzTable,
                                    &diagonalMat);
        hasTask = pDfTaskCtrl->getQueue2(orbInfo_p,
                                         orbInfo_q,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    this->finalize_I2PQ(pI2PQ);
    this->finalize(&diagonalMat);
    this->finalize(pSchwartzTable);

    // set diagonals
    const index_type numOfI = pI2PQ->size();
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PQ)[i].index1();
        const index_type col = (*pI2PQ)[i].index2();
        const double value = diagonalMat.get(row, col);
        (*pDiagonals)[i] = value;
    }
}


void DfCD::calcDiagonals_kernel(const TlOrbitalInfoObject& orbInfo,
                                const std::vector<DfTaskCtrl::Task2>& taskList,
                                TlSparseSymmetricMatrix *pSchwartzTable,
                                TlSparseSymmetricMatrix *pDiagonalMat,
                                PQ_PairArray *pI2PQ)
{
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseSymmetricMatrix local_diagMat(numOfAOs);
        TlSparseSymmetricMatrix local_schwartzTable(numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        assert(0 <= threadID);
        assert(threadID < this->numOfThreads_);

        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            assert(this->pEngines_[threadID] != NULL);
            this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                            0, orbInfo, shellIndexQ,
                                            0, orbInfo, shellIndexP,
                                            0, orbInfo, shellIndexQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                        const int pq_index = p * maxStepsQ + q;
                        const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                        
                        const double value = this->pEngines_[threadID]->value(pqpq_index);

                        // for schwartz
                        maxValue = std::max(maxValue, std::fabs(value));
                        
                        // for I~ to pq table
                        if (std::fabs(value) > tau) {
                            if (value > 0) {
                                local_diagMat.set(indexP, indexQ, value);
                                local_I2PQ.push_back(Index2(indexP, indexQ));
                            } else {
                                this->log_.warn(TlUtils::format("pqpq_value: (%d %d)=% e is spoiled.",
                                                                indexP, indexQ, value));
                            }
                        }
                    }
                }
            }
            local_schwartzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
#pragma omp critical(DfCD__calcDiagonals_kernel_1)
            {
                pI2PQ->insert(pI2PQ->end(),
                              local_I2PQ.begin(), local_I2PQ.end());
            }
#pragma omp critical(DfCD__calcDiagonals_kernel_2)
            {
                pDiagonalMat->merge(local_diagMat);
            }
#pragma omp critical(DfCD__calcDiagonals_kernel_3)
            {
                pSchwartzTable->merge(local_schwartzTable);
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
            *pSchwartzTable = local_schwartzTable;
        }
#endif // _OPENMP
    }
}

void DfCD::calcDiagonalsA_kernel(const TlOrbitalInfoObject& orbInfo_p,
                                 const TlOrbitalInfoObject& orbInfo_q,
                                 const std::vector<DfTaskCtrl::Task2>& taskList,
                                 PQ_PairArray *pI2PQ,
                                 TlSparseMatrix *pSchwartzTable,
                                 TlSparseMatrix *pDiagonalMat)
{
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();
    pDiagonalMat->resize(numOfOrbs_p, numOfOrbs_q);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseMatrix local_diagMat(numOfOrbs_p, numOfOrbs_q);
        TlSparseMatrix local_schwartzTable(numOfOrbs_p, numOfOrbs_q);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbInfo_p.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo_q.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            this->pEngines_[threadID]->calc(0, orbInfo_p, shellIndexP,
                                            0, orbInfo_q, shellIndexQ,
                                            0, orbInfo_p, shellIndexP,
                                            0, orbInfo_q, shellIndexQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    const int pq_index = p * maxStepsQ + q;
                    const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                    
                    const double value = this->pEngines_[threadID]->value(pqpq_index);
                    
                    // for schwartz
                    maxValue = std::max(maxValue, std::fabs(value));
                    
                    // for I~ to pq table
                    if (std::fabs(value) > tau) {
                        if (value > 0) {
                            local_diagMat.set(indexP, indexQ, value);
                            local_I2PQ.push_back(Index2(indexP, indexQ));
                        } else {
                            this->log_.warn(TlUtils::format("pqpq_value: (%d %d)=% e is spoiled.",
                                                            indexP, indexQ, value));
                        }
                    }
                }
            }
            local_schwartzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
#pragma omp critical(DfCD__calcDiagonalsA_kernel_1)
            {
                pI2PQ->insert(pI2PQ->end(),
                              local_I2PQ.begin(), local_I2PQ.end());
            }
#pragma omp critical(DfCD__calcDiagonalsA_kernel_2)
            {
                pDiagonalMat->merge(local_diagMat);
            }
#pragma omp critical(DfCD__calcDiagonalsA_kernel_3)
            {
                pSchwartzTable->merge(local_schwartzTable);
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
            *pSchwartzTable = local_schwartzTable;
        }
#endif // _OPENMP
    }
}


bool DfCD::isAliveBySchwartzCutoff(const index_type shellIndexP,
                                   const index_type shellIndexQ,
                                   const index_type shellIndexR,
                                   const index_type shellIndexS,
                                   const int shellQuartetType,
                                   const TlSparseMatrix& schwarzTable,
                                   const double threshold)
{
    bool answer = false;

    const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
    const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

    if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
        answer = true;

#pragma omp atomic
        ++(this->cutoffAlive_schwartz_[shellQuartetType]);
    }

#pragma omp atomic
    ++(this->cutoffAll_schwartz_[shellQuartetType]);

    return answer;
}


void DfCD::initializeCutoffStats(const int maxShellType)
{
    // clear cutoff stats
    const int numOfShellPairType = maxShellType* maxShellType;
    const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;
    this->cutoffAll_schwartz_.clear();
    this->cutoffAlive_schwartz_.clear();
    this->cutoffAll_schwartz_.resize(numOfShellQuartetType, 0);
    this->cutoffAlive_schwartz_.resize(numOfShellQuartetType, 0);
}


void DfCD::schwartzCutoffReport(const int maxShellType)
{
    // const int maxShellType = this->maxShellType_;
    std::vector<std::string> typeStr4(maxShellType * maxShellType * maxShellType * maxShellType);
    {
        static const char typeChar[] = "SPDFG";
        std::string tmp(4, 'X');
        int index = 0;
        for (int i = 0; i < maxShellType; ++i) {
            tmp[0] = typeChar[i];
            for (int j = 0; j < maxShellType; ++j) {
                tmp[1] = typeChar[j];
                for (int k = 0; k < maxShellType; ++k) {
                    tmp[2] = typeChar[k];
                    for (int l = 0; l < maxShellType; ++l) {
                        tmp[3] = typeChar[l];
                        typeStr4[index] = tmp;
                        ++index;
                    }
                }
            }
        }
    }

    // static const char typeStr4[][5] = {
    //     "SSSS", "SSSP", "SSSD", "SSPS", "SSPP", "SSPD", "SSDS", "SSDP", "SSDD",
    //     "SPSS", "SPSP", "SPSD", "SPPS", "SPPP", "SPPD", "SPDS", "SPDP", "SPDD",
    //     "SDSS", "SDSP", "SDSD", "SDPS", "SDPP", "SDPD", "SDDS", "SDDP", "SDDD",
    //     "PSSS", "PSSP", "PSSD", "PSPS", "PSPP", "PSPD", "PSDS", "PSDP", "PSDD",
    //     "PPSS", "PPSP", "PPSD", "PPPS", "PPPP", "PPPD", "PPDS", "PPDP", "PPDD",
    //     "PDSS", "PDSP", "PDSD", "PDPS", "PDPP", "PDPD", "PDDS", "PDDP", "PDDD",
    //     "DSSS", "DSSP", "DSSD", "DSPS", "DSPP", "DSPD", "DSDS", "DSDP", "DSDD",
    //     "DPSS", "DPSP", "DPSD", "DPPS", "DPPP", "DPPD", "DPDS", "DPDP", "DPDD",
    //     "DDSS", "DDSP", "DDSD", "DDPS", "DDPP", "DDPD", "DDDS", "DDDP", "DDDD",
    // };

    // cutoff report for schwarz
    bool hasCutoffSchwarz = false;
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
            for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                    const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                    if (this->cutoffAll_schwartz_[shellTypeABCD] != 0) {
                        hasCutoffSchwarz = true;
                        break;
                    }
                }
            }
        }
    }
    if (hasCutoffSchwarz == true) {
        this->log_.info("schwarz cutoff report");
        this->log_.info(TlUtils::format("threshold: %e", this->CDAM_tau_));
        this->log_.info("type: alive / all (ratio)");
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
                for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                    const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                    for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                        const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                        
                        if (this->cutoffAll_schwartz_[shellTypeABCD] > 0) {
                            const double ratio = (double)this->cutoffAlive_schwartz_[shellTypeABCD]
                                / (double)this->cutoffAll_schwartz_[shellTypeABCD]
                                * 100.0;
                            this->log_.info(TlUtils::format(" %4s: %12ld / %12ld (%6.2f%%)",
                                                            typeStr4[shellTypeABCD].c_str(),
                                                            this->cutoffAlive_schwartz_[shellTypeABCD],
                                                            this->cutoffAll_schwartz_[shellTypeABCD],
                                                            ratio));
                        }
                    }
                }
            }
        }
    }
}


std::vector<double>
DfCD::getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                             const index_type G_row,
                             const std::vector<index_type>& G_col_list,
                             const PQ_PairArray& I2PQ,
                             const TlSparseSymmetricMatrix& schwartzTable)
{
    this->ERI_cache_.clear();

    const std::vector<Index4> calcList = this->getCalcList(orbInfo, G_row, G_col_list, I2PQ);
    this->calcERIs(orbInfo, calcList, schwartzTable);
    const std::vector<double> answer = this->setERIs(orbInfo, G_row, G_col_list, I2PQ);

    return answer;
}

std::vector<double>
DfCD::getSuperMatrixElementsA(const TlOrbitalInfoObject& orbInfo_p,
                              const TlOrbitalInfoObject& orbInfo_q,
                              const index_type G_row,
                              const std::vector<index_type>& G_col_list,
                              const PQ_PairArray& I2PQ,
                              const TlSparseMatrix& schwartzTable)
{
    this->ERI_cache_A_.clear();

    const std::vector<IndexPair4A> calcList = this->getCalcListA(orbInfo_p, orbInfo_q,
                                                                 G_row, G_col_list, I2PQ);
    this->calcERIsA(orbInfo_p, orbInfo_q,
                    calcList, schwartzTable);
    const std::vector<double> answer = this->setERIsA(orbInfo_p, orbInfo_q,
                                                      G_row, G_col_list, I2PQ);

    return answer;
}


std::vector<DfCD::Index4> 
DfCD::getCalcList(const TlOrbitalInfoObject& orbInfo,
                  const index_type G_row,
                  const std::vector<index_type>& G_col_list,
                  const PQ_PairArray& I2PQ)
{
    std::set<Index4> calcSet;

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();
    const index_type shellIndexP = orbInfo.getShellIndex(indexP);
    const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);

    const index_type numOf_G_cols = G_col_list.size();
#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();
        const index_type shellIndexR = orbInfo.getShellIndex(indexR);
        const index_type shellIndexS = orbInfo.getShellIndex(indexS);

        Index4 index4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList) 
        {
            calcSet.insert(index4);
        }
    }

    std::vector<Index4> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}

std::vector<DfCD::IndexPair4A> 
DfCD::getCalcListA(const TlOrbitalInfoObject& orbInfo_p,
                   const TlOrbitalInfoObject& orbInfo_q,
                   const index_type G_row,
                   const std::vector<index_type>& G_col_list,
                   const PQ_PairArray& I2PQ)
{
    std::set<IndexPair4A> calcSet;

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();
    const index_type shellIndexP = orbInfo_p.getShellIndex(indexP);
    const index_type shellIndexQ = orbInfo_q.getShellIndex(indexQ);

    const index_type numOf_G_cols = G_col_list.size();
#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();
        const index_type shellIndexR = orbInfo_p.getShellIndex(indexR);
        const index_type shellIndexS = orbInfo_q.getShellIndex(indexS);

        IndexPair4A indexPair4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList) 
        {
            calcSet.insert(indexPair4);
        }
    }

    std::vector<IndexPair4A> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}


void DfCD::calcERIs(const TlOrbitalInfoObject& orbInfo,
                    const std::vector<Index4>& calcList,
                    const TlSparseSymmetricMatrix& schwartzTable) 
{
    const int maxShellType = orbInfo.getMaxShellType();
    const double threshold = this->CDAM_tau_;
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < numOfList; ++i) {
            const index_type shellIndexP = calcList[i].index1();
            const index_type shellIndexQ = calcList[i].index2();
            const index_type shellIndexR = calcList[i].index3();
            const index_type shellIndexS = calcList[i].index4();
            
            const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            const int shellTypeS = orbInfo.getShellType(shellIndexS);
            
            const int shellQuartetType =
                ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
            const bool isAlive = this->isAliveBySchwartzCutoff(shellIndexP, shellIndexQ,
                                                               shellIndexR, shellIndexS,
                                                               shellQuartetType,
                                                               schwartzTable,
                                                               threshold);
            if (isAlive == true) {
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                
                this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexQ,
                                                0, orbInfo, shellIndexR,
                                                0, orbInfo, shellIndexS);
                
                const int steps = maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                std::vector<double> buf(steps);
                for (int count = 0; count < steps; ++count) {
                    buf[count] = this->pEngines_[threadID]->value(count);
                }

                local_cache[calcList[i]] = buf;
            }
        }

        // merge cache
#pragma omp critical(DfCD__calcERIs)
        {
            this->ERI_cache_.insert(local_cache.begin(), local_cache.end());
        }
    }
}

void DfCD::calcERIsA(const TlOrbitalInfoObject& orbInfo_p,
                     const TlOrbitalInfoObject& orbInfo_q,
                     const std::vector<IndexPair4A>& calcList,
                     const TlSparseMatrix& schwartzTable) 
{
    const int maxShellType = orbInfo_p.getMaxShellType();
    assert(maxShellType == orbInfo_q.getMaxShellType());
    // const double threshold = this->CDAM_tau_;
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType_A local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < numOfList; ++i) {
            const index_type shellIndexP = calcList[i].index1();
            const index_type shellIndexQ = calcList[i].index2();
            const index_type shellIndexR = calcList[i].index3();
            const index_type shellIndexS = calcList[i].index4();
            
            const int shellTypeP = orbInfo_p.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo_q.getShellType(shellIndexQ);
            const int shellTypeR = orbInfo_p.getShellType(shellIndexR);
            const int shellTypeS = orbInfo_q.getShellType(shellIndexS);
            
            // const int shellQuartetType =
            //     ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
            // const bool isAlive = this->isAliveBySchwartzCutoff(shellIndexP, shellIndexQ,
            //                                                    shellIndexR, shellIndexS,
            //                                                    shellQuartetType,
            //                                                    schwartzTable,
            //                                                    threshold);
            const bool isAlive = true;
            if (isAlive == true) {
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                this->pEngines_[threadID]->calc(0, orbInfo_p, shellIndexP,
                                                0, orbInfo_q, shellIndexQ,
                                                0, orbInfo_p, shellIndexR,
                                                0, orbInfo_q, shellIndexS);
                
                const int steps = maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                std::vector<double> buf(steps);
                for (int count = 0; count < steps; ++count) {
                    buf[count] = this->pEngines_[threadID]->value(count);
                }

                local_cache[calcList[i]] = buf;
            }
        }

        // merge cache
#pragma omp critical(DfCD__calcERIs)
        {
            this->ERI_cache_A_.insert(local_cache.begin(), local_cache.end());
        }
    }
}


std::vector<double>
DfCD::setERIs(const TlOrbitalInfoObject& orbInfo,
              const index_type G_row,
              const std::vector<index_type> G_col_list,
              const PQ_PairArray& I2PQ)
{
    const index_type indexP_orig = I2PQ[G_row].index1();
    const index_type indexQ_orig = I2PQ[G_row].index2();
    const index_type shellIndexP_orig = orbInfo.getShellIndex(indexP_orig);
    const index_type shellIndexQ_orig = orbInfo.getShellIndex(indexQ_orig);

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);
#pragma omp parallel for
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        index_type indexP = indexP_orig;
        index_type indexQ = indexQ_orig;
        index_type shellIndexP = shellIndexP_orig;
        index_type shellIndexQ = shellIndexQ_orig;

        const index_type G_col = G_col_list[i];

        index_type indexR = I2PQ[G_col].index1();
        index_type indexS = I2PQ[G_col].index2();
        index_type shellIndexR = orbInfo.getShellIndex(indexR);
        index_type shellIndexS = orbInfo.getShellIndex(indexS);

        std::vector<double> values;
#pragma omp critical(DfCD__setERIs)
        {
            assert(this->ERI_cache_.find(Index4(shellIndexP, shellIndexQ,
                                                shellIndexR, shellIndexS)) != this->ERI_cache_.end());
            values = this->ERI_cache_[Index4(shellIndexP, shellIndexQ,
                                                 shellIndexR, shellIndexS)];
        }

        if (values.empty() != true) {
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            const int basisTypeR = indexR - shellIndexR;
            const int basisTypeS = indexS - shellIndexS;
            
            // const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            const int shellTypeS = orbInfo.getShellType(shellIndexS);
            // const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            
            const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
            assert(static_cast<int>(values.size()) > index);
            
#pragma omp atomic
            answer[i] += values.at(index);
        }
    }

    return answer;
}

std::vector<double>
DfCD::setERIsA(const TlOrbitalInfoObject& orbInfo_p,
               const TlOrbitalInfoObject& orbInfo_q,
               const index_type G_row,
               const std::vector<index_type> G_col_list,
               const PQ_PairArray& I2PQ)
{
    const index_type indexP_orig = I2PQ[G_row].index1();
    const index_type indexQ_orig = I2PQ[G_row].index2();
    const index_type shellIndexP_orig = orbInfo_p.getShellIndex(indexP_orig);
    const index_type shellIndexQ_orig = orbInfo_q.getShellIndex(indexQ_orig);

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);
#pragma omp parallel for
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        index_type indexP = indexP_orig;
        index_type indexQ = indexQ_orig;
        index_type shellIndexP = shellIndexP_orig;
        index_type shellIndexQ = shellIndexQ_orig;

        const index_type G_col = G_col_list[i];

        index_type indexR = I2PQ[G_col].index1();
        index_type indexS = I2PQ[G_col].index2();
        index_type shellIndexR = orbInfo_p.getShellIndex(indexR);
        index_type shellIndexS = orbInfo_q.getShellIndex(indexS);

        std::vector<double> values;
#pragma omp critical(DfCD__setERIs)
        {
            assert(this->ERI_cache_A_.find(IndexPair4A(shellIndexP, shellIndexQ,
                                                       shellIndexR, shellIndexS)) != this->ERI_cache_A_.end());
            values = this->ERI_cache_A_[IndexPair4A(shellIndexP, shellIndexQ,
                                                    shellIndexR, shellIndexS)];
        }

        assert(values.empty() != true);
        {
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            const int basisTypeR = indexR - shellIndexR;
            const int basisTypeS = indexS - shellIndexS;
            
            // const int shellTypeP = orbInfo_p.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo_q.getShellType(shellIndexQ);
            const int shellTypeR = orbInfo_p.getShellType(shellIndexR);
            const int shellTypeS = orbInfo_q.getShellType(shellIndexS);
            // const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            
            const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
            assert(static_cast<int>(values.size()) > index);
            
#pragma omp atomic
            answer[i] += values.at(index);
        }
    }

    return answer;
}

////////////////////////////////////////////////////////////////////////////////
// for DEBUG 
// 
////////////////////////////////////////////////////////////////////////////////
TlSymmetricMatrix DfCD::getSuperMatrix(const TlOrbitalInfoObject& orbInfo) 
{
    const index_type numOfOrbs = orbInfo.getNumOfOrbitals();

    PQ_PairArray I2PQ;
    TlSparseSymmetricMatrix schwartzTable(numOfOrbs);
    TlVector d; // 対角成分
    this->calcDiagonals(orbInfo,
                        &I2PQ, &schwartzTable, &d);
    this->saveI2PQ(I2PQ, this->getI2pqVtrXCPath());
    const std::size_t N = I2PQ.size();

    TlSymmetricMatrix V(N);
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t i = 0; i < N; ++i) {
            const index_type indexP = I2PQ[i].index1();
            const index_type indexQ = I2PQ[i].index2();
            const index_type shellIndexP = orbInfo.getShellIndex(indexP);
            const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            // const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
            //const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            for (std::size_t j = 0; j <= i; ++j) {
                const index_type indexR = I2PQ[j].index1();
                const index_type indexS = I2PQ[j].index2();
                const index_type shellIndexR = orbInfo.getShellIndex(indexR);
                const index_type shellIndexS = orbInfo.getShellIndex(indexS);
                const int basisTypeR = indexR - shellIndexR;
                const int basisTypeS = indexS - shellIndexS;
                const int shellTypeR = orbInfo.getShellType(shellIndexR);
                const int shellTypeS = orbInfo.getShellType(shellIndexS);
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                
                this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexQ,
                                                0, orbInfo, shellIndexR,
                                                0, orbInfo, shellIndexS);
                const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
                const double value = this->pEngines_[threadID]->value(index);
                V.set(i, j, value);
            }
        }
    }

    return V;
}

TlSymmetricMatrix DfCD::getSuperMatrix(const TlOrbitalInfoObject& orbInfo_p,
                                       const TlOrbitalInfoObject& orbInfo_q) 
{
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();

    PQ_PairArray I2PQ;
    TlSparseMatrix schwartzTable(numOfOrbs_p, numOfOrbs_q);
    TlVector d; // 対角成分
    this->calcDiagonalsA(orbInfo_p, orbInfo_q,
                         &I2PQ, &schwartzTable, &d);
    this->saveI2PQ(I2PQ, this->getI2pqVtrXCPath());
    const std::size_t N = I2PQ.size();

    TlSymmetricMatrix V(N);
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t i = 0; i < N; ++i) {
            const index_type indexP = I2PQ[i].index1();
            const index_type indexQ = I2PQ[i].index2();
            const index_type shellIndexP = orbInfo_p.getShellIndex(indexP);
            const index_type shellIndexQ = orbInfo_q.getShellIndex(indexQ);
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            // const int shellTypeP = orbInfo_p.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo_q.getShellType(shellIndexQ);
            //const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            for (std::size_t j = 0; j <= i; ++j) {
                const index_type indexR = I2PQ[j].index1();
                const index_type indexS = I2PQ[j].index2();
                const index_type shellIndexR = orbInfo_p.getShellIndex(indexR);
                const index_type shellIndexS = orbInfo_q.getShellIndex(indexS);
                const int basisTypeR = indexR - shellIndexR;
                const int basisTypeS = indexS - shellIndexS;
                const int shellTypeR = orbInfo_p.getShellType(shellIndexR);
                const int shellTypeS = orbInfo_q.getShellType(shellIndexS);
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                
                this->pEngines_[threadID]->calc(0, orbInfo_p, shellIndexP,
                                                0, orbInfo_q, shellIndexQ,
                                                0, orbInfo_p, shellIndexR,
                                                0, orbInfo_q, shellIndexS);
                const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
                const double value = this->pEngines_[threadID]->value(index);
                V.set(i, j, value);
            }
        }
    }

    return V;
}

TlMatrix DfCD::calcCholeskyVectors(const TlSymmetricMatrix& V)
{
    const index_type N = V.getNumOfRows();
    TlVector d(N); // 対角成分
    for (index_type i = 0; i < N; ++i) {
        d[i] = V.get(i, i);
    }
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    TlMatrix L(N, N);
    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));

    index_type m = 0;
    while (error > threshold) {
        // pivot
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }

        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        const index_type pivot_m = pivot[m];
        const index_type numOf_G_cols = N -(m+1);

        // CD calc
        const TlVector L_pm = L.getRowVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getRowVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();
            // const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;
            const double l_m_pi = (V.get(pivot_m, pivot[m+1 + i]) - sum_ll) * inv_l_m_pm;

#pragma omp atomic
            L_xm[i] += l_m_pi; // for OpenMP
            
#pragma omp atomic
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            L.set(pivot_i, m, L_xm[i]);
        }

        error = 0.0;
        for (int i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }

        ++m;
    }
    L.resize(N, m);

    return L;
}
