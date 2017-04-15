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
#include "CnError.h"
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
    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut_value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut_value"].getDouble();
    }    
    this->cutoffThreshold_primitive_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_threshold_primitive"].getStr().empty() != true) {
        this->cutoffThreshold_primitive_ = (*pPdfParam)["cutoff_threshold_primitive"].getDouble();
    }    

    this->CDAM_tau_ = 1.0E-10;
    if ((*pPdfParam)["CDAM_tau"].getStr().empty() != true) {
        this->CDAM_tau_ = (*pPdfParam)["CDAM_tau"].getDouble();
    }    

    this->epsilon_ = 1.0E-4;
    if ((*pPdfParam)["CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["CD_epsilon"].getDouble();
    }    

    this->CDAM_tau_K_ = 1.0E-4;
    if ((*pPdfParam)["CDAM_tau_K"].getStr().empty() != true) {
        this->CDAM_tau_K_ = (*pPdfParam)["CDAM_tau_K"].getDouble();
    }    

    this->epsilon_K_ = 1.0E-4;
    if ((*pPdfParam)["CD_epsilon_K"].getStr().empty() != true) {
        this->epsilon_K_ = (*pPdfParam)["CD_epsilon_K"].getDouble();
    }    

    
    if (this->K_engine_ == DfObject::K_ENGINE_CD) {
        this->fastCDK_mode_ = FASTCDK_NONE;
    } else {
        assert(this->K_engine_ == DfObject::K_ENGINE_FASTCDK);
        this->fastCDK_mode_ = FASTCDK_PRODUCTIVE;

        if ((*pPdfParam)["debug/DfCD/FastCDK_mode"].getStr().empty() != true) {
            const std::string fastCDK_mode_str = TlUtils::toUpper((*pPdfParam)["debug/DfCD/FastCDK_mode"].getStr());
            if (fastCDK_mode_str == "PRODUCTIVE_FULL") {
                this->fastCDK_mode_ = FASTCDK_PRODUCTIVE_FULL;
            } else if (fastCDK_mode_str == "FULL_SUPERMATRIX") {
                this->fastCDK_mode_ = FASTCDK_DEBUG_FULL_SUPERMATRIX;
            } else if (fastCDK_mode_str == "SUPERMATRIX") {
                this->fastCDK_mode_ = FASTCDK_DEBUG_SUPERMATRIX;
            }
        }    
    }
    
    this->debugBuildSuperMatrix_ = false;
    if ((*pPdfParam)["debug/DfCD/build_supermatrix"].getStr().empty() != true) {
        this->debugBuildSuperMatrix_ = (*pPdfParam)["debug/DfCD/build_supermatrix"].getBoolean();
    }    

    this->debugCheckCD_ = false;
    if ((*pPdfParam)["debug/DfCD/check_CD"].getStr().empty() != true) {
        this->debugCheckCD_ = (*pPdfParam)["debug/DfCD/check_CD"].getBoolean();
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
    this->log_.info("calc CholeskyVectors (serial)");

    // for J & K
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_set"]);
    if (this->debugBuildSuperMatrix_) {
        // DEBUG code
        this->log_.info("call DEBUG routine:");
        this->log_.info("build L matrix by supermatrix.");

        TlSymmetricMatrix V;
        {
            PQ_PairArray I2PQ;
            this->createEngines<DfEriEngine>();
            V = this->getSuperMatrix(orbInfo, &I2PQ);
            this->destroyEngines();
            V.save("fl_Work/debug_Vjk.mat");
            this->saveI2PQ(I2PQ, this->getI2pqVtrPath());
        }
            
        TlMatrix L = this->calcCholeskyVectors(V);

        // check CD
        if (this->debugCheckCD_) {
            TlMatrix tL = L;
            tL.transpose();
            TlMatrix LL = L * tL;
            LL.save("fl_Work/debug_LL.mat");
        }

        this->saveLjk(TlRowVectorMatrix(L));
        this->debugOutLjk(L); // debug
        this->log_.info("");
    } else {
        // productive code
        this->log_.info("CD (J) routine: start");
        this->createEngines<DfEriEngine>();

        const TlRowVectorMatrix Ljk 
            = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                     this->getI2pqVtrPath(),
                                                     this->epsilon_,
                                                     &DfCD::calcDiagonals,
                                                     &DfCD::getSuperMatrixElements);
        
        this->saveLjk(Ljk);
        // this->debugOutLjk(Ljk.getTlMatrixObject()); // debug

        this->destroyEngines();
        this->log_.info("CD (J) routine: end");
    }

}


// K only
void DfCD::calcCholeskyVectorsForK()
{
    this->log_.info("calc CholeskyVectors for K (serial)");
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_set"]);

    switch (this->fastCDK_mode_) {
    case FASTCDK_DEBUG_FULL_SUPERMATRIX:
        {
            this->log_.info("fast CDK routine (using debug FULL supermatrix) is valided.");

            TlSymmetricMatrix V;
            {
                PQ_PairArray I2PQ;
                
                this->createEngines<DfEriEngine>();
                V = this->getSuperMatrix_K_full(orbInfo, &I2PQ);
                this->destroyEngines();
                
                V.save("fl_Work/debug_Vk.mat");
                this->saveI2PQ(I2PQ, this->getI2prVtrPath());
            }

            TlMatrix L = this->calcCholeskyVectors(V);
            
            // check CD
            if (this->debugCheckCD_) {
                TlMatrix tL = L;
                tL.transpose();
                TlMatrix LL = L * tL;
                LL.save("fl_Work/debug_LL_K.mat");
            }
            
            this->saveLk(TlRowVectorMatrix(L));
            this->debugOutLk(L); // debug
            this->log_.info("");
        }
        break;

    case FASTCDK_DEBUG_SUPERMATRIX:
        {
            this->log_.info("fast CDK routine (using debug supermatrix) is valided.");
            TlSymmetricMatrix V;
            {
                PQ_PairArray I2PQ;
                this->createEngines<DfEriEngine>();
                V = this->getSuperMatrix_K_half(orbInfo, &I2PQ);
                this->destroyEngines();
                
                V.save("fl_Work/debug_Vk.mat");
                this->saveI2PQ(I2PQ, this->getI2prVtrPath());
            }

            TlMatrix L = this->calcCholeskyVectors(V);
            
            if (this->debugCheckCD_) {
                // check CD
                TlMatrix tL = L;
                tL.transpose();
                TlMatrix LL = L * tL;
                LL.save("fl_Work/debug_LL_K.mat");
            }
            
            this->saveLk(TlRowVectorMatrix(L));
            this->debugOutLk(L); // debug
            this->log_.info("");
        }
        break;

    case FASTCDK_PRODUCTIVE_FULL:
        {
            this->log_.info("fast CDK routine(full).");
            this->createEngines<DfEriEngine>();

            const TlRowVectorMatrix Lk 
                = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                         this->getI2prVtrPath(),
                                                         this->epsilon_K_,
                                                         &DfCD::calcDiagonals_K_full,
                                                         &DfCD::getSuperMatrixElements_K_full);
            this->saveLk(Lk);
            // this->debugOutLk(Lk.getTlMatrixObject()); // debug

            if (this->debugCheckCD_) {
                // check CD
                TlMatrix mLk = Lk.getTlMatrixObject();
                TlMatrix tmLk = mLk;
                tmLk.transpose();
                TlMatrix LL = mLk * tmLk;
                LL.save("fl_Work/debug_LL_K.mat");
            }
            
            this->destroyEngines();
            this->log_.info("");
        }
        break;
        
    case FASTCDK_PRODUCTIVE:
        {
            this->log_.info("CD (K) routine: start");
            this->createEngines<DfEriEngine>();

            const TlRowVectorMatrix Lk 
                = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo,
                                                         this->getI2prVtrPath(),
                                                         this->epsilon_K_,
                                                         &DfCD::calcDiagonals_K_half,
                                                         &DfCD::getSuperMatrixElements_K_half);
            this->saveLk(Lk);
            // this->debugOutLk(Lk.getTlMatrixObject()); // debug

            if (this->debugCheckCD_) {
                // check CD
                TlMatrix mLk = Lk.getTlMatrixObject();
                TlMatrix tmLk = mLk;
                tmLk.transpose();
                TlMatrix LL = mLk * tmLk;
                LL.save("fl_Work/debug_LL_K.mat");
            }

            this->destroyEngines();
            this->log_.info("CD (K) routine: end");
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

void DfCD::calcCholeskyVectorsForGridFree()
{
    // for XC(gridfree)
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);

    if (this->isDedicatedBasisForGridFree_) {
        const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                      (*this->pPdfParam_)["basis_set_gridfree"]);
        if (this->debugBuildSuperMatrix_) {
            // DEBUG code
            this->log_.info("call DEBUG routine:");
            this->log_.info("build Lxc matrix by supermatrix.");

            TlSymmetricMatrix V;
            {
                PQ_PairArray I2PQ;
                this->createEngines<DfOverlapEngine>();
                V = this->getSuperMatrix(orbInfo_p, orbInfo_q, &I2PQ);
                this->destroyEngines();
                V.save("fl_Work/debug_Vxc.mat");
                this->saveI2PQ(I2PQ, this->getI2pqVtrXCPath());
            }
            
            TlMatrix L = this->calcCholeskyVectors(V);
            this->saveLxc(TlRowVectorMatrix(L));
            this->debugOutLxc(L); // debug
        } else {
            // productive code
            const TlRowVectorMatrix Lxc 
                = this->calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p,
                                                                     orbInfo_q,
                                                                     this->getI2pqVtrXCPath());
            this->saveLxc(Lxc);
        }
    } else {
        if (this->debugBuildSuperMatrix_) {
            // DEBUG code
            this->log_.info("call DEBUG routine:");
            this->log_.info("build Lxc matrix by supermatrix.");

            TlSymmetricMatrix V;
            {
                PQ_PairArray I2PQ;
                this->createEngines<DfOverlapEngine>();
                V = this->getSuperMatrix(orbInfo_p, &I2PQ);
                this->destroyEngines();
                V.save("fl_Work/debug_Vxc.mat");
                this->saveI2PQ(I2PQ, this->getI2pqVtrXCPath());
            }
            
            TlMatrix L = this->calcCholeskyVectors(V);
            this->saveLxc(TlRowVectorMatrix(L));
            this->debugOutLxc(L); // debug
        } else {
            // productive code
            this->log_.info("build Lxc matrix by on-the-fly method.");
            this->createEngines<DfOverlapEngine>();

            const TlRowVectorMatrix Lxc
                = this->calcCholeskyVectorsOnTheFlyS_new(orbInfo_p,
                                                         this->getI2pqVtrXCPath(),
                                                         this->epsilon_,
                                                         &DfCD::calcDiagonals,
                                                         &DfCD::getSuperMatrixElements);
            
            this->saveLxc(Lxc);
            this->destroyEngines();
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
    this->log_.info(TlUtils::format("save I2PQ database: %s", filepath.c_str()));
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

void DfCD::saveLjk(const TlRowVectorMatrix& Ljk)
{
    this->log_.info("save Ljk");
    const std::string path = DfObject::getLjkMatrixPath();
    Ljk.saveByTlColVectorMatrix(path);
}

void DfCD::saveLk(const TlRowVectorMatrix& Lk)
{
    this->log_.info("save Lk");
    const std::string path = DfObject::getLkMatrixPath();
    Lk.saveByTlColVectorMatrix(path);
}

void DfCD::saveLxc(const TlRowVectorMatrix& Lxc)
{
    this->log_.info("save Lxc");
    const std::string path = this->getLxcMatrixPath();
    Lxc.saveByTlColVectorMatrix(path);
}

void DfCD::debugOutLjk(const TlMatrix& Ljk)
{
    const std::string path = TlUtils::format("%s.debug", DfObject::getLjkMatrixPath().c_str());
    Ljk.save(path);
}

void DfCD::debugOutLk(const TlMatrix& Lk)
{
    const std::string path = TlUtils::format("%s.debug", DfObject::getLkMatrixPath().c_str());
    Lk.save(path);
}

void DfCD::debugOutLxc(const TlMatrix& Lxc)
{
    const std::string path = TlUtils::format("%s.debug", DfObject::getLxcMatrixPath().c_str());
    Lxc.save(path);
}

TlColVectorMatrix DfCD::getLjk()
{
    TlColVectorMatrix Ljk = DfObject::getLjkMatrix<TlColVectorMatrix>();
    return Ljk;
}

TlColVectorMatrix DfCD::getLk()
{
    TlColVectorMatrix Lk = DfObject::getLkMatrix<TlColVectorMatrix>();
    return Lk;
}

TlColVectorMatrix DfCD::getLxc()
{
    TlColVectorMatrix Lxc = DfObject::getLxcMatrix<TlColVectorMatrix>();
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
    this->log_.info("calc J by CD method (serial).");

    // cholesky vector
    const TlColVectorMatrix L = this->getLjk();
    this->log_.info(TlUtils::format("L(J): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();
    
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlVector vP = this->getScreenedDensityMatrix(I2PQ);

    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);

    const index_type numOfI = I2PQ.size();
    TlVector vJ(numOfI);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        const TlVector LI = L.getColVector(I);
        assert(LI.getSize() == vJ.getSize());

        TlVector tmpLI = LI;
        const double qi = tmpLI.dot(vP).sum();

        vJ += qi*LI;
    }

    this->expandJMatrix(vJ, I2PQ, pJ);
    this->finalize(pJ);
}

TlVector DfCD::getScreenedDensityMatrix(const PQ_PairArray& I2PQ)
{
    const TlSymmetricMatrix P = this->getPMatrix();
    const std::size_t numOfI = I2PQ.size();
    TlVector answer(numOfI);
    
    for (std::size_t i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        const double coef = (r != c) ? 2.0 : 1.0;
        answer.set(i, coef * P.get(r, c));
    }

    return answer;
}

TlVector DfCD::getScreenedDensityMatrix(const RUN_TYPE runType, const PQ_PairArray& I2PR)
{
    TlSymmetricMatrix P;
    switch (runType) {
    case RUN_RKS:
        P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
        break;

    case RUN_UKS_ALPHA:
    case RUN_UKS_BETA:
        P = this->getPpqMatrix<TlSymmetricMatrix>(runType, this->m_nIteration -1);
        break;

    case RUN_ROKS_ALPHA:
        {
            P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration -1);
            P += this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_OPEN, this->m_nIteration -1);
        }
        break;

    case RUN_ROKS_BETA:
        {
            P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration -1);
        }
        break;

    default:
        this->log_.critical(TlUtils::format("Program Error: %s:%d", __FILE__, __LINE__));
        CnErr.abort();
    }

    const std::size_t numOfI = I2PR.size();
    TlVector answer(numOfI);
    
    for (std::size_t i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PR[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        // const double coef = (r != c) ? 2.0 : 1.0;
        // const double coef = (r != c) ? 1.0 : 0.5;
        const double coef = 1.0;
        answer.set(i, coef * P.get(r, c));
    }

    return answer;
}


void DfCD::expandJMatrix(const TlVector& vJ, const PQ_PairArray& I2PQ, TlSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    const index_type numOfI = I2PQ.size();
    for (index_type i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        pJ->set(r, c, vJ.get(i));
    }
}


void DfCD::expandKMatrix(const TlVector& vK, const PQ_PairArray& I2PR, TlSymmetricMatrix* pK)
{
    assert(pK != NULL);
    const index_type numOfI = I2PR.size();
    for (index_type i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PR[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        const double coef = (r != c) ? 0.5 : 1.0;
        pK->add(r, c, coef * vK.get(i));
    }
}


void DfCD::getJ_A(TlSymmetricMatrix* pJ)
{
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    const TlColVectorMatrix L = this->getLjk();
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
    switch(this->fastCDK_mode_) {
    case FASTCDK_NONE:
        this->getK_S_woCD(runType, pK);
        break;
        
    case FASTCDK_DEBUG_FULL_SUPERMATRIX:
    case FASTCDK_DEBUG_SUPERMATRIX:
    case FASTCDK_PRODUCTIVE_FULL:
    case FASTCDK_PRODUCTIVE:
        this->getK_S_fast(runType, pK);
        break;

    default:
        this->log_.critical(TlUtils::format("%s: %d: program error.", __FILE__, __LINE__));
        CnErr.abort();
    }
}


TlSymmetricMatrix DfCD::getPMatrix(const RUN_TYPE runType, int itr)
{
    TlSymmetricMatrix P;
    {
        this->log_.info("use density matrix.");
        switch (runType) {
        case RUN_RKS:
            P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, itr);
            break;
            
        case RUN_UKS_ALPHA:
        case RUN_UKS_BETA:
            P = this->getPpqMatrix<TlSymmetricMatrix>(runType, itr);
            break;
            
        case RUN_ROKS_ALPHA:
            {
                P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSED, itr);
                P += this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_OPEN, itr);
            }
            break;
            
        case RUN_ROKS_BETA:
            {
                P = 0.5 * this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSED, itr);
            }
            break;
            
        default:
            this->log_.critical(TlUtils::format("Program Error: %s:%d", __FILE__, __LINE__));
            CnErr.abort();
        }
        this->log_.info("CD: density matrix");
    }

    return P;
}

void DfCD::getK_S_woCD(const RUN_TYPE runType,
                       TlSymmetricMatrix *pK)
{
    this->log_.info("calc K by CD method (serial).");
    //const index_type numOfAOs = this->m_nNumOfAOs;

    const TlColVectorMatrix L = this->getLjk();
    this->log_.info(TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const TlSymmetricMatrix P = this->getPMatrix(runType, this->m_nIteration -1);

    this->log_.info("start loop");
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        const TlSymmetricMatrix l = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(l.getNumOfRows() == this->m_nNumOfAOs);

        TlMatrix X = l * P;
        X *= l;

        *pK += X;
    }
    
    *pK *= -1.0;
    this->log_.info("finalize");
    this->finalize(pK);
}


void DfCD::getK_S_fast(const RUN_TYPE runType,
                       TlSymmetricMatrix* pK)
{
    this->log_.info("calc K(fast) by CD method (serial).");

    // cholesky vector
    const TlColVectorMatrix L = this->getLk();
    this->log_.info(TlUtils::format("L(J): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();
    
    const PQ_PairArray I2PR = this->getI2PQ(this->getI2prVtrPath());
    const TlVector vP = this->getScreenedDensityMatrix(runType, I2PR);

    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);

    const index_type numOfI = I2PR.size();
    TlVector vK(numOfI);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        const TlVector LI = L.getColVector(I);
        assert(LI.getSize() == vK.getSize());

        TlVector tmpLI = LI;
        const double qi = tmpLI.dot(vP).sum();

        vK += qi*LI;
    }
    vK *= -1.0;

    this->expandKMatrix(vK, I2PR, pK);
    this->finalize(pK);
}


void DfCD::getK_A(const RUN_TYPE runType,
                  TlSymmetricMatrix *pK)
{
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set"]);

    const TlColVectorMatrix L = this->getLjk();
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
                                (*this->pPdfParam_)["basis_set"]);
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();

    pM->resize(numOfAOs);

    // cholesky vector
    const TlColVectorMatrix L = this->getLxc();
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
                                  (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"],
                                  (*this->pPdfParam_)["basis_set_gridfree"]);
    // const index_type numOfAOs = orbInfo_p.getNumOfOrbitals();
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    const TlColVectorMatrix L = this->getLxc();
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
        assert(l.getNumOfRows() == orbInfo_p.getNumOfOrbitals());
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
        P  = this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration -1);
        P += this->getPpqMatrix<TlSymmetricMatrix>(RUN_ROKS_OPEN, this->m_nIteration -1);
        break;
        
    default:
        this->log_.critical("program error");
        break;
    }
    return P;
}


TlRowVectorMatrix DfCD::calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
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

    TlRowVectorMatrix L(N, 1, 1, 0, this->isEnableMmap_);
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
    L.reserveColSize(division);
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
            L.reserveColSize(division * progress);
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
        const TlVector L_pm = L.getVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();

            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp critical(DfCD__calcCholeskyVectorsOnTheFlyA)
            {
                L_xm[i] += l_m_pi;
                d[pivot_i] -= l_m_pi * l_m_pi;
            }
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
                         TlVector *pDiagonals)
{
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs + 1) / 2;

    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffThreshold_primitive_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(numOfPQs);
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


void DfCD::calcDiagonals_K_full(const TlOrbitalInfoObject& orbInfo,
                                PQ_PairArray *pI2PR,
                                TlVector *pDiagonals)
{
    const double tau = this->CDAM_tau_K_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffThreshold_primitive_));

    // initialize
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPRs = numOfAOs * numOfAOs;
    pI2PR->clear();
    pI2PR->reserve(numOfPRs);
    TlSparseMatrix diagonalMat(numOfAOs, numOfAOs);

    // task
    this->log_.info("diagonal calculation(K): start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                          orbInfo,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_K_full_kernel(orbInfo,
                                          taskList,
                                          &diagonalMat, pI2PR);
        hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                         orbInfo,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    this->finalize_I2PQ(pI2PR);
    this->finalize(&diagonalMat);

    // set diagonals
    const index_type numOfI = pI2PR->size();
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PR)[i].index1();
        const index_type col = (*pI2PR)[i].index2();
        const double value = diagonalMat.get(row, col);
        (*pDiagonals)[i] = value;
    }
}

void DfCD::calcDiagonals_K_half(const TlOrbitalInfoObject& orbInfo,
                                PQ_PairArray *pI2PR,
                                TlVector *pDiagonals)
{
    const double tau = this->CDAM_tau_K_;
    this->log_.info(TlUtils::format("CDAM tau(K): %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffThreshold_primitive_));

    // initialize
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPRs = numOfAOs * (numOfAOs +1) / 2;
    pI2PR->clear();
    pI2PR->reserve(numOfPRs);
    TlSparseSymmetricMatrix diagonalMat(numOfAOs);

    // task
    this->log_.info("diagonal calculation(K): start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_K_half_kernel(orbInfo,
                                          taskList,
                                          &diagonalMat, pI2PR);
        hasTask = pDfTaskCtrl->getQueue2(orbInfo,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    this->finalize_I2PQ(pI2PR);
    this->finalize(&diagonalMat);

    // set diagonals
    const index_type numOfI = pI2PR->size();
    this->log_.info(TlUtils::format("num of diagonals: %d", numOfI));
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PR)[i].index1();
        const index_type col = (*pI2PR)[i].index2();
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
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffThreshold_primitive_));

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
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        assert(0 <= threadID);
        assert(threadID < this->numOfThreads_);
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
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
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
        }
#endif // _OPENMP
    }
}


void DfCD::calcDiagonals_K_full_kernel(const TlOrbitalInfoObject& orbInfo,
                                       const std::vector<DfTaskCtrl::Task2>& taskList,
                                       TlSparseMatrix *pDiagonalMat,
                                       PQ_PairArray *pI2PR)
{
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs, numOfAOs);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PR;
        TlSparseMatrix local_diagMat(numOfAOs, numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        assert(0 <= threadID);
        assert(threadID < this->numOfThreads_);
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexR = taskList[i].shellIndex2;
            const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;

            assert(this->pEngines_[threadID] != NULL);
            this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                            0, orbInfo, shellIndexP,
                                            0, orbInfo, shellIndexR,
                                            0, orbInfo, shellIndexR);
                
            const int maxStepsRR = maxStepsR * maxStepsR;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                const int pp_index = p * maxStepsP + p;
                for (int r = 0; r < maxStepsR; ++r) {
                    const index_type indexR = shellIndexR + r;
                    const index_type rr_index = r * maxStepsR + r;
                    
                    const int pprr_index = pp_index * maxStepsRR + rr_index;
                        
                    const double value = this->pEngines_[threadID]->value(pprr_index);

                    // for schwartz
                    maxValue = std::max(maxValue, std::fabs(value));
                        
                    // for I~ to pq table
                    if (std::fabs(value) > tau) {
                        if (value > 0) {
                            local_diagMat.set(indexP, indexR, value);
                            local_I2PR.push_back(Index2(indexP, indexR));
                        } else {
                            this->log_.warn(TlUtils::format("pprr_value: (%d %d)=% e is spoiled.",
                                                            indexP, indexR, value));
                        }
                    }
                }
            }
        }

        // add up
#ifdef _OPENMP
        {
#pragma omp critical(DfCD__calcDiagonals_kernel_1)
            {
                pI2PR->insert(pI2PR->end(),
                              local_I2PR.begin(), local_I2PR.end());
            }
#pragma omp critical(DfCD__calcDiagonals_kernel_2)
            {
                pDiagonalMat->merge(local_diagMat);
            }
        }
#else
        {
            *pI2PR = local_I2PR;
            *pDiagonalMat = local_diagMat;
        }
#endif // _OPENMP
    }
}

void DfCD::calcDiagonals_K_half_kernel(const TlOrbitalInfoObject& orbInfo,
                                       const std::vector<DfTaskCtrl::Task2>& taskList,
                                       TlSparseSymmetricMatrix *pDiagonalMat,
                                       PQ_PairArray *pI2PR)
{
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        PQ_PairArray local_I2PR;
        TlSparseSymmetricMatrix local_diagMat(numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        assert(0 <= threadID);
        assert(threadID < this->numOfThreads_);
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexR = taskList[i].shellIndex2;
            const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;

            double maxValue = 0.0;
            {
                assert(this->pEngines_[threadID] != NULL);
                this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexR,
                                                0, orbInfo, shellIndexR);
                
                const int maxStepsRR = maxStepsR * maxStepsR;
                double tmpMaxValue = 0.0;
                for (int p = 0; p < maxStepsP; ++p) {
                    const index_type indexP = shellIndexP + p;
                    const index_type pp_index = p * maxStepsP + p;
                    for (int r = 0; r < maxStepsR; ++r) {
                        const index_type indexR = shellIndexR + r;
                        const index_type rr_index = r * maxStepsR + r;
                        
                        if (indexP >= indexR) {
                            const int pprr_index = pp_index * maxStepsRR + rr_index;
                            
                            const double coef = (indexP != indexR) ? 2.0 : 1.0;
                            const double value = coef * this->pEngines_[threadID]->value(pprr_index);
                            
                            
                            // for I~ to pq table
                            local_diagMat.add(indexP, indexR, value);
                        }
                    }
                }

                maxValue += tmpMaxValue;
            }

            {
                assert(this->pEngines_[threadID] != NULL);
                this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexR,
                                                0, orbInfo, shellIndexP,
                                                0, orbInfo, shellIndexR);
                
                double tmpMaxValue = 0.0;
                for (int p = 0; p < maxStepsP; ++p) {
                    const index_type indexP = shellIndexP + p;
                    for (int r = 0; r < maxStepsR; ++r) {
                        const index_type indexR = shellIndexR + r;
                        
                        if (indexP > indexR) { // indexP != indexR
                            const int prpr_index = ((p * maxStepsR + r) * maxStepsP + p) * maxStepsR + r;
                            
                            const double coef = (indexP != indexR) ? 2.0 : 1.0;
                            const double value = coef * this->pEngines_[threadID]->value(prpr_index);
                            
                            
                            // for I~ to pq table
                            local_diagMat.add(indexP, indexR, value);
                        }
                    }
                }

                maxValue += tmpMaxValue;
            }

            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int r = 0; r < maxStepsR; ++r) {
                    const index_type indexR = shellIndexR + r;

                    if (indexP >= indexR) {
                        if (std::fabs(local_diagMat.get(indexP, indexR)) > tau) {
                            local_I2PR.push_back(Index2(indexP, indexR));
                        }
                    }
                }
            }
        }

        // add up
#ifdef _OPENMP
        {
#pragma omp critical(DfCD__calcDiagonals_kernel_1)
            {
                pI2PR->insert(pI2PR->end(),
                              local_I2PR.begin(), local_I2PR.end());
            }
#pragma omp critical(DfCD__calcDiagonals_kernel_2)
            {
                pDiagonalMat->merge(local_diagMat);
            }
        }
#else
        {
            *pI2PR = local_I2PR;
            *pDiagonalMat = local_diagMat;
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
        assert(0 <= threadID);
        assert(threadID < this->numOfThreads_);
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbInfo_p.getShellType(shellIndexP);
            const int shellTypeQ = orbInfo_q.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            assert(this->pEngines_[threadID] != NULL);
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
                             const PQ_PairArray& I2PQ)
{
    this->ERI_cache_.clear();

    const int start = 0;
    const int end = G_col_list.size();

    const std::vector<IndexPair4S> calcList = this->getCalcList(orbInfo, G_row,
                                                                G_col_list, start, end, I2PQ);
    this->calcERIs(orbInfo, calcList);
    const std::vector<double> answer = this->setERIs(orbInfo, G_row,
                                                     G_col_list, start, end, I2PQ);

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


std::vector<DfCD::IndexPair4S> 
DfCD::getCalcList(const TlOrbitalInfoObject& orbInfo,
                  const index_type G_row,
                  const std::vector<index_type>& G_col_list,
                  const index_type start,
                  const index_type end,
                  const PQ_PairArray& I2PQ)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    std::set<IndexPair4S> calcSet;

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();
    const index_type shellIndexP = orbInfo.getShellIndex(indexP);
    const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);

#pragma omp parallel for schedule(runtime)
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();
        const index_type shellIndexR = orbInfo.getShellIndex(indexR);
        const index_type shellIndexS = orbInfo.getShellIndex(indexS);

        IndexPair4S index4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList) 
        {
            calcSet.insert(index4);
        }
    }

    std::vector<IndexPair4S> calcList(calcSet.size());
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
                    const std::vector<IndexPair4S>& calcList) 
{
    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
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
            
            {
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
    assert(orbInfo_p.getMaxShellType() == orbInfo_q.getMaxShellType());

    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType_A local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
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
              const index_type start,
              const index_type end,
              const PQ_PairArray& I2PQ)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);
#pragma omp parallel for
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();

        double value = 0.0;
        if (this->getCachedValue(orbInfo, indexP, indexQ, indexR, indexS, this->ERI_cache_, &value)) { 
            answer[i] = value;
        } else {
            CnErr.abort(TlUtils::format("%s: %d: not found value in cache.", __FILE__, __LINE__));
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
            
#pragma omp critical(DfCD__setERIs_set_answer)
            {
                answer[i] += values.at(index);
            }
        }
    }

    return answer;
}

////////////////////////////////////////////////////////////////////////////////
// for DEBUG 
// 
////////////////////////////////////////////////////////////////////////////////
TlSymmetricMatrix DfCD::getSuperMatrix(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PQ) 
{
    // const index_type numOfOrbs = orbInfo.getNumOfOrbitals();

    PQ_PairArray I2PQ;
    TlVector d; // 対角成分
    this->calcDiagonals(orbInfo,
                        &I2PQ, &d);
    if (pI2PQ != NULL) {
        *pI2PQ = I2PQ;
    }
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

TlSymmetricMatrix DfCD::getSuperMatrix_K_full(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PR) 
{
    const index_type numOfOrbs = orbInfo.getNumOfOrbitals();

    PQ_PairArray I2PR;
    TlVector d; // 対角成分
    this->calcDiagonals_K_full(orbInfo,
                               &I2PR, &d);
    if (pI2PR != NULL) {
        *pI2PR = I2PR;
    }
    const std::size_t N = I2PR.size();

    TlSymmetricMatrix V(N);

    const int v2size = numOfOrbs * (numOfOrbs +1) / 2;
    TlSymmetricMatrix V2(v2size);
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t i = 0; i < N; ++i) {
            const index_type indexP = I2PR[i].index1();
            const index_type indexR = I2PR[i].index2();
            const index_type shellIndexP = orbInfo.getShellIndex(indexP);
            const index_type shellIndexR = orbInfo.getShellIndex(indexR);
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeR = indexR - shellIndexR;
            // const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            //const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            
            // for (std::size_t j = 0; j <= i; ++j) {
            for (std::size_t j = 0; j < N; ++j) {
                const index_type indexQ = I2PR[j].index1();
                const index_type indexS = I2PR[j].index2();
                const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
                const index_type shellIndexS = orbInfo.getShellIndex(indexS);
                const int basisTypeQ = indexQ - shellIndexQ;
                const int basisTypeS = indexS - shellIndexS;
                const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
                const int shellTypeS = orbInfo.getShellType(shellIndexS);
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsS = 2 * shellTypeS + 1;

                //
                {
                    this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexQ,
                                                    0, orbInfo, shellIndexR,
                                                    0, orbInfo, shellIndexS);
                    const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
                    const double value = this->pEngines_[threadID]->value(index);
                    // V.add(i, j, value);
                    V.set(i, j, value);

                    {
                        index_type p = indexP;
                        index_type q = indexQ;
                        index_type r = indexR;
                        index_type s = indexS;
                        
                        if (p < r) {
                            std::swap(p, r);
                        }
                        index_type pr = p * (p +1) /2 + r;
                        if (q < s) {
                            std::swap(q, s);
                        }
                        index_type qs = q * (q +1) /2 + s;

                        assert(pr < v2size);
                        assert(qs < v2size);
                        if (pr >= qs) {
                            V2.add(pr, qs, value);
                        }
                    }
                }
            }
        }
    }

    V2.save("V2.mat");
    return V;
}

TlSymmetricMatrix DfCD::getSuperMatrix_K_half(const TlOrbitalInfoObject& orbInfo, PQ_PairArray* pI2PR) 
{
    // const index_type numOfOrbs = orbInfo.getNumOfOrbitals();

    PQ_PairArray I2PR;
    TlVector d; // 対角成分
    this->calcDiagonals_K_half(orbInfo,
                               &I2PR,
                               &d);
    if (pI2PR != NULL) {
        *pI2PR = I2PR;
    }
    const std::size_t N = I2PR.size();

    TlSymmetricMatrix V(N);
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t i = 0; i < N; ++i) {
            const index_type indexP = I2PR[i].index1();
            const index_type indexR = I2PR[i].index2();
            const index_type shellIndexP = orbInfo.getShellIndex(indexP);
            const index_type shellIndexR = orbInfo.getShellIndex(indexR);
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeR = indexR - shellIndexR;
            const int shellTypeP = orbInfo.getShellType(shellIndexP);
            const int shellTypeR = orbInfo.getShellType(shellIndexR);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;

            // 非対角項
            for (std::size_t j = 0; j < i; ++j) {
                const index_type indexQ = I2PR[j].index1();
                const index_type indexS = I2PR[j].index2();
                const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
                const index_type shellIndexS = orbInfo.getShellIndex(indexS);
                const int basisTypeQ = indexQ - shellIndexQ;
                const int basisTypeS = indexS - shellIndexS;
                const int shellTypeQ = orbInfo.getShellType(shellIndexQ);
                const int shellTypeS = orbInfo.getShellType(shellIndexS);
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsS = 2 * shellTypeS + 1;

                {
                    this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexQ,
                                                    0, orbInfo, shellIndexR,
                                                    0, orbInfo, shellIndexS);
                    const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;

                    const double coef = (indexP != indexR) ? 2.0 : 1.0;
                    const double value = coef * this->pEngines_[threadID]->value(index);
                    V.add(i, j, value);
                }

                if (indexQ != indexS) {
                    this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexS,
                                                    0, orbInfo, shellIndexR,
                                                    0, orbInfo, shellIndexQ);
                    const int index = ((basisTypeP * maxStepsS + basisTypeS) * maxStepsR + basisTypeR) * maxStepsQ + basisTypeQ;
                    
                    const double coef = (indexP != indexR) ? 2.0 : 1.0;
                    const double value = coef * this->pEngines_[threadID]->value(index);
                    V.add(i, j, value);
                }
            }

            // 対角項
            {
                {
                    this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexR,
                                                    0, orbInfo, shellIndexR);
                    const int index = ((basisTypeP * maxStepsP + basisTypeP) * maxStepsR + basisTypeR) * maxStepsR + basisTypeR;

                    double coef = 1.0;
                    coef *= (indexP != indexR) ? 2.0 : 1.0;
                    const double value = coef * this->pEngines_[threadID]->value(index);

                    V.add(i, i, value);
                }

                if (indexP != indexR) {
                    this->pEngines_[threadID]->calc(0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexR,
                                                    0, orbInfo, shellIndexP,
                                                    0, orbInfo, shellIndexR);
                    const int index = ((basisTypeP * maxStepsR + basisTypeR) * maxStepsP + basisTypeP) * maxStepsR + basisTypeR;

                    double coef = 1.0;
                    coef *= (indexP != indexR) ? 2.0 : 1.0;
                    const double value = coef * this->pEngines_[threadID]->value(index);

                    V.add(i, i, value);
                }

            }
        }
    }

    return V;
}


TlSymmetricMatrix DfCD::getSuperMatrix(const TlOrbitalInfoObject& orbInfo_p,
                                       const TlOrbitalInfoObject& orbInfo_q,
                                       PQ_PairArray* pI2PQ) 
{
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();

    PQ_PairArray I2PQ;
    TlSparseMatrix schwartzTable(numOfOrbs_p, numOfOrbs_q);
    TlVector d; // 対角成分
    this->calcDiagonalsA(orbInfo_p, orbInfo_q,
                         &I2PQ, &schwartzTable, &d);
    if (pI2PQ != NULL) {
        *pI2PQ = I2PQ;
    }
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

#pragma omp critical(DfCD__calcCholeskyVectors)
            {
                L_xm[i] += l_m_pi; // for OpenMP
                d[pivot_i] -= l_m_pi * l_m_pi;
            }
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

// ---------------------------------------------------------------------------
// void DfCD::calcDiagonals_J(const TlOrbitalInfoObject& orbInfo,
//                            PQ_PairArray* pI2PQ,
//                            SparseSymmetricMatrix* pSchwartzTable,
//                            TlVector* pDiagonals)
// {
// }

// void DfCD::calcDiagonals_K(const TlOrbitalInfoObject& orbInfo,
//                            PQ_PairArray* pI2PQ,
//                            SparseSymmetricMatrix* pSchwartzTable,
//                            TlVector* pDiagonals)
// {
// }

// K full --------------------------------------------------------------
std::vector<double>
DfCD::getSuperMatrixElements_K_full(const TlOrbitalInfoObject& orbInfo,
                                    const index_type G_row,
                                    const std::vector<index_type>& G_col_list,
                                    const PQ_PairArray& I2PQ)
{
    this->ERI_cache_.clear();

    const int start = 0;
    const int end = G_col_list.size();
    
    const std::vector<IndexPair4S> calcList = this->getCalcList_K_full(orbInfo, G_row,
                                                                       G_col_list, start, end,
                                                                       I2PQ);
    this->calcERIs_K(orbInfo, calcList);
    const std::vector<double> answer = this->setERIs_K_full(orbInfo, G_row,
                                                            G_col_list, start, end, I2PQ);

    return answer;
}

std::vector<DfCD::IndexPair4S> 
DfCD::getCalcList_K_full(const TlOrbitalInfoObject& orbInfo,
                         const index_type G_row,
                         const std::vector<index_type>& G_col_list,
                         const index_type start,
                         const index_type end,
                         const PQ_PairArray& I2PR)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    std::set<IndexPair4S> calcSet;

    const index_type indexP = I2PR[G_row].index1();
    const index_type indexR = I2PR[G_row].index2();
    const index_type shellIndexP = orbInfo.getShellIndex(indexP);
    const index_type shellIndexR = orbInfo.getShellIndex(indexR);

#pragma omp parallel for schedule(runtime)
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexQ = I2PR[G_col].index1();
        const index_type indexS = I2PR[G_col].index2();
        const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
        const index_type shellIndexS = orbInfo.getShellIndex(indexS);

        IndexPair4S index4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList_K_full) 
        {
            calcSet.insert(index4);
        }
    }

    std::vector<IndexPair4S> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}


std::vector<double>
DfCD::setERIs_K_full(const TlOrbitalInfoObject& orbInfo,
                     const index_type G_row,
                     const std::vector<index_type> G_col_list,
                     const index_type start,
                     const index_type end,
                     const PQ_PairArray& I2PR)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    const index_type indexP = I2PR[G_row].index1();
    const index_type indexR = I2PR[G_row].index2();

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);

#pragma omp parallel for
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];
        index_type indexQ = I2PR[G_col].index1();
        index_type indexS = I2PR[G_col].index2();

        double value = 0.0;
        if (this->getCachedValue(orbInfo, indexP, indexQ, indexR, indexS, this->ERI_cache_, &value)) {
            answer[i] = value;
        } else {
            CnErr.abort(TlUtils::format("%s: %d: not found value in cache.", __FILE__, __LINE__));
        }
    }

    return answer;
}


// K half --------------------------------------------------------------
std::vector<double>
DfCD::getSuperMatrixElements_K_half(const TlOrbitalInfoObject& orbInfo,
                                    const index_type G_row,
                                    const std::vector<index_type>& G_col_list,
                                    const PQ_PairArray& I2PR)
{
    this->ERI_cache_.clear();

    const int start = 0;
    const int end = G_col_list.size();

    const std::vector<IndexPair4S> calcList = this->getCalcList_K_half(orbInfo, G_row,
                                                                       G_col_list, start, end, I2PR);
    this->calcERIs_K(orbInfo, calcList);
    const std::vector<double> answer = this->setERIs_K_half(orbInfo, G_row,
                                                            G_col_list, start, end, I2PR);

    return answer;
}


std::vector<DfCD::IndexPair4S> 
DfCD::getCalcList_K_half(const TlOrbitalInfoObject& orbInfo,
                         const index_type G_row,
                         const std::vector<index_type>& G_col_list,
                         const index_type start,
                         const index_type end,
                         const PQ_PairArray& I2PR)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    std::set<IndexPair4S> calcSet;

    const index_type indexP = I2PR[G_row].index1();
    const index_type indexR = I2PR[G_row].index2();
    const index_type shellIndexP = orbInfo.getShellIndex(indexP);
    const index_type shellIndexR = orbInfo.getShellIndex(indexR);

#pragma omp parallel for schedule(runtime)
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexQ = I2PR[G_col].index1();
        const index_type indexS = I2PR[G_col].index2();
        const index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
        const index_type shellIndexS = orbInfo.getShellIndex(indexS);

        if ((indexP + indexR == indexQ + indexS) &&
            (indexP * indexR == indexQ * indexS)) {
            // 対角項
            IndexPair4S index4(shellIndexP, shellIndexP, shellIndexR, shellIndexR);
#pragma omp critical(DfCD__getCalcList_K_half) 
            {
                calcSet.insert(index4);
            }
            if (indexP != indexR) {
                IndexPair4S index4(shellIndexP, shellIndexR, shellIndexP, shellIndexR);
#pragma omp critical(DfCD__getCalcList_K_half)
                {
                    calcSet.insert(index4);
                }
            }
        } else {
            // 非対角項
            {
                IndexPair4S index4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList_K_half) 
                {
                    calcSet.insert(index4);
                }
            }
            if (indexQ != indexS) {
                IndexPair4S index4(shellIndexP, shellIndexS, shellIndexR, shellIndexQ);
#pragma omp critical(DfCD__getCalcList_K_half) 
                {
                    calcSet.insert(index4);
                }
            }
        }
    }

    std::vector<IndexPair4S> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}


void DfCD::calcERIs_K(const TlOrbitalInfoObject& orbInfo,
                      const std::vector<IndexPair4S>& calcList)
{
    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEngines_[threadID]->setPrimitiveLevelThreshold(this->cutoffThreshold_primitive_);
        
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
            
            {
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

std::vector<double>
DfCD::setERIs_K_half(const TlOrbitalInfoObject& orbInfo,
                     const index_type G_row,
                     const std::vector<index_type> G_col_list,
                     const index_type start,
                     const index_type end,
                     const PQ_PairArray& I2PQ)
{
    assert(0 <= start);
    assert(end <= static_cast<index_type>(G_col_list.size()));
    assert(start <= end);

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexR = I2PQ[G_row].index2();

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);

#pragma omp parallel for
    for (index_type i = start; i < end; ++i) {
        const index_type G_col = G_col_list[i];
        index_type indexQ = I2PQ[G_col].index1();
        index_type indexS = I2PQ[G_col].index2();

        double value = 0.0;
        double tmp = 0.0;
        if (G_row == G_col) {
            // 対角項
            {
                if (this->getCachedValue(orbInfo, indexP, indexP, indexR, indexR, this->ERI_cache_, &tmp)) {
                    const double coef = (indexP != indexR) ? 2.0 : 1.0;
#pragma omp critical(DfCD__setERIs_set_answer)
                    {
                        value += coef * tmp;
                    }
                } else {
                    CnErr.abort("not found value in cache1.");
                }
            }
            if (indexP != indexR) {
                if (this->getCachedValue(orbInfo, indexP, indexR, indexP, indexR, this->ERI_cache_, &tmp)) {
                    // const double coef = (indexP != indexR) ? 2.0 : 1.0;
                    const double coef = 2.0;
#pragma omp critical(DfCD__setERIs_set_answer)
                    {
                        value += coef * tmp;
                    }
                } else {
                    CnErr.abort("not found value in cache2.");
                }
            }
        } else {
            // 非対角項
            if (this->getCachedValue(orbInfo, indexP, indexQ, indexR, indexS, this->ERI_cache_, &tmp)) {
                const double coef = (indexP != indexR) ? 2.0 : 1.0;
#pragma omp critical(DfCD__setERIs_set_answer)
                {
                    value += coef * tmp;
                }
            } else {
                CnErr.abort("not found value in cache3.");
            }
            
            if (indexQ != indexS) {
                if (this->getCachedValue(orbInfo, indexP, indexS, indexR, indexQ, this->ERI_cache_, &tmp)) {
                    const double coef = (indexP != indexR) ? 2.0 : 1.0;
#pragma omp critical(DfCD__setERIs_set_answer)
                    {
                        value += coef * tmp;
                    }
                } else {
                    CnErr.abort("not found value in cache4.");
                }
            }
        }

        answer[i] = value;

        // debug
        // {
        //     index_type I = 0;
        //     assert(this->get_I_index(this->debug_I2PQ_, indexP, indexR, &I));
        //     index_type J = 0;
        //     assert(this->get_I_index(this->debug_I2PQ_, indexQ, indexS, &J));

        //     const double check = this->debug_V_.get(I, J);
        //     if (std::fabs(check - value) > 1.0E-5) {
        //         this->log_.info(TlUtils::format("(%d, %d|%d, %d) I=%d, J=%d: c=% e, v= %e",
        //                                         indexP, indexQ, indexR, indexS, I, J,
        //                                         check, value));
        //         CnErr.abort();
        //     }
        //     // else {
        //     //     this->log_.info(TlUtils::format("(%d, %d|%d, %d) I=%d, J=%d: c=% e, v= %e",
        //     //                                     indexP, indexQ, indexR, indexS, I, J,
        //     //                                     check, value));
        //     // }
        // }

    }

    return answer;
}

bool DfCD::getCachedValue(const TlOrbitalInfoObject& orbInfo,
                          index_type indexP,
                          index_type indexQ,
                          index_type indexR,
                          index_type indexS,
                          const ERI_CacheType& cache,
                          double* pValue)
{
    assert(pValue != NULL);
    bool answer = false;
    
    index_type shellIndexP = orbInfo.getShellIndex(indexP);
    index_type shellIndexQ = orbInfo.getShellIndex(indexQ);
    index_type shellIndexR = orbInfo.getShellIndex(indexR);
    index_type shellIndexS = orbInfo.getShellIndex(indexS);

    if (shellIndexP < shellIndexQ) {
        std::swap(shellIndexP, shellIndexQ);
        std::swap(indexP, indexQ);
    }
    if (shellIndexR < shellIndexS) {
        std::swap(shellIndexR, shellIndexS);
        std::swap(indexR, indexS);
    }
    if (Index2(shellIndexP, shellIndexQ) < Index2(shellIndexR, shellIndexS)) {
        std::swap(shellIndexP, shellIndexR);
        std::swap(shellIndexQ, shellIndexS);
        std::swap(indexP, indexR);
        std::swap(indexQ, indexS);
    }

    ERI_CacheType::const_iterator it = cache.find(IndexPair4S(shellIndexP, shellIndexQ, shellIndexR, shellIndexS));

    if (it != cache.end()) {
        std::vector<double> values = it->second;
        
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

        *pValue = values.at(index);
        answer = true;
    }
    //  else {
    //     this->log_.critical(TlUtils::format("not found: (%d %d|%d %d)", shellIndexP, shellIndexQ, shellIndexR, shellIndexS));
    // }
    
    return answer;
}


TlRowVectorMatrix DfCD::calcCholeskyVectorsOnTheFlyS_new(const TlOrbitalInfoObject& orbInfo,
                                                         const std::string& I2PQ_path,
                                                         const double threshold,
                                                         void(DfCD::*calcDiagonalsFunc)(
                                                             const TlOrbitalInfoObject&,
                                                             PQ_PairArray*,
                                                             TlVector*),
                                                         std::vector<double>(DfCD::*getSuperMatrixElements)(
                                                             const TlOrbitalInfoObject&,
                                                             const index_type,
                                                             const std::vector<index_type>&,
                                                             const PQ_PairArray&))
{
    this->log_.info("call on-the-fly Cholesky Decomposition routine (symmetric)");
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
    TlRowVectorMatrix L(numOfPQtilde, 1, 1, 0, this->isEnableMmap_);

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
        std::vector<double> G_pm;
        const index_type numOf_G_cols = numOfPQtilde -(numOfCDVcts +1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[(numOfCDVcts +1) +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            G_pm = (this->*getSuperMatrixElements)(orbInfo,
                                                   pivot_m, G_col_list, I2PQ);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        const TlVector L_pm = L.getVector(pivot_m);
        assert(L_pm.getSize() == numOfCDVcts +1);

        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[(numOfCDVcts +1) +i]; // from (m+1) to N
            TlVector L_pi = L.getVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();
            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp atomic
            L_xm[i] += l_m_pi;
#pragma omp atomic
            diagonals[pivot_i] -= l_m_pi * l_m_pi;
        }

        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[(numOfCDVcts +1) +i]; // from (m+1) to N
            L.set(pivot_i, numOfCDVcts, L_xm[i]);
        }

        error = diagonals[pivot[numOfCDVcts]];
        ++numOfCDVcts;
    }
    L.resize(numOfPQtilde, numOfCDVcts);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", numOfCDVcts));

    return L;
}

bool DfCD::get_I_index(const PQ_PairArray& I2PQ, const index_type p, index_type q,
                       index_type* pI)
{
    assert(pI != NULL);

    bool answer = false;
    const Index2 pq(p, q);
    const int numOfIs = I2PQ.size();
    for (int i = 0; i < numOfIs; ++i) {
        if (I2PQ[i] == pq) {
            *pI = i;
            answer = true;
            break;
        }
    }

    return answer;
}


