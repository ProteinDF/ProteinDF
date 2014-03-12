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

#ifndef DFGRIDFREEXC_H
#define DFGRIDFREEXC_H

#include "DfObject.h"
#include "DfXCFunctional.h"
#include "DfOverlapX.h"
#include "DfXMatrix.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "DfFunctional.h"

class DfGridFreeXC : public DfObject 
{
public:
    struct IndexPair2 {
    public:
        explicit IndexPair2(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
            if (this->index1_ > this->index2_) {
                std::swap(this->index1_, this->index2_);
            }
        }
        
        std::size_t index() const {
            // 'U' format
            assert(this->index1_ <= this->index2_);
            return this->index1_ + this->index2_ * (this->index2_ +1) / 2;
        }
        
        bool operator<(const IndexPair2& rhs) const {
            return (this->index() < rhs.index());
        }
        
        index_type index1() const {
            return this->index1_;
        }
        
        index_type index2() const {
            return this->index2_;
        }
        
    private:
        index_type index1_;
        index_type index2_;
    };
    
    struct IndexPair4 {
    public:
        explicit IndexPair4(index_type i1 =0, index_type i2 =0,
                            index_type i3 =0, index_type i4 =0) 
            : indexPair1_(i1, i2), indexPair2_(i3, i4) {
            if (this->indexPair1_.index() > this->indexPair2_.index()) {
                std::swap(this->indexPair1_, this->indexPair2_);
            }
        }
        
        std::size_t index() const {
            std::size_t pair_index1 = this->indexPair1_.index();
            std::size_t pair_index2 = this->indexPair2_.index();
            assert(pair_index1 <= pair_index2);
            return pair_index1 + pair_index2 * (pair_index2 +1) / 2;
        }
        
        bool operator<(const IndexPair4& rhs) const {
            return (this->index() < rhs.index());
        }
        
        bool operator==(const IndexPair4& rhs) const {
            return (this->index() == rhs.index());
        }
        
        index_type index1() const {
            return this->indexPair1_.index1();
        }
        
        index_type index2() const {
            return this->indexPair1_.index2();
        }
        
        index_type index3() const {
            return this->indexPair2_.index1();
        }
        
        index_type index4() const {
            return this->indexPair2_.index2();
        }
    
    private:
        IndexPair2 indexPair1_;
        IndexPair2 indexPair2_;
    };

    typedef std::vector<IndexPair2> PQ_PairArray;
    
public:
    DfGridFreeXC(TlSerializeData* pPdfParam);
    virtual ~DfGridFreeXC();
    
public:
    virtual void preprocessBeforeSCF();
    void buildFxc();

protected:
    DfOverlapX* getDfOverlapObject();
    DfXMatrix* getDfXMatrixObject();

    template<class DfOverlapClass,
             class DfXMatrixClass,
             class SymmetricMatrixType, class MatrixType>
    void preprocessBeforeSCF_templ();


    virtual void buildFxc_LDA();

    virtual void buildFxc_GGA();

protected:
    static const int MAX_SHELL_TYPE;
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;

    struct ShellPair {
    public:
        ShellPair(index_type index1 =0, index_type index2 =0) : shellIndex1(index1), shellIndex2(index2) {
        }
        
    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    typedef std::vector<ShellPair> ShellPairArray;
    typedef std::vector<ShellPairArray> ShellPairArrayTable;
    
protected:
    // virtual void getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    // virtual void getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    
    // TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);
    // void getM_part(const TlOrbitalInfoObject& orbitalInfo,
    //                const std::vector<DfTaskCtrl::Task4>& taskList,
    //                const TlMatrixObject& P, TlMatrixObject* pM);
    void getM_part(const TlOrbitalInfoObject& orbitalInfo_PQ,
                   const TlOrbitalInfoObject& orbitalInfo_RS,
                   const std::vector<DfTaskCtrl::Task4>& taskList,
                   const TlMatrixObject& P, TlMatrixObject* pM);
    // void storeM(const index_type shellIndexP, const int maxStepsP,
    //             const index_type shellIndexQ, const int maxStepsQ,
    //             const index_type shellIndexR, const int maxStepsR,
    //             const index_type shellIndexS, const int maxStepsS,
    //             const DfOverlapEngine& engine,
    //             const TlMatrixObject& P,
    //             TlMatrixObject* pM);
    void storeM_A(const index_type shellIndexP, const int maxStepsP,
                  const index_type shellIndexQ, const int maxStepsQ,
                  const index_type shellIndexR, const int maxStepsR,
                  const index_type shellIndexS, const int maxStepsS,
                  const DfOverlapEngine& engine,
                  const TlMatrixObject& P,
                  TlMatrixObject* pM);

    virtual void createEngines();
    virtual void destroyEngines();
    virtual DfTaskCtrl* getDfTaskCtrlObject() const;
    virtual void finalize(TlSymmetricMatrix* pMtx);

    void get_F_lamda(const TlVector lamda,
                     TlMatrixObject* pF_lamda,
                     TlMatrixObject* pE_lamda);

    // void get_F_lamda_GGA(const TlVector lambda_f,
    //                      const TlVector lambda_g,
    //                      TlSymmetricMatrix* pE_f,
    //                      TlSymmetricMatrix* pE_g,
    //                      TlSymmetricMatrix* pF_f_rho,
    //                      TlSymmetricMatrix* pF_g_rho,
    //                      TlSymmetricMatrix* pF_f_gaa,
    //                      TlSymmetricMatrix* pF_g_gaa);

    void getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);
    ShellPairArrayTable getShellPairArrayTable(const ShellArrayTable& shellArrayTable);

protected:
    template<class SymmetricMatrixType>
    SymmetricMatrixType getPMatrix(const RUN_TYPE runType);

    template<class DfOverlapClass,
             class DfXMatrixClass,
             class SymmetricMatrixType, class MatrixType>
    void buildFxc_LDA_method();

    template<class DfOverlapClass,
             class DfCD_class,
             class SymmetricMatrixType, class MatrixType>
    void buildFxc_LDA_runtype(const RUN_TYPE runType);


    template<class DfOverlapClass,
             class DfXMatrixClass,
             class SymmetricMatrixType, class MatrixType>
    void buildFxc_GGA_method();

    template<class DfOverlapClass,
             class DfCD_class,
             class SymmetricMatrixType, class MatrixType>
    void buildFxc_GGA_runtype(const RUN_TYPE runType);

public:
    TlMatrix getForce();

    // virtual void calcCholeskyVectors_onTheFly();

protected:
    TlMatrix selectGradMat(const TlMatrix& input, const int atomIndex);

    // void calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
    //                    PQ_PairArray *pI2PQ,
    //                    TlVector *pDiagonals);
    // void calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
    //                           TlSparseSymmetricMatrix *pSchwartzTable,
    //                           TlSparseSymmetricMatrix *pDiagonalMat,
    //                           PQ_PairArray *pI2PQ);
    // void saveI2PQ(const PQ_PairArray& I2PQ);
    // void saveL(const TlMatrix& L);

    // std::vector<double>
    // getSuperMatrixElements(const index_type G_row,
    //                        const std::vector<index_type>& G_col_list,
    //                        const PQ_PairArray& I2PQ,
    //                        const TlSparseSymmetricMatrix& schwartzTable);

    // std::vector<DfGridFreeXC::IndexPair4> 
    // getCalcList(const index_type G_row,
    //             const std::vector<index_type>& G_col_list,
    //             const PQ_PairArray& I2PQ);

    // void calcElements(const std::vector<IndexPair4>& calcList,
    //                   const TlSparseSymmetricMatrix& schwartzTable);

    // bool isAliveBySchwartzCutoff(const index_type shellIndexP,
    //                              const index_type shellIndexQ,
    //                              const index_type shellIndexR,
    //                              const index_type shellIndexS,
    //                              const int shellQuartetType,
    //                              const TlSparseSymmetricMatrix& schwarzTable,
    //                              const double threshold);
    // void initializeCutoffStats();
    // void schwartzCutoffReport();
    // mutable std::vector<unsigned long> cutoffAll_schwartz_;
    // mutable std::vector<unsigned long> cutoffAlive_schwartz_;

    // std::vector<double>
    // setElements(const index_type G_row,
    //             const std::vector<index_type> G_col_list,
    //             const PQ_PairArray& I2PQ);

    // void getM_byCD(TlSymmetricMatrix* pM);
    //TlSymmetricMatrix getPMatrix();
    TlMatrix getL();
    PQ_PairArray getI2PQ();
    void divideCholeskyBasis(const index_type numOfCBs,
                             index_type *pStart, index_type *pEnd);
    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const PQ_PairArray& I2PQ);

    DfFunctional_GGA* getFunctionalGGA();
    
protected:
    static const double ONE_THIRD; // = 1.0 / 3.0

    DfOverlapEngine* pOvpEngines_;
    TlOrbitalInfo orbitalInfo_;

    index_type numOfPQs_;
    double tau_;
    double epsilon_;

    /// キャッシュ
    typedef std::map<IndexPair4, std::vector<double> > ElementsCacheType;
    ElementsCacheType elements_cache_;

    ///
    bool debugSaveM_;

    /// 規格直交化ルーチンにCanonical Orthogonalizeを使う
    /// true: canonical
    /// false: lowdin
    bool isCanonicalOrthogonalize_;
    // int GF_mode_;
};


template<class DfOverlapClass,
         class DfXMatrixClass,
         class SymmetricMatrixType, class MatrixType>
void DfGridFreeXC::preprocessBeforeSCF_templ()
{
    const DfXCFunctional xcFunc(this->pPdfParam_);
    const DfXCFunctional::FUNCTIONAL_TYPE funcType = xcFunc.getFunctionalType();
    bool isGGA = (funcType == DfXCFunctional::GGA);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]); // orbital用
    const TlOrbitalInfo orbitalInfo_GF((*(this->pPdfParam_))["coordinates"],
                                       (*(this->pPdfParam_))["basis_sets_GF"]); // GridFree用

    DfOverlapClass dfOvp(this->pPdfParam_);
    if (this->isDedicatedBasisForGridFree_) {
        {
            this->log_.info("build S_gf matrix");
            SymmetricMatrixType gfS;
            dfOvp.getOvpMat(orbitalInfo_GF, &gfS);
            
            this->log_.info("build V matrix");
            
            DfXMatrixClass dfXMatrix(this->pPdfParam_);
            MatrixType gfV;
            if (this->isCanonicalOrthogonalize_) {
                this->log_.info("orthogonalize method: canoncal");
                dfXMatrix.canonicalOrthogonalize(gfS, &gfV, NULL);
            } else {
                this->log_.info("orthogonalize method: lowdin");
                dfXMatrix.lowdinOrthogonalize(gfS, &gfV, NULL);
            }
            this->log_.info("save V matrix");
            DfObject::saveGfVMatrix(gfV);
        }
        
        {
            this->log_.info("build S~ matrix: start");
            MatrixType gfStilde;
            dfOvp.getTransMat(this->orbitalInfo_,
                              orbitalInfo_GF,
                              &gfStilde);
            this->log_.info("build S~ matrix: save");
            DfObject::saveGfStildeMatrix(gfStilde);
            
            this->log_.info("build (S~)^-1 matrix: start");
            SymmetricMatrixType Sinv = DfObject::getSpqMatrix<SymmetricMatrixType>();
            Sinv.inverse();
            
            gfStilde.transpose();
            const MatrixType gfOmega = gfStilde * Sinv;
            this->log_.info("build (S~)^-1 matrix: save");
            DfObject::saveGfOmegaMatrix(gfOmega);
            this->log_.info("build (S~)^-1 matrix: finish");
        }
    }

    // for GGA
    if (isGGA) {
        this->log_.info("build gradient matrix: start");
        MatrixType Gx, Gy, Gz;
        if (this->isDedicatedBasisForGridFree_) {
            dfOvp.getGradient(orbitalInfo_GF, &Gx, &Gy, &Gz);
        } else {
            dfOvp.getGradient(orbitalInfo, &Gx, &Gy, &Gz);
        }
        this->log_.info("build gradient matrix: save");
        this->saveDipoleVelocityIntegralsXMatrix(Gx);
        this->saveDipoleVelocityIntegralsYMatrix(Gy);
        this->saveDipoleVelocityIntegralsZMatrix(Gz);
        this->log_.info("build gradient matrix: finish");
    }
}

template<class SymmetricMatrixType>
SymmetricMatrixType DfGridFreeXC::getPMatrix(const RUN_TYPE runType)
{
    SymmetricMatrixType P = DfObject::getPpqMatrix<SymmetricMatrixType>(runType,
                                                                        this->m_nIteration -1);
    if (runType ==RUN_RKS) {
        P *= 0.5;
    }

    return P;
}

template<class DfOverlapClass,
         class DfCD_class,
         class SymmetricMatrixType, class MatrixType>
void DfGridFreeXC::buildFxc_LDA_method()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:   
        this->buildFxc_LDA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_RKS);
        break;

    case METHOD_UKS:
        this->buildFxc_LDA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_UKS_ALPHA);
        this->buildFxc_LDA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->log_.critical("NOT in support");
        break;

    default:
        this->log_.critical("wrong program");
        break;
    }
}


template<class DfOverlapClass,
         class DfCD_class,
         class SymmetricMatrixType, class MatrixType>
void DfGridFreeXC::buildFxc_LDA_runtype(const RUN_TYPE runType)
{
    this->log_.info("build Fxc by grid-free method: functional type is LDA.");

    std::string basisset_param = "basis_sets";
    if (this->isDedicatedBasisForGridFree_) {
        basisset_param = "basis_sets_GF";
    }
    const TlOrbitalInfo orbitalInfo_GF((*(this->pPdfParam_))["coordinates"],
                                       (*(this->pPdfParam_))[basisset_param]);

    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfGfOrbs = orbitalInfo_GF.getNumOfOrbitals();
    this->log_.info(TlUtils::format("AOs = %d", numOfAOs));
    this->log_.info(TlUtils::format("auxAOs for GF = %d", numOfGfOrbs));

    const SymmetricMatrixType P = this->getPMatrix<SymmetricMatrixType>(runType);

    SymmetricMatrixType M;
    if (this->XC_engine_ == XC_ENGINE_GRIDFREE_CD) {
        this->log_.info("begin to create M matrix based on CD.");
        {
            DfCD_class dfCD(this->pPdfParam_);
            dfCD.getM(P, &M);
        }
    } else {
        this->log_.info("begin to create M matrix using 4-center overlap.");
        DfOverlapClass dfOvp(this->pPdfParam_);
        if (this->isDedicatedBasisForGridFree_) {
            dfOvp.getM_A(P, &M);
        } else {
            dfOvp.getM(P, &M);
        }
    }
    if (this->debugSaveM_) {
        M.save(TlUtils::format("fl_Work/debug_M.%d.mat", this->m_nIteration));
    }

    // DEBUG
    // {
    //     this->getM_exact(P, &M);
    //     M.save("M_exact.mat");
    // }

    this->log_.info("begin to generate Fxc using grid-free method.");
    // tV * S * V == I
    MatrixType S;
    MatrixType V;
    if (this->isDedicatedBasisForGridFree_) {
        S = DfObject::getGfStildeMatrix<MatrixType>();
        V = DfObject::getGfVMatrix<MatrixType>();
    } else {
        S = DfObject::getSpqMatrix<SymmetricMatrixType>();
        V = DfObject::getXMatrix<MatrixType>();
    }
    MatrixType tV = V;
    tV.transpose();
    
    SymmetricMatrixType M_tilda = tV * M * V;
    
    MatrixType U;
    TlVector lambda;
    M_tilda.diagonal(&lambda, &U);
    
    SymmetricMatrixType F_lambda(lambda.getSize());
    SymmetricMatrixType E_lambda(lambda.getSize());
    this->get_F_lamda(lambda, &F_lambda, &E_lambda);
    
    MatrixType SVU = S * V * U;
    MatrixType UVS = SVU;
    UVS.transpose();
    
    // save
    {
        SymmetricMatrixType Fxc = SVU * F_lambda * UVS;
        DfObject::saveFxcMatrix(runType, this->m_nIteration, Fxc);
    }
    {
        SymmetricMatrixType Exc = SVU * E_lambda * UVS;
        DfObject::saveExcMatrix(runType, this->m_nIteration, Exc);
    }
}


template<class DfOverlapClass,
         class DfCD_class,
         class SymmetricMatrixType, class MatrixType>
void DfGridFreeXC::buildFxc_GGA_method()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:   
        this->buildFxc_GGA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_RKS);
        break;

    case METHOD_UKS:
        this->buildFxc_GGA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_UKS_ALPHA);
        this->buildFxc_GGA_runtype<DfOverlapClass, DfCD_class,
                                   SymmetricMatrixType, MatrixType>(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->log_.critical("NOT in support");
        break;

    default:
        this->log_.critical("wrong program");
        break;
    }
}


template<class DfOverlapClass,
         class DfCD_class,
         class SymmetricMatrixType, class MatrixType>
void DfGridFreeXC::buildFxc_GGA_runtype(const RUN_TYPE runType)
{
    this->log_.info("build Fxc by grid-free method: functional type is GGA.");

    std::string basisset_param = "basis_sets";
    if (this->isDedicatedBasisForGridFree_) {
        basisset_param = "basis_sets_GF";
    }
    const TlOrbitalInfo orbitalInfo_GF((*(this->pPdfParam_))["coordinates"],
                                       (*(this->pPdfParam_))[basisset_param]);
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfGfOrbs = orbitalInfo_GF.getNumOfOrbitals();
    this->log_.info(TlUtils::format("AOs = %d", numOfAOs));
    this->log_.info(TlUtils::format("auxAOs for GF = %d", numOfGfOrbs));

    // RKS
    // const SymmetricMatrixType PA = 0.5 * DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS,
    //                                                                              this->m_nIteration -1);
    const SymmetricMatrixType P = this->getPMatrix<SymmetricMatrixType>(runType);
    assert(P.getNumOfRows() == numOfAOs);

    SymmetricMatrixType M;
    if (this->XC_engine_ == XC_ENGINE_GRIDFREE_CD) {
        this->log_.info("begin to create M matrix based on CD.");
        {
            DfCD_class dfCD(this->pPdfParam_);
            dfCD.getM(P, &M);
            // M *= 2.0;
        }
    } else {
        this->log_.info("begin to create M matrix using 4-center overlap.");
        DfOverlapClass dfOvp(this->pPdfParam_);
        if (this->isDedicatedBasisForGridFree_) {
            dfOvp.getM_A(P, &M);
        } else {
            dfOvp.getM(P, &M);
        }
    }
    if (this->debugSaveM_) {
        M.save(TlUtils::format("fl_Work/debug_M.%d.mat", this->m_nIteration));
    }

    this->log_.info("begin to generate Fxc using grid-free method.");
    // M~(=V^t * M * V) および SVU(=SVU, but now only SV)の作成
    MatrixType S;
    MatrixType V;
    if (this->isDedicatedBasisForGridFree_) {
        S = DfObject::getGfStildeMatrix<MatrixType>();
        V = DfObject::getGfVMatrix<MatrixType>();
    } else {
        S = DfObject::getSpqMatrix<SymmetricMatrixType>();
        V = DfObject::getXMatrix<MatrixType>();
    }

    const index_type numOfGFOrthNormBasis = V.getNumOfCols();
    this->log_.info(TlUtils::format("orthonormal basis = %d", numOfGFOrthNormBasis));
    MatrixType Vt = V;
    Vt.transpose();

    MatrixType St = S;
    St.transpose();

    SymmetricMatrixType Mtilde = Vt * M * V;
    //Mtilde.save("Mtilde.mat");

    TlVector lambda;
    MatrixType U;
    Mtilde.diagonal(&lambda, &U);
    //lambda.save("lambda.vct");
    //U.save("U.mat");
    MatrixType Ut = U;
    Ut.transpose();
    assert(lambda.getSize() == numOfGFOrthNormBasis);
    assert(U.getNumOfRows() == numOfGFOrthNormBasis);
    assert(U.getNumOfCols() == numOfGFOrthNormBasis);

    // M[rho^(-1/3)]
    SymmetricMatrixType Mtilde_13;
    {
        SymmetricMatrixType lambda_13(numOfGFOrthNormBasis);
        for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
            double v_13 = 0.0;
            const double v = lambda.get(i);
            if (v > 1.0E-16) {
                v_13 = std::pow(v, - ONE_THIRD);
            }
            lambda_13.set(i, i, v_13);
        }
        Mtilde_13 = U * lambda_13 * Ut;
    }
    //Mtilde_13.save("Mtilde_13.mat");

    // GGA用gradient
    MatrixType Gx = this->getDipoleVelocityIntegralsXMatrix<MatrixType>();
    MatrixType Gy = this->getDipoleVelocityIntegralsYMatrix<MatrixType>();
    MatrixType Gz = this->getDipoleVelocityIntegralsZMatrix<MatrixType>();
    Gx *= -1.0;
    Gy *= -1.0;
    Gz *= -1.0;
    const MatrixType DX = Vt * Gx * V;
    const MatrixType DY = Vt * Gy * V;
    const MatrixType DZ = Vt * Gz * V;
    // DX.save("DX.mat");
    // DY.save("DY.mat");
    // DZ.save("DZ.mat");

    MatrixType DXt = DX;
    DXt.transpose();
    MatrixType DYt = DY;
    DYt.transpose();
    MatrixType DZt = DZ;
    DZt.transpose();
    const MatrixType RTX = 3.0*(DXt * Mtilde_13 + Mtilde_13 * DX);
    const MatrixType RTY = 3.0*(DYt * Mtilde_13 + Mtilde_13 * DY);
    const MatrixType RTZ = 3.0*(DZt * Mtilde_13 + Mtilde_13 * DZ);
    MatrixType RTXt = RTX;
    RTXt.transpose();
    MatrixType RTYt = RTY;
    RTYt.transpose();
    MatrixType RTZt = RTZ;
    RTZt.transpose();

    // RX2 := M[{nabla rho / rho^(-4/3)}^2]
    const SymmetricMatrixType RX2 = RTXt*RTX + RTYt*RTY + RTZt*RTZ;
    // RTX.save("RTX.mat");
    // RTY.save("RTY.mat");
    // RTZ.save("RTZ.mat");
    // RX2.save("RX2.mat");

    // RZ2 := [nabla rho / rho^(4/3)] * (DX + DY + DZ)
    const MatrixType RZ2 = RTX*DX + RTY*DY + RTZ*DZ;
    //RZ2.save("RZ2.mat");
    MatrixType RZ2t = RZ2;
    RZ2t.transpose();

    TlVector x2;
    MatrixType Ux2;
    RX2.diagonal(&x2, &Ux2);
    //x2.save("x2.vct");
    //Ux2.save("Ux2.mat");
    MatrixType Ux2t = Ux2;
    Ux2t.transpose();

    // ------------------
    assert(lambda.getSize() == numOfGFOrthNormBasis);
    // TlVector rhoAs(numOfGfOrbs);
    // TlVector xAs(numOfGfOrbs);
    TlVector rhoAs(numOfGFOrthNormBasis);
    TlVector xAs(numOfGFOrthNormBasis);
    //for (index_type i = 0; i < numOfGfOrbs; ++i) {
    for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
        const double rho_value = lambda[i];
        const double rho = (rho_value > 1.0E-16) ? rho_value : 0.0;
        rhoAs[i] = rho;
        
        const double x2_value = x2.get(i);
        const double x = (x2_value > 1.0E-16) ? std::sqrt(x2_value) : 0.0;
        xAs[i] = x;
    }
    const TlVector rhoBs = rhoAs;
    const TlVector xBs = xAs;
    
    DfFunctional_GGA* pFunc = this->getFunctionalGGA();
    // Fxc -------------------------------------------------------------
    SymmetricMatrixType FxcA(numOfAOs); // alpha spin
    SymmetricMatrixType FxcB(numOfAOs); // beta spin
    {
        DerivativeFunctionalSets dfs = pFunc->getDerivativeFunctional_GF(rhoAs, rhoBs, xAs, xBs);
        
        TlVector rhoAA43(numOfGFOrthNormBasis);
        TlVector rhoAB43(numOfGFOrthNormBasis);
        TlVector rhoBB43(numOfGFOrthNormBasis);
        for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
            const double rhoA = lambda[i];
            if (rhoA > 1.0E-16) {
                const double rhoA43 = std::pow(rhoA, 4.0/3.0);
                const double rhoB43 = rhoA43; // RKS
                rhoAA43[i] = rhoA43;
                rhoBB43[i] = rhoB43;
                rhoAB43[i] = std::sqrt(rhoA43 * rhoB43);
            }
        }

        MatrixType FxcA_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        MatrixType FxcB_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        SymmetricMatrixType diag_RAR(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_RAX(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_RBR(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_RBX(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GAAR(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GAAX(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GABR(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GABX(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GBBR(numOfGFOrthNormBasis);
        SymmetricMatrixType diag_GBBX(numOfGFOrthNormBasis);
        const int numOfTerms = pFunc->getNumOfDerivativeFunctionalTerms();
        for (int term = 0; term < numOfTerms; ++term) {
            for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
                diag_RAR(i, i) = dfs.rFrRhoA_R(term, i);
                diag_RAX(i, i) = dfs.rFrRhoA_X(term, i);
                diag_RBR(i, i) = dfs.rFrRhoB_R(term, i);
                diag_RBX(i, i) = dfs.rFrRhoB_X(term, i);

                diag_GAAR(i, i) = dfs.rFrGAA_R(term, i) * rhoAA43[i];
                diag_GAAX(i, i) = dfs.rFrGAA_X(term, i);
                diag_GABR(i, i) = dfs.rFrGAB_R(term, i) * rhoAB43[i];
                diag_GABX(i, i) = dfs.rFrGAB_X(term, i);
                diag_GBBR(i, i) = dfs.rFrGBB_R(term, i) * rhoBB43[i];
                diag_GBBX(i, i) = dfs.rFrGBB_X(term, i);
            }

            // alpha spin ------------
            {
                const SymmetricMatrixType Fxc_RR = U * diag_RAR * Ut;
                const SymmetricMatrixType Fxc_RX = Ux2 * diag_RAX * Ux2t;
                MatrixType Fxc_tilde_term1 = 0.5 * (Fxc_RR * Fxc_RX + Fxc_RX * Fxc_RR);
                FxcA_tilde += Fxc_tilde_term1;
            }

            MatrixType Fxc_GAA;
            {
                const MatrixType Fxc_GAAR = U * diag_GAAR * Ut;
                const MatrixType Fxc_GAAX = Ux2 * diag_GAAX * Ux2t;
                Fxc_GAA = 2.0 * 0.5 * (Fxc_GAAR * Fxc_GAAX + Fxc_GAAX * Fxc_GAAR);
            }

            MatrixType Fxc_GAB;
            {
                const SymmetricMatrixType Fxc_GABR = U * diag_GABR * Ut;
                const SymmetricMatrixType Fxc_GABX = Ux2 * diag_GABX * Ux2t;
                Fxc_GAB = 0.5 * (Fxc_GABR * Fxc_GABX + Fxc_GABX * Fxc_GABR);
            }

            MatrixType Fxc_GBB;
            {
                const SymmetricMatrixType Fxc_GBBR = U * diag_GBBR * Ut;
                const SymmetricMatrixType Fxc_GBBX = Ux2 * diag_GBBX * Ux2t;
                Fxc_GBB = 2.0 * 0.5 * (Fxc_GBBR * Fxc_GBBX + Fxc_GBBX * Fxc_GBBR);
            }

            {
                MatrixType FxcA_term2 = Fxc_GAA + Fxc_GAB;
                MatrixType FxcA_term2t = FxcA_term2;
                FxcA_term2t.transpose();
                MatrixType FxcA_tilde_term2 = FxcA_term2t * RZ2 + RZ2t * FxcA_term2;
                FxcA_tilde += FxcA_tilde_term2;
            }
            {
                MatrixType FxcB_term2 = Fxc_GBB + Fxc_GAB;
                MatrixType FxcB_term2t = FxcB_term2;
                FxcB_term2t.transpose();
                MatrixType FxcB_tilde_term2 = FxcB_term2t * RZ2 + RZ2t * FxcB_term2;
                FxcB_tilde += FxcB_tilde_term2;
            }
        }
        FxcA = S * V * FxcA_tilde * Vt * St;
        FxcB = S * V * FxcB_tilde * Vt * St;
    }
    //DfObject::saveFxcMatrix(RUN_RKS, this->m_nIteration, SymmetricMatrixType(FxcA));
    DfObject::saveFxcMatrix(runType, this->m_nIteration, SymmetricMatrixType(FxcA));

    // Exc -------------------------------------------------------------
    SymmetricMatrixType ExcA(numOfAOs);
    //SymmetricMatrixType ExcB(numOfAOs);
    {
        const FunctionalSets fs = pFunc->getFunctional_GF(rhoAs, rhoBs, xAs, xBs);

        MatrixType ExcA_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        //MatrixType ExcB_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        {
            SymmetricMatrixType diag_AR(numOfGFOrthNormBasis);
            SymmetricMatrixType diag_AX(numOfGFOrthNormBasis);
            //SymmetricMatrixType diag_BR(numOfGFOrthNormBasis);
            //SymmetricMatrixType diag_BX(numOfGFOrthNormBasis);
            const int numOfTerms = pFunc->getNumOfFunctionalTerms();
            for (int term = 0; term < numOfTerms; ++term) {
                for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
                    const double rho = rhoAs[i] + rhoBs[i];
                    if (rho > 1.0E-16) {
                        const double inv_rho = 1.0 / rho;
                        diag_AR(i, i) = fs.FA_termR(term, i) * inv_rho;
                        diag_AX(i, i) = fs.FA_termX(term, i);
                        // diag_BR(i, i) = fs.FB_termR(term, i) * inv_rho;
                        // diag_BX(i, i) = fs.FB_termX(term, i);
                    }
                }
                
                const SymmetricMatrixType ExcA_R = U * diag_AR * Ut;
                const SymmetricMatrixType ExcA_X = Ux2 * diag_AX * Ux2t;
                const MatrixType ExcA_tilde_term = ExcA_R * ExcA_X + ExcA_X * ExcA_R;
                // const SymmetricMatrixType ExcB_R = U * diag_BR * Ut;
                // const SymmetricMatrixType ExcB_X = Ux2 * diag_BX * Ux2t;
                // const MatrixType ExcB_tilde_term = ExcB_R * ExcB_X + ExcB_X * ExcB_R;
                
                ExcA_tilde += ExcA_tilde_term;
                // ExcB_tilde += ExcB_tilde_term;
            }
        }

        ExcA = S * V * ExcA_tilde * Vt * St;
        // ExcB = S * V * ExcB_tilde * Vt * St;
        // Exc *= 2.0; // means RKS
    }
    //DfObject::saveExcMatrix(RUN_RKS, this->m_nIteration, SymmetricMatrixType(ExcA));
    DfObject::saveExcMatrix(runType, this->m_nIteration, SymmetricMatrixType(ExcA));

    delete pFunc;
    pFunc = NULL;
}

#endif // DFGRIDFREEXC_H

