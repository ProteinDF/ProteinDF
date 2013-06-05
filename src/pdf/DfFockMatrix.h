#ifndef DFFOCKMATRIX_H
#define DFFOCKMATRIX_H

#include <cassert>
#include "DfObject.h"
#include "DfXCFunctional.h"
#include "CnError.h"
#include "TlVector.h"

class DfFockMatrix : public DfObject {
public:
    DfFockMatrix(TlSerializeData* pPdfParam);
    virtual ~DfFockMatrix();

public:
    void DfFockMatrixMain();

protected:
    virtual void mainDIRECT_RKS();
    virtual void mainDIRECT_UKS();
    virtual void mainDIRECT_ROKS();

    virtual void setXC_RI(RUN_TYPE nRunType, TlSymmetricMatrix& F);
    virtual void setXC_DIRECT(RUN_TYPE nRunType, TlSymmetricMatrix& F);
    virtual void setCoulomb(METHOD_TYPE nMethodType, TlSymmetricMatrix& F);
    virtual TlSymmetricMatrix getFpqMatrix(RUN_TYPE nRunType, int nIteration);

    virtual TlVector getRho(RUN_TYPE nRunType, int nIteration);
    virtual TlVector getMyu(RUN_TYPE nRunType, int nIteration);
    //virtual void saveFpqMatrix(RUN_TYPE nRunType, const TlSymmetricMatrix& F);

protected:
    template<typename SymmetricMatrixType>
    void mainDIRECT_RKS();

    template<typename SymmetricMatrixType>
    void mainDIRECT_UKS();

    template<typename MatrixType, typename SymmetricMatrixType, typename VectorType, class DfEriClass, class DfOverlapClass>
    void mainDIRECT_ROKS();

    template<typename SymmetricMatrixType, typename VectorType, class DfOverlapClass>
    void setXC_RI(const RUN_TYPE nRunType, SymmetricMatrixType& F);

    template<typename SymmetricMatrixType>
    void setXC_DIRECT(const RUN_TYPE nRunType, SymmetricMatrixType& F);

    template<typename SymmetricMatrixType, typename VectorType, class DfEriClass>
    void setCoulomb(const METHOD_TYPE nMethodType, SymmetricMatrixType& F);

    template<typename SymmetricMatrixType>
    void setHpq(const RUN_TYPE nRunType, SymmetricMatrixType& F);

protected:
    bool isUseNewEngine_;
};

// =====================================================================
// template
template<typename SymmetricMatrixType>
void DfFockMatrix::mainDIRECT_RKS()
{
    SymmetricMatrixType F(this->m_nNumOfAOs);

    {
        SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
        F += Hpq;
    }
    if (this->m_nNumOfDummyAtoms > 0) {
        SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();
        const int chargeExtrapolateNumber = std::max(this->chargeExtrapolateNumber_, 1);
        const int times = std::min(this->m_nIteration, chargeExtrapolateNumber);
        F += static_cast<int>(times) * Hpq2;
    }

    {
        SymmetricMatrixType J = DfObject::getJMatrix<SymmetricMatrixType>(this->m_nIteration);
        F += J;
    }

    {
        SymmetricMatrixType Fxc = DfObject::getFxcMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
        F += Fxc;
    }

    {
        const DfXCFunctional dfXCFunctional(this->pPdfParam_);
        if (dfXCFunctional.isHybridFunctional() == true) {
            const double coef = dfXCFunctional.getFockExchangeCoefficient();
            this->log_.info(TlUtils::format("coefficient of K: %f", coef));
            SymmetricMatrixType K = DfObject::getHFxMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
            K *= coef;
            F += K;
        }
    }

    DfObject::saveFpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration, F);
}

template<typename SymmetricMatrixType>
void DfFockMatrix::mainDIRECT_UKS()
{
    SymmetricMatrixType FA(this->m_nNumOfAOs);
    SymmetricMatrixType FB(this->m_nNumOfAOs);

    {
        SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
        FA += Hpq;
        FB += Hpq;
    }
    if (this->m_nNumOfDummyAtoms > 0) {
        SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();
        const int chargeExtrapolateNumber = std::max(this->chargeExtrapolateNumber_, 1);
        const int times = std::min(this->m_nIteration, chargeExtrapolateNumber);
        FA += static_cast<int>(times) * Hpq2;
        FB += static_cast<int>(times) * Hpq2;
    }

    {
        SymmetricMatrixType J = DfObject::getJMatrix<SymmetricMatrixType>(this->m_nIteration);
        FA += J;
        FB += J;
    }

    {
        SymmetricMatrixType FxcA = DfObject::getFxcMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration);
        FA += FxcA;
    }
    {
        SymmetricMatrixType FxcB = DfObject::getFxcMatrix<SymmetricMatrixType>(RUN_UKS_BETA, this->m_nIteration);
        FB += FxcB;
    }

    {
        const DfXCFunctional dfXCFunctional(this->pPdfParam_);
        if (dfXCFunctional.isHybridFunctional() == true) {
            const double coef = dfXCFunctional.getFockExchangeCoefficient();
            this->log_.info(TlUtils::format("coefficient of K: %f", coef));
            {
                SymmetricMatrixType KA = DfObject::getHFxMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA,
                                                                                     this->m_nIteration);
                KA *= coef;
                FA += KA;
            }
            {
                SymmetricMatrixType KB = DfObject::getHFxMatrix<SymmetricMatrixType>(RUN_UKS_BETA,
                                                                                     this->m_nIteration);
                KB *= coef;
                FB += KB;
            }
        }
    }

    DfObject::saveFpqMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration, FA);
    DfObject::saveFpqMatrix<SymmetricMatrixType>(RUN_UKS_BETA,  this->m_nIteration, FB);
}


template<typename MatrixType, typename SymmetricMatrixType, typename VectorType, class DfEriClass, class DfOverlapClass>
void DfFockMatrix::mainDIRECT_ROKS() 
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    SymmetricMatrixType Fc(numOfAOs);
    SymmetricMatrixType Fo(numOfAOs);

    {
        SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
        Fc += Hpq;
    }
    if (this->m_nNumOfDummyAtoms > 0) {
        SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();
        const int chargeExtrapolateNumber = std::max(this->chargeExtrapolateNumber_, 1);
        const int times = std::min(this->m_nIteration, chargeExtrapolateNumber);
        Fc += static_cast<int>(times) * Hpq2;
    }

    {
        SymmetricMatrixType J = DfObject::getJMatrix<SymmetricMatrixType>(this->m_nIteration);
        Fo += 0.5 * Fc;
        Fc += J;
    }

    {
        SymmetricMatrixType Fxc_o = DfObject::getFxcMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN, this->m_nIteration);
        Fo += Fxc_o;
    }
    {
        SymmetricMatrixType Fxc_c = DfObject::getFxcMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, this->m_nIteration);
        Fc += Fxc_c;
    }

    {
        const DfXCFunctional dfXCFunctional(this->pPdfParam_);
        if (dfXCFunctional.isHybridFunctional() == true) {
            const double coef = dfXCFunctional.getFockExchangeCoefficient();
            this->log_.info(TlUtils::format("coefficient of K: %f", coef));
            {
                SymmetricMatrixType Ko = DfObject::getHFxMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN, this->m_nIteration);
                Ko *= coef;
                Fo += Ko;
            }
            {
                SymmetricMatrixType Kc = DfObject::getHFxMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, this->m_nIteration);
                Kc *= coef;
                Fc += Kc;
            }
        }
    }

    // -------------------------------------------------------------------------
    MatrixType SDc, DcS, SDo, DoS;
    {
        const SymmetricMatrixType S = DfObject::getSpqMatrix<SymmetricMatrixType>();
        const SymmetricMatrixType Dc = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, this->m_nIteration -1);
        const SymmetricMatrixType Do = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN,  this->m_nIteration -1);
        SDc = S * Dc;
        DcS = Dc * S;
        SDo = S * Do;
        DoS = Do * S;
    }

    SymmetricMatrixType F(numOfAOs);
    {
        SymmetricMatrixType E(numOfAOs);
        for (index_type i = 0; i < numOfAOs; ++i) {
            E.set(i, i, 1.0);
        }
        const MatrixType E_SDc = E - SDc;
        const MatrixType E_DcS = E - DcS;
        const MatrixType E_SDo = E - SDo;
        const MatrixType E_DoS = E - DoS;
        
        F += E_SDo * Fc * E_DoS;
        F += E_SDc * Fo * E_DcS;
    }

    {
        const SymmetricMatrixType FcFo = Fc - Fo;
        F += SDc * FcFo * DoS;
        F += SDo * FcFo * DcS;
    }

    DfObject::saveFpqMatrix(RUN_ROKS, this->m_nIteration, F);
}


// template<typename MatrixType, typename SymmetricMatrixType, typename VectorType, class DfEriClass, class DfOverlapClass>
// void DfFockMatrix::mainDIRECT_ROKS()
// {
//     this->logger("ROKS Direct scheme method is employed\n");
//     assert(this->m_nMethodType == METHOD_ROKS);

//     this->logger("construct Kohn-Sham matrix F1\n");
//     SymmetricMatrixType F(this->m_nNumOfAOs);
//     VectorType Rho, Myua, Myub;

//     if (this->m_bMemorySave == false) {
//         F.load("fl_Work/fl_Mtr_F1pqtmp");
//     } else {
//         {
//             VectorType Rhoa, Rhob;
//             Rhoa.load(this->getRhoPath(RUN_UKS_ALPHA, this->m_nIteration));
//             Rhob.load(this->getRhoPath(RUN_UKS_BETA,  this->m_nIteration));
//             assert(Rhoa.getSize() == Rhob.getSize());
//             assert(Rhoa.getSize() == static_cast<TlVector::size_type>(this->m_nNumOfAux));

//             Rho = Rhoa + Rhob;
//         }

//         // read Myu
//         Myua.load(this->getMyuPath(RUN_UKS_ALPHA, this->m_nIteration));
//         Myub.load(this->getMyuPath(RUN_UKS_BETA,  this->m_nIteration));
//         assert(Myua.getSize() == Myub.getSize());

//         // read previous Rou
//         VectorType prevRho;
//         if (this->m_nIteration > 1) {
//             VectorType prevRhoa, prevRhob;
//             prevRhoa.load(this->getRhoPath(RUN_UKS_ALPHA, this->m_nIteration -1));
//             prevRhoa.load(this->getRhoPath(RUN_UKS_BETA,  this->m_nIteration -1));
//             assert(prevRhoa.getSize() == prevRhob.getSize());
//             assert(prevRhoa.getSize() == static_cast<TlVector::size_type>(this->m_nNumOfAux));

//             prevRho = prevRhoa + prevRhob;
//         }

//         // read previous Myu
//         VectorType prevMyua, prevMyub;
//         if (this->m_nIteration > 1) {
//             Myua.load(this->getMyuPath(RUN_UKS_ALPHA, this->m_nIteration -1));
//             Myub.load(this->getMyuPath(RUN_UKS_BETA,  this->m_nIteration -1));
//             assert(Myua.getSize() == Myub.getSize());
//             assert(Myua.getSize() == static_cast<TlVector::size_type>(this->m_nNumOfAux));
//         }

//         // calculate delta rou and delta myu
//         {
//             Rho -= prevRho;
//             Myua -= prevMyua;
//             Myub -= prevMyub;

//             // integral object generated then add coulomb contribution of rho to F ?
//             {
//                 DfEriClass dferi(this->pPdfParam_);
//                 dferi.getdeltaHpqA(Rho, F);
//             }

//             // integral object generated then add xc potential contribution of myu-alpha to F ?
//             {
//                 DfOverlapClass dfovr(this->pPdfParam_);
//                 VectorType Myu;
//                 Myu = 0.5 * Myua;
//                 dfovr.getdeltaHpqG(Myu, F);

//                 // integral object generated then add xc potential contribution of myu-beta to F ?
//                 Myu = 0.5 * Myub;
//                 dfovr.getdeltaHpqG(Myu,F);
//             }
//         }
//     }

//     if (this->m_nIteration > 1) {
//         // add previous F-matrix to deltaH
//         SymmetricMatrixType prevFpq;
//         prevFpq.load("fl_Work/fl_Mtr_F1pq.matrix."
//                      + DfObject::m_sRunTypeSuffix[RUN_ROKS]
//                      + TlUtils::xtos(this->m_nIteration -1));
//         assert(prevFpq.getNumOfRows() == this->m_nNumOfAOs);

//         F += prevFpq;

//         // add dummy charge
//         if ((this->m_nNumOfDummyAtoms != 0) && (this->isRestart_ == false) &&
//             (this->m_nIteration <= this->chargeExtrapolateNumber_ +1)) {
//             SymmetricMatrixType Hpq2 = this->getHpq2Matrix<SymmetricMatrixType>();

//             F += Hpq2;
//         }
//     } else {
//         // at first iteration, add one electron part to deltaH
//         SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();

//         if (Hpq.getNumOfRows() != this->m_nNumOfAOs || Hpq.getNumOfCols() != this->m_nNumOfAOs) {
//             CnErr.abort("DfFockmatrix", "mainDirect", "", "program error");
//         }

//         F += Hpq;
//     }

//     F.save("fl_Work/fl_Mtr_F1pq.matrix."
//            + DfObject::m_sRunTypeSuffix[RUN_ROKS]
//            + TlUtils::xtos(this->m_nIteration));

//     // construct kohn-sham matrix of open part ( F2 )
//     this->logger("construct Kohn-Sham matrix F2\n");

//     if (this->m_bMemorySave == false) {
//         F.load("fl_Work/fl_Mtr_F2pqtmp");
//     } else {
//         F = SymmetricMatrixType(this->m_nNumOfAOs);

//         // integral object generated then add coulomb contribution of rho to F ?
//         {
//             VectorType Rhoa = 0.5 * Rho;

//             DfEriClass dferi(this->pPdfParam_);
//             dferi.getdeltaHpqA(Rhoa, F);
//         }

//         // integral object generated then add xc potential contribution of myu-alpha to F ?
//         {
//             VectorType Myu = 0.5 * Myua;

//             DfOverlapClass dfovr(this->pPdfParam_);
//             dfovr.getdeltaHpqG(Myu, F);
//         }
//     }

//     if (this->m_nIteration > 1) {
//         // add previous F-matrix to deltaH
//         SymmetricMatrixType Fpq;
//         Fpq.load("fl_Work/fl_Mtr_F2pq.matrix."
//                  + DfObject::m_sRunTypeSuffix[RUN_ROKS]
//                  + TlUtils::xtos(this->m_nIteration -1));
//         assert(Fpq.getNumOfRows() == this->m_nNumOfAOs);

//         F += Fpq;

//         // add dummy charge
//         if ((this->m_nNumOfDummyAtoms != 0) && (this->isRestart_ == false)
//             && (this->m_nIteration <= this->chargeExtrapolateNumber_ +1)) {
//             SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();

//             F += Hpq2;
//         }

//     } else {
//         // at first iteration, add one electron part to deltaH
//         SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();

//         if (Hpq.getNumOfRows() != this->m_nNumOfAOs || Hpq.getNumOfCols() != this->m_nNumOfAOs) {
//             CnErr.abort("DfFockmatrix", "mainDirect", "", "program error");
//         }

//         Hpq *= 0.5;

//         F += Hpq;
//     }

//     // write F2 matrix to a file
//     F.save("fl_Work/fl_Mtr_F2pq.matrix."
//            + DfObject::m_sRunTypeSuffix[RUN_ROKS]
//            + TlUtils::xtos(this->m_nIteration));

//     // construct united kohn-sham matrix from F1 and F2
//     {
//         // S_1, S_1^dagger と S_2, S_2^dagger を作る
//         {
//             SymmetricMatrixType W = DfObject::getSpqMatrix<SymmetricMatrixType>();

//             SymmetricMatrixType X;
//             X.load("fl_Work/fl_Mtr_P1pq.matrix.roks" + TlUtils::xtos(this->m_nIteration -1));

//             MatrixType A = W * X;
//             A.save("fl_Work/DfFockmatrix.fl_Matrix.S_1");

//             A.transpose();
//             A.save("fl_Work/DfFockmatrix.fl_Matrix.S_1^dagger");

//             SymmetricMatrixType Y;
//             Y.load("fl_Work/fl_Mtr_P2pq.matrix.roks" + TlUtils::xtos(this->m_nIteration -1));

//             A =  W * Y;
//             A.save("fl_Work/DfFockmatrix.fl_Matrix.S_2");

//             A.transpose();
//             A.save("fl_Work/DfFockmatrix.fl_Matrix.S_2^dagger");
//         }

//         // (1 - S_2) F_1 (1 - S_2^dagger) の計算
//         {
//             MatrixType A(this->m_nNumOfAOs, this->m_nNumOfAOs);

//             // (1 - S_2) F_1 の計算
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_2");

//                 SymmetricMatrixType Z;
//                 Z.load("fl_Work/fl_Mtr_F1pq.matrix.rok" + TlUtils::xtos(this->m_nIteration));

//                 B *= -1.0;
//                 for (int i=0; i < this->m_nNumOfAOs; i++) {
//                     B(i, i) += 1.0;
//                 }

//                 A = B * Z;
//             }

//             // "A" (1 - S_2^dagger) の計算
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_2^dagger");

//                 B *= -1.0;
//                 for (int i=0; i < this->m_nNumOfAOs; i++) {
//                     B(i, i) += 1.0;
//                 }

//                 A *= B;
//             }

//             A.save("fl_Work/DfFockmatrix.fl_Matrix.(1-S_2)F_1(1-S_2^dagger)");
//         }

//         // (1 - S_1) F_2 (1 - S_1^dagger) の計算
//         {
//             MatrixType A(this->m_nNumOfAOs, this->m_nNumOfAOs);

//             // (1 - S_1) F_2 の計算
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_1");

//                 SymmetricMatrixType Z;
//                 Z.load("fl_Work/fl_Mtr_F2pq.matrix.roks" + TlUtils::xtos(this->m_nIteration));

//                 B *= -1.0;
//                 for (int i=0; i < this->m_nNumOfAOs; i++) {
//                     B(i, i) += 1.0;
//                 }

//                 A = B * Z;
//             }

//             // "A" (1 - S_1^dagger) の計算
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_1^dagger");

//                 B *= -1.0;
//                 for (int i=0; i < this->m_nNumOfAOs; i++) {
//                     B(i, i) += 1.0;
//                 }

//                 A *= B;
//             }

//             A.save("fl_Work/DfFockmatrix.fl_Matrix.(1-S_1)F_2(1-S_1^dagger)");
//         }


//         // S_1 (F_1 - F_2) S_2^dagger と {S_1 (F_1 - F_2) S_2^dagger}^dagger を作る
//         {

//             // (F_1 - F_2) の計算
//             SymmetricMatrixType Z;
//             Z.load("fl_Work/fl_Mtr_F1pq.matrix.roks" + TlUtils::xtos(this->m_nIteration));
//             {

//                 SymmetricMatrixType Y;
//                 Y.load("fl_Work/fl_Mtr_F2pq.matrix.roks" + TlUtils::xtos(this->m_nIteration));

//                 Z -= Y;
//             }

//             // S_1 "Z" の計算
//             MatrixType A;
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_1");

//                 A = B * Z;
//             }

//             // "A" S_2^dagger の計算
//             {
//                 MatrixType B;
//                 B.load("fl_Work/DfFockmatrix.fl_Matrix.S_2^dagger");

//                 A *= B;
//             }

//             A.save("fl_Work/DfFockmatrix.fl_Matrix.S_1(F_1-F_2)S_2^dagger");

//             A.transpose();
//             A.save("fl_Work/DfFockmatrix.fl_Matrix.(S_1(F_1-F_2)S_2^dagger)^dagger");
//         }

//         // F = "(1 - S_2) F_1 (1 - S_2^dagger)" + "(1 - S_1) F_2 (1 - S_1^dagger)"
//         //   + "S_1 (F_1 - F_2) S_2^dagger" + "{S_1 (F_1 - F_2) S_2^dagger}^dagger" の計算
//         {
//             MatrixType A;
//             A.load("fl_Work/DfFockmatrix.fl_Matrix.(1-S_2)F_1(1-S_2^dagger)");

//             SymmetricMatrixType F(A);

//             A.load("fl_Work/DfFockmatrix.fl_Matrix.(1-S_1)F_2(1-S_1^dagger)");
//             F += A;

//             A.load("fl_Work/DfFockmatrix.fl_Matrix.S_1(F_1-F_2)S_2^dagger");
//             F += A;

//             A.load("fl_Work/DfFockmatrix.fl_Matrix.(S_1(F_1-F_2)S_2^dagger)^dagger");
//             F += A;

//             F.save(this->getFpqMatrixPath(RUN_ROKS, this->m_nIteration));
//         }
//     }
// }

template<typename SymmetricMatrixType, typename VectorType, class DfOverlapClass>
void DfFockMatrix::setXC_RI(const RUN_TYPE nRunType, SymmetricMatrixType& F)
{
    // RI 法
    VectorType Myu = this->getMyu(nRunType, this->m_nIteration);
    if (this->m_nIteration >= 2) {
        const VectorType prevMyu = this->getMyu(nRunType, this->m_nIteration -1);
        Myu -= prevMyu;
    }

    DfOverlapClass dfOverlap(this->pPdfParam_);
    dfOverlap.get_pqg(Myu, &F);
}

template<typename SymmetricMatrixType>
void DfFockMatrix::setXC_DIRECT(const RUN_TYPE nRunType, SymmetricMatrixType& F)
{
    F = this->getFxcMatrix<SymmetricMatrixType>(nRunType, this->m_nIteration);

    if (this->m_nIteration >= 2) {
        // this->setHpq() で前回のFockを足し込むので、
        // ２回転目以降は差分を求める必要がある
        F -= this->getFxcMatrix<SymmetricMatrixType>(nRunType, this->m_nIteration -1);
    }
}

template<typename SymmetricMatrixType, typename VectorType, class DfEriClass>
void DfFockMatrix::setCoulomb(const METHOD_TYPE nMethodType, SymmetricMatrixType& F)
{
    VectorType Rho;
    switch (nMethodType) {
    case METHOD_RKS:
        Rho = this->getRho(RUN_RKS, this->m_nIteration);
        break;
    case METHOD_UKS:
        Rho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration);
        Rho += this->getRho(RUN_UKS_BETA, this->m_nIteration);
        break;
    default:
        std::cerr << "unsopported. sorry." << std::endl;
        CnErr.abort();
        break;
    }

    if (this->m_nIteration >= 2) {
        VectorType prevRho;
        switch (nMethodType) {
        case METHOD_RKS:
            prevRho = this->getRho(RUN_RKS, this->m_nIteration -1);
            break;
        case METHOD_UKS:
            prevRho = this->getRho(RUN_UKS_ALPHA, this->m_nIteration -1);
            prevRho += this->getRho(RUN_UKS_BETA, this->m_nIteration -1);
            break;
        default:
            std::cerr << "unsopported. sorry." << std::endl;
            CnErr.abort();
            break;
        }

        Rho -= prevRho;
    }

    DfEriClass dfEri(this->pPdfParam_);
    dfEri.getJ(Rho, &F);
}


template<typename SymmetricMatrixType>
void DfFockMatrix::setHpq(const RUN_TYPE runType, SymmetricMatrixType& F)
{
    // one electron part =================================================
    if (this->m_nIteration == 1) {
        // at first iteration, add one electron part to deltaH
        {
            SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
            assert(Hpq.getNumOfRows() == this->m_nNumOfAOs);
            F += Hpq;
        }

        if ((this->m_nNumOfDummyAtoms != 0) && (this->chargeExtrapolateNumber_ == 0)) {
            SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();
            F += Hpq2;
        }
    } else {
        // add previous F-matrix to deltaH
        {
            SymmetricMatrixType prevFpq = DfObject::getFpqMatrix<SymmetricMatrixType>(runType, this->m_nIteration -1);
            if (prevFpq.getNumOfRows() != this->m_nNumOfAOs) {
                std::string msg = TlUtils::format("dimension of previous Fock matrix is not equal to the current: %d != %d",
                                                  prevFpq.getNumOfRows(), this->m_nNumOfAOs);
                std::cerr << msg << std::endl;
            }
            assert(prevFpq.getNumOfRows() == this->m_nNumOfAOs);
            F += prevFpq;
        }

        // add dummy charge
        if ((this->m_nNumOfDummyAtoms != 0) && (this->m_nIteration <= this->chargeExtrapolateNumber_ +1)) {
            SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();
            F += Hpq2;
            this->logger(TlUtils::format("Added Dummy Charge / %d\n", this->chargeExtrapolateNumber_));
        }
    }
}

#endif // DFFOCKMATRIX_H
