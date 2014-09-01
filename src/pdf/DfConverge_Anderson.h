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

#ifndef DFCONVERGE_ANDERSON_H
#define DFCONVERGE_ANDERSON_H

#include "DfConverge_Damping.h"

class DfConverge_Anderson : public DfConverge_Damping {
public:
    DfConverge_Anderson(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Anderson();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

protected:
    template<class VectorType>
    void convergeRhoTilde(const DfObject::RUN_TYPE runType);

    template<class SymmetricMatrixType, class VectorType>
    void convergeKSMatrix(const DfObject::RUN_TYPE runType);

    template<class SymmetricMatrixType, class VectorType>
    void convergePMatrix(const DfObject::RUN_TYPE runType);

protected:
    template<class VectorType>
    VectorType anderson(const VectorType& X1, const VectorType& X2,
                        const VectorType& Y0, const VectorType& Y1);

protected:
    template<class SymmetricMatrixType, class VectorType>
    VectorType getVectorOfKSMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const;

    template<class SymmetricMatrixType, class VectorType>
    void writeKSMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration,
                                 const VectorType& v) const;

    template<class SymmetricMatrixType, class VectorType>
    VectorType getVectorOfPMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const;

    template<class SymmetricMatrixType, class VectorType>
    void writePMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration,
                                const VectorType& v) const;

protected:
    int m_nStartIterationOfAnderson;
    double m_dDampingFactorOfAnderson;
};

template<class VectorType>
void DfConverge_Anderson::convergeRhoTilde(const DfObject::RUN_TYPE runType)
{
    const int nIteration = this->m_nIteration;
    const int nAndersonStartIterationNumber = this->m_nStartIterationOfAnderson;

    this->log_.info("Anderson's convergence method applied to rho~:");
    this->log_.info(TlUtils::format(" start-number = %2d", nAndersonStartIterationNumber));

    if (this->m_nIteration >= (nAndersonStartIterationNumber -1)) {
        VectorType nextY1;
        nextY1.load(this->getRhoPath(runType, nIteration));
        nextY1.save("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));
    }

    if (this->m_nIteration >= nAndersonStartIterationNumber) {
        // do Anderson
        VectorType Y0;
        Y0.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));

        VectorType X1;
        X1.load(this->getRhoPath(runType, nIteration -1));

        VectorType Y1;
        Y1.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration -1));

        VectorType X2;
        X2.load(this->getRhoPath(runType, nIteration -2));

        VectorType X0 = this->anderson<VectorType>(X1, X2, Y0, Y1);
        X0.save(this->getRhoPath(runType, nIteration));
    } else {
        // do damping
        DfConverge_Damping::convergeRhoTilde<VectorType>(runType);
    }
}

template<class SymmetricMatrixType, class VectorType>
void DfConverge_Anderson::convergeKSMatrix(const DfObject::RUN_TYPE runType)
{
    const int nIteration = this->m_nIteration;
    const int nAndersonStartIterationNumber = this->m_nStartIterationOfAnderson;

    this->log_.info("Anderson's convergence method applied to KS Matrix:");
    this->log_.info(TlUtils::format(" start-number = %2d", nAndersonStartIterationNumber));

    if (this->m_nIteration >= (nAndersonStartIterationNumber -1)) {
        VectorType nextY1 = this->getVectorOfKSMatrix<SymmetricMatrixType, VectorType>(runType, nIteration);
        nextY1.save("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));
    }

    if (this->m_nIteration >= nAndersonStartIterationNumber) {
        // do Anderson
        VectorType Y0;
        Y0.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));

        VectorType X1 = this->getVectorOfKSMatrix<SymmetricMatrixType, VectorType>(runType, nIteration -1);

        VectorType Y1;
        Y1.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration -1));

        VectorType X2 = this->getVectorOfKSMatrix<SymmetricMatrixType, VectorType>(runType, nIteration -2);

        VectorType X0 = this->anderson(X1, X2, Y0, Y1);
        this->writeKSMatrixFromVector<SymmetricMatrixType, VectorType>(runType, nIteration, X0);
    } else {
        // do damping
        DfConverge_Damping::convergeKSMatrix<SymmetricMatrixType>(runType);
    }
}

template<class SymmetricMatrixType, class VectorType>
void DfConverge_Anderson::convergePMatrix(const DfObject::RUN_TYPE runType)
{
    // 密度行列はnIteration回目はまだ作成されていない。
    //const int itr = this->m_nIteration -1;
    const int nAndersonStartIterationNumber = this->m_nStartIterationOfAnderson;

    this->log_.info("Anderson's convergence method applied to P Matrix:");
    this->log_.info(TlUtils::format(" start-number = %2d", nAndersonStartIterationNumber));

    if (this->m_nIteration >= (nAndersonStartIterationNumber -1)) {
        VectorType nextY1 = this->getVectorOfPMatrix<SymmetricMatrixType, VectorType>(runType, this->m_nIteration -1);
        nextY1.save("fl_Work/fl_Vct_anderson." + TlUtils::xtos(this->m_nIteration -1));
    }

    if (this->m_nIteration >= nAndersonStartIterationNumber) {
        // do Anderson
        VectorType Y0;
        Y0.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(this->m_nIteration -1));

        VectorType X1 = this->getVectorOfPMatrix<SymmetricMatrixType, VectorType>(runType, this->m_nIteration -2);
        X1.save("fl_Work/fl_Vct_andersonX1." + TlUtils::xtos(this->m_nIteration -2));

        VectorType Y1;
        Y1.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(this->m_nIteration -2));

        VectorType X2 = this->getVectorOfPMatrix<SymmetricMatrixType, VectorType>(runType, this->m_nIteration -3);
        X2.save("fl_Work/fl_Vct_andersonX2." + TlUtils::xtos(this->m_nIteration -3));

        VectorType X0 = this->anderson(X1, X2, Y0, Y1);
        this->writePMatrixFromVector<SymmetricMatrixType, VectorType>(runType, this->m_nIteration -1, X0);
    } else {
        // do damping
        DfConverge_Damping::convergePMatrix<SymmetricMatrixType>(runType);
    }
}



/// Anderson's damping method
///
/// @param[in] X1 (n-1)回目最終的に更新された物理量: X^(n-1)
/// @param[in] X2 (n-2)回目最終的に更新された物理量: X^(n-2)
/// @param[in] Y0 n回目のSCFで決まった物理量: Y^(n)
/// @param[in] Y1 (n-1)回目のSCFで決まった物理量: Y^(n-1)
///
// X0: X^(n)
// X1: X~(n-1)
// 変数末尾の数字の大きさに注意すること
template<class VectorType>
VectorType DfConverge_Anderson::anderson(const VectorType& X1, const VectorType& X2,
                                         const VectorType& Y0, const VectorType& Y1)
{
    const double beta = this->m_dDampingFactorOfAnderson; // damping factor of Anderson's method
    this->log_.info(TlUtils::format(" beta = %f", beta));

    // theta を求める
    double theta = 0.0;
    {
        // r^(n-1) = Y^(n) - X^(n-1)
        const VectorType r1 = Y0 - X1;
        r1.save("fl_Work/fl_Vct_andersonr1." + TlUtils::xtos(this->m_nIteration));

        // r^(n-2) = Y^(n-1) - X^(n-2)
        const VectorType r2 = Y1 - X2;
        r2.save("fl_Work/fl_Vct_andersonr2." + TlUtils::xtos(this->m_nIteration));

        const VectorType r12 = r1 - r2;
        r12.save("fl_Work/fl_Vct_andersonr12." + TlUtils::xtos(this->m_nIteration));

        double t1 = r1 * r12;
        double t2 = r12 * r12;
        this->log_.info(TlUtils::format(" t1 = %f", t1));
        this->log_.info(TlUtils::format(" t2 = %f", t2));
        theta = t1 / t2;
    }
    this->log_.info(TlUtils::format(" theta = %f", theta));

    const double theta_rest = 1.0 - theta;
    // U^(n-1) = (1-theta)*X^(n-1) + theta * X^(n-2)
    //const TlVector U1 = theta_rest * X1 + theta * Y0;
    const VectorType U1 = theta_rest * X1 + theta * X2;

    // V^(n) = (1-theta)*Y^(n) + theta*Y^(n-1)
    const VectorType V0 = theta_rest * Y0 + theta * Y1;

    // X^(n) = (1-theta)*U^(n-1) + theta*V^(n)
    const VectorType X0 = (1.0 - beta) * U1 + beta * V0;

    return X0;
}

template<class SymmetricMatrixType, class VectorType>
VectorType DfConverge_Anderson::getVectorOfKSMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const
{
    SymmetricMatrixType KS;
    KS.load(this->getFpqMatrixPath(runType, nIteration));
    assert(KS.getNumOfRows() == this->m_nNumOfAOs);

    return KS.getVector();
}

template<class SymmetricMatrixType, class VectorType>
void DfConverge_Anderson::writeKSMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration,
                                                  const VectorType& v) const
{
    SymmetricMatrixType KS(v, this->m_nNumOfAOs);
    KS.save(this->getFpqMatrixPath(runType, nIteration));
}

template<class SymmetricMatrixType, class VectorType>
VectorType DfConverge_Anderson::getVectorOfPMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const
{
    SymmetricMatrixType P;
    P.load(this->getPpqMatrixPath(runType, nIteration));
    assert(P.getNumOfRows() == this->m_nNumOfAOs);

    return P.getVector();
}

template<class SymmetricMatrixType, class VectorType>
void DfConverge_Anderson::writePMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration,
                                                 const VectorType& v) const
{
    SymmetricMatrixType P(v, this->m_nNumOfAOs);
    P.save(this->getPpqMatrixPath(runType, nIteration));
}

#endif //DFCONVERGE_ANDERSON_H
