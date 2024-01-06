#ifndef DF_CONVERGE_DAMPING_ANDERSON_TEMPLATE_H
#define DF_CONVERGE_DAMPING_ANDERSON_TEMPLATE_H

#include "df_converge_damping_template.h"

template <class Matrix, class SymmetricMatrix, class Vector>
class DfConverge_Damping_Anderson_Template : public DfConverge_Damping_Template<SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_Anderson_Template(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Anderson_Template();

protected:
    virtual void convergeRhoTilde(const DfObject::RUN_TYPE runType);
    virtual void convergeKSMatrix(const DfObject::RUN_TYPE runType);
    virtual void convergePMatrix(const DfObject::RUN_TYPE runType);

protected:
    Vector anderson(const Vector& X1, const Vector& X2, const Vector& Y0, const Vector& Y1);

protected:
    Vector getVectorOfKSMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const;
    void writeKSMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration, const Vector& v) const;

    Vector getVectorOfSymmetricMatrix(const std::string& path) const;
    void saveSymmetricMatrixFromVector(const Vector& v, const std::string& path) const;

protected:
    int startIterationOfAnderson_;
    double dampingFactorOfAnderson_;
};

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::DfConverge_Damping_Anderson_Template(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Template<SymmetricMatrix, Vector>(pPdfParam) {
}

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::~DfConverge_Damping_Anderson_Template() {
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    const int nIteration = this->m_nIteration;
    const int startIterationOfAnderson = this->startIterationOfAnderson_;

    this->log_.info("Anderson's convergence method applied to rho~:");
    this->log_.info(
        TlUtils::format(" start-number = %2d", startIterationOfAnderson));

    if (this->m_nIteration >= (startIterationOfAnderson - 1)) {
        Vector nextY1;
        nextY1.load(this->getRhoPath(runType, nIteration));
        nextY1.save("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));
    }

    if (this->m_nIteration >= startIterationOfAnderson) {
        // do Anderson
        Vector Y0;
        Y0.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));

        Vector X1;
        X1.load(this->getRhoPath(runType, nIteration - 1));

        Vector Y1;
        Y1.load("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration - 1));

        Vector X2;
        X2.load(this->getRhoPath(runType, nIteration - 2));

        Vector X0 = this->anderson(X1, X2, Y0, Y1);
        X0.save(this->getRhoPath(runType, nIteration));
    } else {
        // do damping
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    const int nIteration = this->m_nIteration;
    const int startIterationOfAnderson = this->startIterationOfAnderson_;

    this->log_.info("Anderson's convergence method applied to KS Matrix:");
    this->log_.info(TlUtils::format(" start-number = %2d", startIterationOfAnderson));

    if (this->m_nIteration >= (startIterationOfAnderson - 1)) {
        Vector nextY1 = this->getVectorOfKSMatrix(runType, nIteration);
        nextY1.save("fl_Work/fl_Vct_anderson." + TlUtils::xtos(nIteration));
    }

    if (this->m_nIteration >= startIterationOfAnderson) {
        // do Anderson
        const Vector X1 = this->getVectorOfSymmetricMatrix(DfObject::getFpqMatrixPath(runType, this->m_nIteration - 1));
        const Vector Y1 = this->getVectorOfSymmetricMatrix(DfObject::getFpqMatrixPath(runType, this->m_nIteration - 1));

        const Vector X2 = this->getVectorOfSymmetricMatrix(DfObject::getFpqMatrixPath(runType, this->m_nIteration - 2));
        const Vector Y2 = this->getVectorOfSymmetricMatrix(DfObject::getFpqMatrixPath(runType, this->m_nIteration - 2));

        const Vector X0 = this->anderson(X1, X2, Y1, Y2);
        this->saveSymmetricMatrixFromVector(X0, DfObject::getFpqMatrixPath(runType, this->m_nIteration));
    } else {
        // do damping
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::convergePMatrix(const DfObject::RUN_TYPE runType) {
    // 密度行列はnIteration回目はまだ作成されていない。
    // const int itr = this->m_nIteration -1;
    const int startIterationOfAnderson = this->startIterationOfAnderson_;

    this->log_.info("Anderson's convergence method applied to P Matrix:");
    this->log_.info(TlUtils::format(" start-number = %2d", startIterationOfAnderson));

    if (this->m_nIteration >= startIterationOfAnderson) {
        // do Anderson
        const Vector X1 = this->getVectorOfSymmetricMatrix(DfObject::getPInMatrixPath(runType, this->m_nIteration - 1));
        const Vector Y1 = this->getVectorOfSymmetricMatrix(DfObject::getPOutMatrixPath(runType, this->m_nIteration - 1));

        const Vector X2 = this->getVectorOfSymmetricMatrix(DfObject::getPInMatrixPath(runType, this->m_nIteration - 2));
        const Vector Y2 = this->getVectorOfSymmetricMatrix(DfObject::getPOutMatrixPath(runType, this->m_nIteration - 2));

        const Vector X0 = this->anderson(X1, X2, Y1, Y2);
        this->saveSymmetricMatrixFromVector(X0, DfObject::getPInMatrixPath(runType, this->m_nIteration));
    } else {
        // do damping
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix(runType);
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
template <class Matrix, class SymmetricMatrix, class Vector>
Vector DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::anderson(const Vector& X1,
                                                                                       const Vector& X2,
                                                                                       const Vector& Y1,
                                                                                       const Vector& Y2) {
    const double beta = this->dampingFactorOfAnderson_;  // damping factor of Anderson's method
    this->log_.info(TlUtils::format(" beta = %f", beta));

    // theta を求める
    double theta = 0.0;
    {
        // r^(n-1) = Y^(n) - X^(n-1)
        const Vector r1 = Y1 - X1;
        // r1.save("fl_Work/fl_Vct_andersonr1." + TlUtils::xtos(this->m_nIteration));

        // r^(n-2) = Y^(n-1) - X^(n-2)
        const Vector r2 = Y2 - X2;
        // r2.save("fl_Work/fl_Vct_andersonr2." + TlUtils::xtos(this->m_nIteration));

        const Vector r12 = r1 - r2;
        // r12.save("fl_Work/fl_Vct_andersonr12." + TlUtils::xtos(this->m_nIteration));

        const double t1 = r1 * r12;
        const double t2 = r12 * r12;
        this->log_.info(TlUtils::format(" t1 = %f", t1));
        this->log_.info(TlUtils::format(" t2 = %f", t2));
        theta = t1 / t2;
    }
    this->log_.info(TlUtils::format(" theta = %f", theta));

    const double theta_rest = 1.0 - theta;
    // U^(n-1) = (1-theta)*X^(n-1) + theta * X^(n-2)
    // const TlDenseVector_Lapack U1 = theta_rest * X1 + theta * Y1;
    const Vector U1 = theta_rest * X1 + theta * X2;

    // V^(n) = (1-theta)*Y^(n) + theta*Y^(n-1)
    const Vector V1 = theta_rest * Y1 + theta * Y2;

    // X^(n) = (1-theta)*U^(n-1) + theta*V^(n)
    const Vector X0 = (1.0 - beta) * U1 + beta * V1;

    return X0;
}

template <class Matrix, class SymmetricMatrix, class Vector>
Vector DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::getVectorOfKSMatrix(const DfObject::RUN_TYPE runType, const int nIteration) const {
    SymmetricMatrix KS;
    KS.load(this->getFpqMatrixPath(runType, nIteration));
    assert(KS.getNumOfRows() == this->m_nNumOfAOs);

    Vector v;
    KS.dump(&v);

    return v;
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::writeKSMatrixFromVector(const DfObject::RUN_TYPE runType, const int nIteration, const Vector& v) const {
    SymmetricMatrix KS(this->m_nNumOfAOs);
    KS.restore(v);

    KS.save(this->getFpqMatrixPath(runType, nIteration));
}

template <class Matrix, class SymmetricMatrix, class Vector>
Vector DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::getVectorOfSymmetricMatrix(const std::string& path) const {
    SymmetricMatrix M;
    M.load(path);
    assert(M.getNumOfRows() == this->m_nNumOfAOs);
    assert(M.getNumOfCols() == this->m_nNumOfAOs);

    Vector v;
    M.dump(&v);

    return v;
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Anderson_Template<Matrix, SymmetricMatrix, Vector>::saveSymmetricMatrixFromVector(const Vector& v, const std::string& path) const {
    SymmetricMatrix M(this->m_nNumOfAOs);
    M.restore(v);

    M.save(path);
}

#endif  // DF_CONVERGE_DAMPING_ANDERSON_TEMPLATE_H
