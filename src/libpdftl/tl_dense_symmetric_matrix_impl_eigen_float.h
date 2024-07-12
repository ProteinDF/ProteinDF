#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"

class TlSparseSymmetricMatrix_ImplEigen;
class TlDenseSymmetricMatrix_ImplViennaCL;

class TlDenseSymmetricMatrix_ImplEigenFloat : public TlDenseGeneralMatrix_ImplEigenFloat {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_ImplEigenFloat(const TlMatrixObject::index_type dim = 0,
                                                   double const* const pBuf = NULL);
    TlDenseSymmetricMatrix_ImplEigenFloat(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs);
    TlDenseSymmetricMatrix_ImplEigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseSymmetricMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& sm);

#ifdef HAVE_VIENNACL
    TlDenseSymmetricMatrix_ImplEigenFloat(
        const TlDenseSymmetricMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlDenseSymmetricMatrix_ImplEigenFloat();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    // virtual double sum() const;
    // virtual double getRMS() const;
    // virtual double getMaxAbsoluteElement(
    //     TlMatrixObject::index_type* outRow,
    //     TlMatrixObject::index_type* outCol) const;
    //
    // const TlDenseSymmetricMatrixMatrix_Eigen_new& dotInPlace(
    //     const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_ImplEigenFloat transpose() const;
    TlDenseSymmetricMatrix_ImplEigenFloat inverse() const;
    bool eig(TlDenseVector_ImplEigenFloat* pEigVal,
             TlDenseGeneralMatrix_ImplEigenFloat* pEigVec) const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    virtual void vtr2mat(double const* const pBuf);

protected:
    virtual TlMatrixObject::size_type getNumOfElements() const;
    float* data();
    const float* data() const;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_EigenFloat;

    friend class TlDenseGeneralMatrix_ImplEigenFloat;
    friend class TlSparseSymmetricMatrix_ImplEigenFloat;

    friend class TlDenseSymmetricMatrix_ImplViennaCL;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

    // DM(S) * DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseSymmetricMatrix_ImplEigenFloat& dms,
        const TlDenseVector_ImplEigenFloat& dv);
    // DV * DM(S)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& dv,
        const TlDenseSymmetricMatrix_ImplEigenFloat& dms);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

TlDenseSymmetricMatrix_ImplEigenFloat operator*(
    const double coef, const TlDenseSymmetricMatrix_ImplEigenFloat& DM);
TlDenseSymmetricMatrix_ImplEigenFloat operator*(
    const TlDenseSymmetricMatrix_ImplEigenFloat& DM, const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H
