#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_eigen.h"

class TlDenseSymmetricMatrix_ImplEigenFloat;
class TlDenseVector_ImplEigen;
class TlSparseSymmetricMatrix_ImplEigen;
class TlDenseSymmetricMatrix_ImplViennaCL;

class TlDenseSymmetricMatrix_ImplEigen : public TlDenseGeneralMatrix_ImplEigen {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_ImplEigen(const TlMatrixObject::index_type dim = 0,
                                              double const* const pBuf = NULL);
    TlDenseSymmetricMatrix_ImplEigen(const TlDenseSymmetricMatrix_ImplEigen& rhs);
    TlDenseSymmetricMatrix_ImplEigen(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs);
    TlDenseSymmetricMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplEigen& rhs);
    TlDenseSymmetricMatrix_ImplEigen(const TlSparseSymmetricMatrix_ImplEigen& sm);

#ifdef HAVE_VIENNACL
    TlDenseSymmetricMatrix_ImplEigen(const TlDenseSymmetricMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlDenseSymmetricMatrix_ImplEigen();

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
    TlDenseSymmetricMatrix_ImplEigen transpose() const;
    TlDenseSymmetricMatrix_ImplEigen inverse() const;
    bool eig(TlDenseVector_ImplEigen* pEigVal,
             TlDenseGeneralMatrix_ImplEigen* pEigVec) const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    virtual void vtr2mat(double const* const pBuf);

protected:
    virtual TlMatrixObject::size_type getNumOfElements() const;
    double* data();
    const double* data() const;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_Eigen;

    friend class TlDenseGeneralMatrix_ImplEigen;
    friend class TlSparseSymmetricMatrix_ImplEigen;

    friend class TlDenseSymmetricMatrix_ImplViennaCL;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ImplEigen operator*(
        const TlDenseGeneralMatrix_ImplEigen& mat1,
        const TlDenseSymmetricMatrix_ImplEigen& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigen operator*(
        const TlDenseSymmetricMatrix_ImplEigen& mat1,
        const TlDenseGeneralMatrix_ImplEigen& mat2);

    // DM(S) * DV
    friend TlDenseVector_ImplEigen operator*(
        const TlDenseSymmetricMatrix_ImplEigen& dms,
        const TlDenseVector_ImplEigen& dv);
    // DV * DM(S)
    friend TlDenseVector_ImplEigen operator*(
        const TlDenseVector_ImplEigen& dv,
        const TlDenseSymmetricMatrix_ImplEigen& dms);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

TlDenseSymmetricMatrix_ImplEigen operator*(
    const double coef, const TlDenseSymmetricMatrix_ImplEigen& DM);
TlDenseSymmetricMatrix_ImplEigen operator*(
    const TlDenseSymmetricMatrix_ImplEigen& DM, const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
