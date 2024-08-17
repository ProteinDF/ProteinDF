#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_FLOAT_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_FLOAT_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include "tl_dense_general_matrix_impl_viennacl_float.h"

class TlDenseGeneralMatrix_ImplViennaCLFloat;
class TlDenseVector_ImplViennaCLFloat;
class TlDenseSymmetricMatrix_ImplEigenFloat;
class TlSparseSymmetricMatrix_ImplViennaCLFloat;

class TlDenseSymmetricMatrix_ImplViennaCLFloat : public TlDenseGeneralMatrix_ImplViennaCLFloat {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_ImplViennaCLFloat(const TlMatrixObject::index_type dim = 0);
    TlDenseSymmetricMatrix_ImplViennaCLFloat(const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ImplViennaCLFloat(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlDenseSymmetricMatrix_ImplViennaCLFloat(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs);
#endif  // HAVE_EIGEN
    virtual ~TlDenseSymmetricMatrix_ImplViennaCLFloat();

    virtual void vtr2mat(const std::vector<double>& vtr);

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const float value);

    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const float value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    // virtual float sum() const;
    // virtual float getRMS() const;
    // virtual float getMaxAbsoluteElement(TlMatrixObject::index_type* outRow,
    //                              TlMatrixObject::index_type* outCol) const;
    //
    // TlDenseSymmetetricMatrix_ImplViennaCL& dotInPlace(
    //     const TlDenseSymmetetricMatrix_ImplViennaCL& rhs);
    TlDenseSymmetricMatrix_ImplViennaCLFloat transpose() const;
    TlDenseSymmetricMatrix_ImplViennaCLFloat inverse() const;

    bool eig(TlDenseVector_ImplViennaCLFloat* pEigval,
             TlDenseGeneralMatrix_ImplViennaCLFloat* pEigvec) const;
    // bool eig_powerIteration(TlDenseVector_ImplViennaCLFloat* pEigval,
    //                         TlDenseGeneralMatrix_ImplViennaCLFloat* pEigvec) const;
    bool eig_QR(TlDenseVector_ImplViennaCLFloat* pEigval,
                TlDenseGeneralMatrix_ImplViennaCLFloat* pEigvec) const;

    // ---------------------------------------------------------------------------
    // private
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_ImplViennaCLFloat;
    friend class TlDenseSymmetricMatrix_ImplEigenFloat;

    friend TlDenseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat1,
        const TlDenseSymmetricMatrix_ImplViennaCLFloat& mat2);
    friend TlDenseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlDenseSymmetricMatrix_ImplViennaCLFloat& mat1,
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat2);
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_FLOAT_H
