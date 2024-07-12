#ifndef TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H
#define TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H

// #include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"

class TlSparseSymmetricMatrix_ImplViennaCL;

class TlSparseSymmetricMatrix_ImplEigenFloat : public TlSparseGeneralMatrix_ImplEigenFloat {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseSymmetricMatrix_ImplEigenFloat(const TlMatrixObject::index_type dim = 0);
    TlSparseSymmetricMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& rhs);
    TlSparseSymmetricMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
    TlSparseSymmetricMatrix_ImplEigenFloat(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs);

#ifdef HAVE_VIENNACL
    TlSparseSymmetricMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseSymmetricMatrix_ImplEigenFloat();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    virtual double get(const TlMatrixObject::index_type row,
                       const TlMatrixObject::index_type col) const;

    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void mul(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
    TlSparseSymmetricMatrix_ImplEigenFloat& operator+=(const TlSparseSymmetricMatrix_ImplEigenFloat& sm);
    TlSparseSymmetricMatrix_ImplEigenFloat& operator-=(const TlSparseSymmetricMatrix_ImplEigenFloat& sm);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ImplEigenFloat;
    friend class TlSparseSymmetricMatrix_ImplViennaCL;

    // DM(G) = DM(G) * SM(S)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat2);
    // DM(G) = SM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat1,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& sm1,
        const TlSparseSymmetricMatrix_ImplEigenFloat& sm2);

    friend TlDenseVector_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vtr);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_FLOAT_H
