#ifndef TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_FLOAT_H
#define TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_FLOAT_H

// #define VIENNACL_WITH_EIGEN 1
#include <Eigen/Sparse>

#include "tl_sparse_matrix_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_eigen_float.h"
class TlSparseSymmetricMatrix_ImplEigenFloat;
class TlSparseGeneralMatrix_ImplViennaCLFloat;

class TlSparseGeneralMatrix_ImplEigenFloat : public TlSparseMatrix_ImplObject {
public:
    typedef Eigen::SparseMatrix<float> MatrixDataType;
    typedef Eigen::Map<MatrixDataType> MapType;
    typedef Eigen::Map<const MatrixDataType> MapTypeConst;

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseGeneralMatrix_ImplEigenFloat(const TlMatrixObject::index_type row = 0,
                                                  const TlMatrixObject::index_type col = 0);
    TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
    TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& rhs);
    TlSparseGeneralMatrix_ImplEigenFloat(const MatrixDataType& rhs);
    TlSparseGeneralMatrix_ImplEigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs);

#ifdef HAVE_VIENNACL
    TlSparseGeneralMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseGeneralMatrix_ImplEigenFloat();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    virtual TlMatrixObject::index_type getNumOfRows() const;
    virtual TlMatrixObject::index_type getNumOfCols() const;
    virtual void resize(const TlMatrixObject::index_type row,
                        const TlMatrixObject::index_type col);

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
    TlSparseGeneralMatrix_ImplEigenFloat& operator=(
        const TlSparseGeneralMatrix_ImplEigenFloat& rhs);

    TlSparseGeneralMatrix_ImplEigenFloat& operator+=(
        const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
    TlSparseGeneralMatrix_ImplEigenFloat& operator-=(
        const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
    TlSparseGeneralMatrix_ImplEigenFloat& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
    bool load(const std::string& path);
    bool save(const std::string& path) const;

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    mutable MatrixDataType matrix_;

#if __cplusplus >= 201103L
    mutable std::mutex matrix_mutex_;
#endif  // __cplusplus

    static const double reference_;  // a meaningful non zero reference value
    static const double epsilon_;    // machine epsilon

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend TlDenseGeneralMatrix_ImplEigenFloat;
    friend TlSparseSymmetricMatrix_EigenFloat;
    friend TlSparseSymmetricMatrix_ImplEigenFloat;

    friend class TlSparseGeneralMatrix_ImplViennaCLFloat;

    // DM(G) = DM(G) * SM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat2);
    // DM(G) = SM(G) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat2);
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

    // DV = SM(G) * DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vtr);
    // DV = DV * SM(G)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& vtr,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_FLOAT_H
