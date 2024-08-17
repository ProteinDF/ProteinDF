#ifndef TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H
#define TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_EIGEN
#define VIENNACL_WITH_EIGEN 1
#include <Eigen/Sparse>
#endif  // HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>

#include "tl_sparse_matrix_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlSparseSymmetricMatrix_ViennaCLFloat;
class TlSparseSymmetricMatrix_ImplViennaCLFloat;
class TlDenseGeneralMatrix_ImplViennaCLFloat;
class TlDenseVector_ImplViennaCLFloat;
class TlSparseGeneralMatrix_ImplEigenFloat;

class TlSparseGeneralMatrix_ImplViennaCLFloat : public TlSparseMatrix_ImplObject {
public:
    typedef viennacl::compressed_matrix<float> MatrixDataType;
#ifdef HAVE_EIGEN
    typedef Eigen::SparseMatrix<float> EigenMatrixDataType;
#endif  // HAVE_EIGEN

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseGeneralMatrix_ImplViennaCLFloat(
        const TlMatrixObject::index_type row = 0,
        const TlMatrixObject::index_type col = 0);
    TlSparseGeneralMatrix_ImplViennaCLFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ImplViennaCLFloat(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ImplViennaCLFloat(const MatrixDataType& rhs);
    TlSparseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlSparseGeneralMatrix_ImplViennaCLFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseGeneralMatrix_ImplViennaCLFloat();

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
    TlSparseGeneralMatrix_ImplViennaCLFloat& operator=(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);

    TlSparseGeneralMatrix_ImplViennaCLFloat& operator+=(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ImplViennaCLFloat& operator-=(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ImplViennaCLFloat& operator*=(const float coef);

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
    //   bool load(const std::string& path);
    //   bool save(const std::string& path) const;

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    mutable MatrixDataType matrix_;

#if __cplusplus >= 201103L
    mutable std::mutex matrix_mutex_;
#endif  // __cplusplus

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend TlSparseSymmetricMatrix_ViennaCLFloat;
    friend TlSparseSymmetricMatrix_ImplViennaCLFloat;
    friend TlSparseGeneralMatrix_ImplEigenFloat;

    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& sm1,
        const TlSparseGeneralMatrix_ImplViennaCLFloat& sm2);

    // DV = SM(G) * DV
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vtr);
    // DV = DV * SM(G)
    //   friend TlDenseVector_ImplViennaCLFloat operator*(
    //       const TlDenseVector_ImplViennaCLFloat& vtr,
    //       const TlSparseGeneralMatrix_ImplViennaCLFloat& mat);
};

#endif  // TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H
