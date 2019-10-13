#ifndef TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_H
#define TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_H

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

class TlSparseSymmetricMatrix_ViennaCL;
class TlSparseSymmetricMatrix_ImplViennaCL;
class TlDenseGeneralMatrix_ImplViennaCL;
class TlDenseVector_ImplViennaCL;
class TlSparseGeneralMatrix_ImplEigen;

class TlSparseGeneralMatrix_ImplViennaCL : public TlSparseMatrix_ImplObject {
   public:
    typedef viennacl::compressed_matrix<double> MatrixDataType;
#ifdef HAVE_EIGEN
    typedef Eigen::SparseMatrix<double> EigenMatrixDataType;
#endif  // HAVE_EIGEN

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseGeneralMatrix_ImplViennaCL(
        const TlMatrixObject::index_type row = 0,
        const TlMatrixObject::index_type col = 0);
    TlSparseGeneralMatrix_ImplViennaCL(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);
    TlSparseGeneralMatrix_ImplViennaCL(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
    TlSparseGeneralMatrix_ImplViennaCL(const MatrixDataType& rhs);
    TlSparseGeneralMatrix_ImplViennaCL(
        const TlDenseGeneralMatrix_ImplViennaCL& rhs);

#ifdef HAVE_EIGEN
    TlSparseGeneralMatrix_ImplViennaCL(
        const TlSparseGeneralMatrix_ImplEigen& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseGeneralMatrix_ImplViennaCL();

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
    TlSparseGeneralMatrix_ImplViennaCL& operator=(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);

    TlSparseGeneralMatrix_ImplViennaCL& operator+=(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);
    TlSparseGeneralMatrix_ImplViennaCL& operator-=(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);
    TlSparseGeneralMatrix_ImplViennaCL& operator*=(const double coef);

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
    friend TlSparseSymmetricMatrix_ViennaCL;
    friend TlSparseSymmetricMatrix_ImplViennaCL;
    friend TlSparseGeneralMatrix_ImplEigen;

    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_ImplViennaCL operator*(
        const TlSparseGeneralMatrix_ImplViennaCL& sm1,
        const TlSparseGeneralMatrix_ImplViennaCL& sm2);

    // DV = SM(G) * DV
    friend TlDenseVector_ImplViennaCL operator*(
        const TlSparseGeneralMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vtr);
    // DV = DV * SM(G)
    //   friend TlDenseVector_ImplViennaCL operator*(
    //       const TlDenseVector_ImplViennaCL& vtr,
    //       const TlSparseGeneralMatrix_ImplViennaCL& mat);
};

#endif  // TL_SPARSE_GENERAL_MATRIX_IMPL_VIENNACL_H
