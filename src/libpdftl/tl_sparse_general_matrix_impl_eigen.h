#ifndef TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_H
#define TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_H

// #define VIENNACL_WITH_EIGEN 1
#include <Eigen/Sparse>
#include "tl_sparse_matrix_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlDenseGeneralMatrix_ImplEigen;
class TlDenseVector_ImplEigen;
class TlSparseSymmetricMatrix_Eigen;
class TlSparseSymmetricMatrix_ImplEigen;
class TlSparseGeneralMatrix_ImplViennaCL;

class TlSparseGeneralMatrix_ImplEigen : public TlSparseMatrix_ImplObject {
 public:
  typedef Eigen::SparseMatrix<double> MatrixDataType;
  typedef Eigen::Map<MatrixDataType> MapType;
  typedef Eigen::Map<const MatrixDataType> MapTypeConst;

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlSparseGeneralMatrix_ImplEigen(
      const TlMatrixObject::index_type row = 0,
      const TlMatrixObject::index_type col = 0);
  TlSparseGeneralMatrix_ImplEigen(const TlSparseGeneralMatrix_ImplEigen& rhs);
  TlSparseGeneralMatrix_ImplEigen(const TlSparseSymmetricMatrix_ImplEigen& rhs);
  TlSparseGeneralMatrix_ImplEigen(const MatrixDataType& rhs);
  TlSparseGeneralMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplEigen& rhs);

  #ifdef HAVE_VIENNACL
  TlSparseGeneralMatrix_ImplEigen(const TlSparseGeneralMatrix_ImplViennaCL& rhs);
  #endif // HAVE_VIENNACL

  virtual ~TlSparseGeneralMatrix_ImplEigen();

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
  TlSparseGeneralMatrix_ImplEigen& operator=(
      const TlSparseGeneralMatrix_ImplEigen& rhs);

  TlSparseGeneralMatrix_ImplEigen& operator+=(
      const TlSparseGeneralMatrix_ImplEigen& rhs);
  TlSparseGeneralMatrix_ImplEigen& operator-=(
      const TlSparseGeneralMatrix_ImplEigen& rhs);
  TlSparseGeneralMatrix_ImplEigen& operator*=(const double coef);

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
  friend TlDenseGeneralMatrix_ImplEigen;
  friend TlSparseSymmetricMatrix_Eigen;
  friend TlSparseSymmetricMatrix_ImplEigen;

  friend class TlSparseGeneralMatrix_ImplViennaCL;

  // DM(G) = DM(G) * SM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat1,
      const TlSparseGeneralMatrix_ImplEigen& mat2);
  // DM(G) = SM(G) * DM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlSparseGeneralMatrix_ImplEigen& mat1,
      const TlDenseGeneralMatrix_ImplEigen& mat2);

  // SM(G) = SM(G) * SM(G)
  friend TlSparseGeneralMatrix_ImplEigen operator*(
      const TlSparseGeneralMatrix_ImplEigen& mat1,
      const TlSparseGeneralMatrix_ImplEigen& mat2);
  // SM(G) = SM(G) * SM(S)
  friend TlSparseGeneralMatrix_ImplEigen operator*(
      const TlSparseGeneralMatrix_ImplEigen& mat1,
      const TlSparseSymmetricMatrix_ImplEigen& mat2);
  // SM(G) = SM(S) * SM(G)
  friend TlSparseGeneralMatrix_ImplEigen operator*(
      const TlSparseSymmetricMatrix_ImplEigen& mat1,
      const TlSparseGeneralMatrix_ImplEigen& mat2);
  // SM(G) = SM(S) * SM(S)
  friend TlSparseGeneralMatrix_ImplEigen operator*(
      const TlSparseSymmetricMatrix_ImplEigen& sm1,
  const TlSparseSymmetricMatrix_ImplEigen& sm2);

  // DV = SM(G) * DV
  friend TlDenseVector_ImplEigen operator*(
      const TlSparseGeneralMatrix_ImplEigen& mat,
      const TlDenseVector_ImplEigen& vtr);

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_SPARSE_GENERAL_MATRIX_IMPL_EIGEN_H
