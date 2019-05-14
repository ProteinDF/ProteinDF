#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_EIGEN_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_EIGEN_H

#include <Eigen/Core>
#include "tl_dense_matrix_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlDenseGeneralMatrix_ImplEigen;
class TlDenseSymmetricMatrix_ImplEigen;
class TlDenseVector_ImplEigen;
class TlSparseGeneralMatrix_ImplEigen;
class TlSparseSymmetricMatrix_ImplEigen;
class TlDenseGeneralMatrix_ImplViennaCL;

class TlDenseGeneralMatrix_ImplEigen : public TlDenseMatrix_ImplObject {
 public:
  typedef Eigen::MatrixXd MatrixDataType;
  typedef Eigen::Map<MatrixDataType> MapType;
  typedef Eigen::Map<const MatrixDataType> MapTypeConst;

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_ImplEigen(
      const TlMatrixObject::index_type row = 0,
      const TlMatrixObject::index_type col = 0);
  TlDenseGeneralMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen(const TlDenseSymmetricMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen(const MatrixDataType& rhs);
  TlDenseGeneralMatrix_ImplEigen(const TlSparseGeneralMatrix_ImplEigen& sm);

#ifdef HAVE_VIENNACL
  TlDenseGeneralMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

  void vtr2mat(const std::vector<double>& vtr);

  virtual ~TlDenseGeneralMatrix_ImplEigen();

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

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_ImplEigen& operator=(
      const TlDenseGeneralMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen& operator+=(
      const TlDenseGeneralMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen& operator-=(
      const TlDenseGeneralMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen& operator*=(const double coef);
  TlDenseGeneralMatrix_ImplEigen& operator/=(const double coef);
  TlDenseGeneralMatrix_ImplEigen& operator*=(
      const TlDenseGeneralMatrix_ImplEigen& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual double sum() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow,
      TlMatrixObject::index_type* outCol) const;
  virtual void transposeInPlace();

  const TlDenseGeneralMatrix_ImplEigen& dotInPlace(
      const TlDenseGeneralMatrix_ImplEigen& rhs);
  TlDenseGeneralMatrix_ImplEigen transpose() const;
  TlDenseGeneralMatrix_ImplEigen inverse() const;

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------

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
  friend class TlDenseSymmetricMatrix_ImplEigen;
  friend class TlSparseGeneralMatrix_ImplEigen;
  friend class TlDenseGeneralMatrix_ImplViennaCL;

  // DM(G) = DM(G) * DM(S)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat1,
      const TlDenseSymmetricMatrix_ImplEigen& mat2);
  // DM(G) = DM(S) * DM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlDenseSymmetricMatrix_ImplEigen& mat1,
      const TlDenseGeneralMatrix_ImplEigen& mat2);

  // DV = DM(G) * DV
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat,
      const TlDenseVector_ImplEigen& vec);
  // DV = DV * DM(G)
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseVector_ImplEigen& vec,
      const TlDenseGeneralMatrix_ImplEigen& mat);

  // v.s. sparse matrix
  // DM(G) = DM(G) * SM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat1,
      const TlSparseGeneralMatrix_ImplEigen& mat2);
  // DM(G) = SM(G) * DM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlSparseGeneralMatrix_ImplEigen& mat1,
      const TlDenseGeneralMatrix_ImplEigen& mat2);

  // DM(G) = DM(G) * SM(S)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat1,
      const TlSparseSymmetricMatrix_ImplEigen& mat2);
  // DM(G) = SM(S) * DM(G)
  friend TlDenseGeneralMatrix_ImplEigen operator*(
      const TlSparseSymmetricMatrix_ImplEigen& mat1,
      const TlDenseGeneralMatrix_ImplEigen& mat2);

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

TlDenseGeneralMatrix_ImplEigen operator*(
    const double coef, const TlDenseGeneralMatrix_ImplEigen& DM);
TlDenseGeneralMatrix_ImplEigen operator*(
    const TlDenseGeneralMatrix_ImplEigen& DM, const double coef);

#endif  // TL_DENSE_GENERAL_MATRIX_IMPLE_EIGEN_H
