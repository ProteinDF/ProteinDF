#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include "tl_dense_matrix_impl_object.h"
#include <viennacl/matrix.hpp>

#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif  // HAVE_EIGEN

class TlDenseSymmetricMatrix_ImplViennaCL;
class TlDenseVector_ImplViennaCL;
class TlDenseGeneralMatrix_ImplEigen;
class TlSparseGeneralMatrix_ImplViennaCL;

class TlDenseGeneralMatrix_ImplViennaCL : public TlDenseMatrix_ImplObject {
 public:
  typedef viennacl::vector<double> VectorDataType;
  typedef viennacl::matrix<double> MatrixDataType;

#ifdef HAVE_EIGEN
  typedef Eigen::MatrixXd EigenMatrixDataType;
  typedef Eigen::Map<MatrixDataType> EigenMapType;
  typedef Eigen::Map<const MatrixDataType> EigenMapTypeConst;
#endif  // HAVE_EIGEN

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_ImplViennaCL(
      const TlMatrixObject::index_type row = 0,
      const TlMatrixObject::index_type col = 0);
  TlDenseGeneralMatrix_ImplViennaCL(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ImplViennaCL(
      const TlDenseSymmetricMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ImplViennaCL(const TlSparseGeneralMatrix_ImplViennaCL& rhs);
#ifdef HAVE_EIGEN
  TlDenseGeneralMatrix_ImplViennaCL(const TlDenseGeneralMatrix_ImplEigen& rhs);
#endif  // HAVE_EIGEN

    void vtr2mat(const std::vector<double>& vtr);

  virtual ~TlDenseGeneralMatrix_ImplViennaCL();

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
  TlDenseGeneralMatrix_ImplViennaCL& operator=(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);

  TlDenseGeneralMatrix_ImplViennaCL& operator+=(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ImplViennaCL& operator-=(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ImplViennaCL& operator*=(const double coef);
  TlDenseGeneralMatrix_ImplViennaCL& operator/=(const double coef);
  TlDenseGeneralMatrix_ImplViennaCL& operator*=(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual double sum() const;
  // virtual double trace() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow,
      TlMatrixObject::index_type* outCol) const;
  virtual void transposeInPlace();

  TlDenseGeneralMatrix_ImplViennaCL& dotInPlace(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ImplViennaCL transpose() const;
  TlDenseGeneralMatrix_ImplViennaCL inverse() const;

  TlDenseGeneralMatrix_ImplViennaCL& reverseColumns();

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  MatrixDataType matrix_;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_ImplViennaCL;
  friend class TlDenseGeneralMatrix_ImplEigen;

  // DM(G) = DM(G) * DM(S)
  friend TlDenseGeneralMatrix_ImplViennaCL operator*(
      const TlDenseGeneralMatrix_ImplViennaCL& mat1,
      const TlDenseSymmetricMatrix_ImplViennaCL& mat2);
  // DM(G) = DM(S) * DM(G)
  friend TlDenseGeneralMatrix_ImplViennaCL operator*(
      const TlDenseSymmetricMatrix_ImplViennaCL& mat1,
      const TlDenseGeneralMatrix_ImplViennaCL& mat2);
  
  // DV = DM(G) * DV
  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseGeneralMatrix_ImplViennaCL& mat,
      const TlDenseVector_ImplViennaCL& vec);
  // DV = DV * DM(G)
  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseVector_ImplViennaCL& vec,
      const TlDenseGeneralMatrix_ImplViennaCL& mat);
};

// DM(G) = double * DM(G)
TlDenseGeneralMatrix_ImplViennaCL operator*(const double coef, const TlDenseGeneralMatrix_ImplViennaCL& dm);

#endif  // TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_H
