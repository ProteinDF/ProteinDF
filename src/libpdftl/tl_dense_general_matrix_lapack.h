#ifndef TL_DENSE_GENERAL_MATRIX_LAPACK_H
#define TL_DENSE_GENERAL_MATRIX_LAPACK_H

#include "tl_dense_general_matrix_impl_lapack.h"
#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_Lapack;
class TlDenseVector_Lapack;

class TlDenseGeneralMatrix_Lapack : public TlDenseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_Lapack(
      const TlMatrixObject::index_type row = 1,
      const TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack(const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_ImplLapack& rhs);

  void vtr2mat(const std::vector<double>& vtr);

  virtual ~TlDenseGeneralMatrix_Lapack();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlMatrixObject::size_type getNumOfElements() const;

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_Lapack& operator=(
      const TlDenseGeneralMatrix_Lapack& rhs);

  const TlDenseGeneralMatrix_Lapack operator+(
      const TlDenseGeneralMatrix_Lapack& rhs) const;
  const TlDenseGeneralMatrix_Lapack operator-(
      const TlDenseGeneralMatrix_Lapack& rhs) const;
  const TlDenseGeneralMatrix_Lapack operator*(const double coef) const;
  const TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseGeneralMatrix_Lapack& rhs) const;

  TlDenseGeneralMatrix_Lapack& operator+=(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack& operator-=(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack& operator*=(const double coef);
  TlDenseGeneralMatrix_Lapack& operator/=(const double coef);
  TlDenseGeneralMatrix_Lapack& operator*=(
      const TlDenseGeneralMatrix_Lapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double sum() const;
  double getRMS() const;

  const TlDenseGeneralMatrix_Lapack& dotInPlace(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack transpose() const;
  TlDenseGeneralMatrix_Lapack inverse() const;

  // solve Ax = B
  TlDenseGeneralMatrix_Lapack getLeastSquaresSolution(
      const TlDenseGeneralMatrix_Lapack& B) const;

 public:
  double* data();
  const double* data() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  void dump(TlDenseVector_Lapack* v) const;
  void restore(const TlDenseVector_Lapack& v);

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Lapack;

  // matrix(sym) x matrix(gen)
  friend TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseSymmetricMatrix_Lapack& rhs1,
      const TlDenseGeneralMatrix_Lapack& rhs2);
  // matrix(gen) x matrix(sym)
  friend TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseGeneralMatrix_Lapack& rhs1,
      const TlDenseSymmetricMatrix_Lapack& rhs2);

  friend TlDenseVector_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1,
                                        const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator*(
      const TlDenseVector_Lapack& rhs1,
      const TlDenseGeneralMatrix_Lapack& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_LAPACK_H
