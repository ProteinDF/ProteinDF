#ifndef TL_DENSE_GENERAL_MATRIX_SCALAPACK_H
#define TL_DENSE_GENERAL_MATRIX_SCALAPACK_H

#include "tl_dense_general_matrix_impl_scalapack.h"
#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_Scalapack;
class TlDenseVector_Scalapack;
class TlSparseMatrix;

class TlDenseGeneralMatrix_Scalapack : public TlDenseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_Scalapack(
      const TlMatrixObject::index_type row = 1,
      const TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_Scalapack(const TlDenseGeneralMatrix_Scalapack& rhs);
  TlDenseGeneralMatrix_Scalapack(const TlDenseSymmetricMatrix_Scalapack& rhs);
  TlDenseGeneralMatrix_Scalapack(const TlDenseGeneralMatrix_ImplScalapack& rhs);

  virtual ~TlDenseGeneralMatrix_Scalapack();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlMatrixObject::size_type getNumOfElements() const;
  virtual double getLocal(TlMatrixObject::index_type row,
                          TlMatrixObject::index_type col) const;

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_Scalapack& operator=(
      const TlDenseGeneralMatrix_Scalapack& rhs);

  const TlDenseGeneralMatrix_Scalapack operator+(
      const TlDenseGeneralMatrix_Scalapack& rhs) const;
  const TlDenseGeneralMatrix_Scalapack operator-(
      const TlDenseGeneralMatrix_Scalapack& rhs) const;
  const TlDenseGeneralMatrix_Scalapack operator*(const double coef) const;
  const TlDenseGeneralMatrix_Scalapack operator*(
      const TlDenseGeneralMatrix_Scalapack& rhs) const;

  TlDenseGeneralMatrix_Scalapack& operator+=(
      const TlDenseGeneralMatrix_Scalapack& rhs);
  TlDenseGeneralMatrix_Scalapack& operator-=(
      const TlDenseGeneralMatrix_Scalapack& rhs);
  TlDenseGeneralMatrix_Scalapack& operator*=(const double coef);
  TlDenseGeneralMatrix_Scalapack& operator/=(const double coef);
  TlDenseGeneralMatrix_Scalapack& operator*=(
      const TlDenseGeneralMatrix_Scalapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double sum() const;
  double getRMS() const;

  TlDenseGeneralMatrix_Scalapack transpose() const;

  TlDenseGeneralMatrix_Scalapack dot(const TlDenseGeneralMatrix_Scalapack& rhs) const;
  const TlDenseGeneralMatrix_Scalapack& dotInPlace(
      const TlDenseGeneralMatrix_Scalapack& rhs);

  TlDenseGeneralMatrix_Scalapack inverse() const;

  bool getSparseMatrix(TlSparseMatrix* pMatrix, bool isFinalize = false) const;
  void mergeSparseMatrix(const TlSparseMatrix& M);

  std::vector<TlMatrixObject::index_type> getRowIndexTable() const;
  std::vector<TlMatrixObject::index_type> getColIndexTable() const;
  void getLocalMatrix(TlDenseGeneralMatrixObject* pOutputMatrix) const;

 public:
  // double* data();
  // const double* data() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  static void setUsingPartialIO(bool isUsePIO);

  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Scalapack;

  // matrix(sym) x matrix(gen)
  friend TlDenseGeneralMatrix_Scalapack operator*(
      const TlDenseSymmetricMatrix_Scalapack& rhs1,
      const TlDenseGeneralMatrix_Scalapack& rhs2);
  // matrix(gen) x matrix(sym)
  friend TlDenseGeneralMatrix_Scalapack operator*(
      const TlDenseGeneralMatrix_Scalapack& rhs1,
      const TlDenseSymmetricMatrix_Scalapack& rhs2);

  // matrix x vector
  friend TlDenseVector_Scalapack operator*(
      const TlDenseGeneralMatrix_Scalapack& rhs1,
      const TlDenseVector_Scalapack& rhs2);
  friend TlDenseVector_Scalapack operator*(
      const TlDenseVector_Scalapack& rhs1,
      const TlDenseGeneralMatrix_Scalapack& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_SCALAPACK_H
