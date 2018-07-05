#ifndef TL_DENSE_VECTOR_IMPL_SCALAPACK_H
#define TL_DENSE_VECTOR_IMPL_SCALAPACK_H

#include <valarray>
#include <vector>
#include "tl_dense_vector_impl_object.h"

class TlDenseGeneralMatrix_ImplScalapack;
class TlDenseVector_ImplLapack;

class TlDenseVector_ImplScalapack : public TlDenseVector_ImplObject {
 protected:
  // typedef std::valarray<double> DataType;

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseVector_ImplScalapack(
      TlDenseVectorObject::index_type size = 1);
  TlDenseVector_ImplScalapack(const TlDenseVector_ImplScalapack& rhs);
  TlDenseVector_ImplScalapack(const TlDenseVector_ImplLapack& rhs);
  virtual ~TlDenseVector_ImplScalapack();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlDenseVectorObject::size_type getSize() const;
  virtual void resize(const TlDenseVectorObject::index_type newSize);

  virtual double get(const TlDenseVectorObject::index_type i) const;
  virtual void set(const TlDenseVectorObject::index_type i, const double value);
  virtual void add(const TlDenseVectorObject::index_type i, const double value);
  virtual void mul(const TlDenseVectorObject::index_type i, const double value);

  std::vector<double> getVector() const;

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ImplScalapack& operator=(
      const TlDenseVector_ImplScalapack& rhs);

  TlDenseVector_ImplScalapack& operator+=(
      const TlDenseVector_ImplScalapack& rhs);
  TlDenseVector_ImplScalapack& operator-=(
      const TlDenseVector_ImplScalapack& rhs);
  TlDenseVector_ImplScalapack& operator*=(const double rhs);
  TlDenseVector_ImplScalapack& operator/=(const double rhs);

  double operator*(const TlDenseVector_ImplScalapack& rhs) const;

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual void sortByGreater();

  TlDenseVector_ImplScalapack& dotInPlace(
      const TlDenseVector_ImplScalapack& rhs);

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  virtual void initialize();
  virtual int getIndex(int nGlobalRow, int nGlobalCol) const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  int m_nContext;
  int m_pDESC[9];
  int m_nRows;    // 大域行列の行数
  int m_nCols;    // 大域行列の列数
  int m_nMyRows;  // ローカル行列の行数
  int m_nMyCols;  // ローカル行列の列数

  int m_nRank;         // プロセスのランク
  int m_nProc;         // 総プロセス数
  int m_nProcGridRow;  // プロセスグリッドの行数
  int m_nProcGridCol;  // プロセスグリッドの列数
  int m_nMyProcRow;    // プロセスグリッドにおける自分の行数
  int m_nMyProcCol;    // プロセスグリッドにおける自分の列数
  int m_nBlockSize;    // ブロックサイズ

  // index table
  // usage: m_RowIndexTable[local_index] = global_index
  std::vector<int> m_RowIndexTable;
  std::vector<int> m_ColIndexTable;

  double* vector_;

 protected:
  // static int systemBlockSize_;

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_ImplScalapack;
  friend class TlDenseSymmetricMatrix_ImplScalapack;

  friend TlDenseVector_ImplScalapack operator*(
      const TlDenseGeneralMatrix_ImplScalapack& A,
      const TlDenseVector_ImplScalapack& X);
  friend TlDenseVector_ImplScalapack operator*(
      const TlDenseVector_ImplScalapack& X,
      const TlDenseGeneralMatrix_ImplScalapack& A);
};

#endif  // TL_DENSE_VECTOR_IMPLE_SCALAPACK_H
