#include "tl_dense_vector_impl_scalapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include "TlCommunicate.h"
#include "scalapack.h"
#include "tl_scalapack_context.h"

// int TlDenseVector_ImplScalapack::systemBlockSize_ = 64;

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVectorObject::index_type nSize)
    : m_nContext(0),
      m_nRows(nSize),
      m_nCols(1),
      m_nBlockSize(TlScalapackContext::getBlockSize()),
      vector_(NULL) {
  this->initialize();
}

TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVector_ImplScalapack& rhs)
    : m_nContext(0),
      m_nRows(rhs.m_nRows),
      m_nCols(rhs.m_nCols),
      m_nBlockSize(rhs.m_nBlockSize),
      vector_(NULL) {
  this->initialize();
  this->vector_ = rhs.vector_;
}

// TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
//     const std::vector<double>& rhs,
//     const TlDenseVectorObject::size_type globalSize)
//     : m_nContext(0),
//       m_nRows(globalSize),
//       m_nCols(1),
//       m_nBlockSize(TlDenseVector_ImplScalapack::systemBlockSize_),
//       vector_(NULL) {
//   this->initialize();
//   const size_t copySize = std::min(rhs.size(), this->vector_.size());
//   std::memcpy(&(this->vector_[0]), &(rhs[0]), copySize);
// }

TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVector_ImplLapack& rhs)
    : m_nContext(0),
      m_nRows(rhs.getSize()),
      m_nCols(1),
      m_nBlockSize(TlScalapackContext::getBlockSize()),
      vector_(NULL) {
  this->initialize();

  TlDenseVectorObject::size_type size = rhs.getSize();
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    const int index = this->getIndex(i, 0);
    if (index != -1) {
      this->vector_[index] = rhs.get(i);
    }
  }
}

TlDenseVector_ImplScalapack::~TlDenseVector_ImplScalapack() {
  if (this->vector_ != NULL) {
    delete[] this->vector_;
    this->vector_ = NULL;
  }
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplScalapack::getSize() const {
  return this->m_nRows;
}

void TlDenseVector_ImplScalapack::resize(
    const TlDenseVectorObject::index_type size) {
  assert(size > 0);

  TlDenseVector_ImplScalapack tmp = *this;
  this->m_nRows = size;
  this->m_nCols = 1;
  this->initialize();

  const TlDenseVectorObject::index_type copySize =
      std::min(size, tmp.getSize());
  for (TlDenseVectorObject::index_type i = 0; i < copySize; ++i) {
    this->set(i, tmp.get(i));
  }
}

double TlDenseVector_ImplScalapack::get(
    const TlDenseVectorObject::size_type index) const {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  // DataType& dataTmp = const_cast<DataType&>(this->vector_);

  double dAnswer = 0.0;
  const int nGlobalRow = index + 1;
  const int nGlobalCol = 0 + 1;
  pdelget_("A", " ", &dAnswer, this->vector_, &nGlobalRow, &nGlobalCol,
           this->m_pDESC);

  return dAnswer;
}

void TlDenseVector_ImplScalapack::set(const TlDenseVectorObject::size_type i,
                                      const double dValue) {
  const int index = this->getIndex(i, 0);
  if (index != -1) {
    this->vector_[index] = dValue;
  }
}

void TlDenseVector_ImplScalapack::add(
    const TlDenseVectorObject::size_type index, const double value) {
  const int localIndex = this->getIndex(index, 0);
  if (localIndex != -1) {
    this->vector_[localIndex] += value;
  }
}

void TlDenseVector_ImplScalapack::mul(
    const TlDenseVectorObject::size_type index, const double value) {
  const int localIndex = this->getIndex(index, 0);
  if (localIndex != -1) {
    this->vector_[localIndex] *= value;
  }
}

std::vector<double> TlDenseVector_ImplScalapack::getVector() const {
  std::vector<double> ans(this->getSize(), 0.0);
  const int nRows = this->m_nMyRows;
  const int nCols = this->m_nMyCols;
  const int nGlobalRows = this->m_nRows;
  const int nGlobalCols = this->m_nCols;

  for (int r = 0; r < nRows; ++r) {
    const int nGlobalRowIndex = this->m_RowIndexTable[r];
    if (nGlobalRowIndex >= nGlobalRows) {
      continue;
    }

    for (int c = 0; c < nCols; ++c) {
      const int nGlobalColIndex = this->m_ColIndexTable[c];
      if (nGlobalColIndex >= nGlobalCols) {
        continue;
      }
      const int nGlobalIndex = nGlobalRowIndex * nGlobalCols + nGlobalColIndex;
      const int index = r + nRows * c;  // row-major
      ans[nGlobalIndex] = this->vector_[index];
    }
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(ans);

  return ans;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator=(
    const TlDenseVector_ImplScalapack& rhs) {
  if (&rhs != this) {
    this->m_nRows = rhs.m_nRows;
    this->m_nCols = rhs.m_nCols;
    this->m_nBlockSize = rhs.m_nBlockSize;

    this->initialize();
    std::copy(rhs.vector_, rhs.vector_ + (this->m_nMyRows * this->m_nMyCols),
              this->vector_);
  }

  return *this;
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator+=(
    const TlDenseVector_ImplScalapack& rhs) {
  assert(this->m_nRows == rhs.m_nRows);
  assert(this->m_nCols == rhs.m_nCols);
  assert(this->m_nBlockSize == rhs.m_nBlockSize);

  const TlDenseVectorObject::index_type localSize =
      this->m_nMyRows * this->m_nMyCols;
  std::transform(this->vector_, this->vector_ + localSize, rhs.vector_,
                 this->vector_, std::plus<double>());

  return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator-=(
    const TlDenseVector_ImplScalapack& rhs) {
  assert(this->m_nRows == rhs.m_nRows);
  assert(this->m_nCols == rhs.m_nCols);
  assert(this->m_nBlockSize == rhs.m_nBlockSize);

  const TlDenseVectorObject::index_type localSize =
      this->m_nMyRows * this->m_nMyCols;
  std::transform(this->vector_, this->vector_ + localSize, rhs.vector_,
                 this->vector_, std::minus<double>());

  return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator*=(
    const double coef) {
  const TlDenseVectorObject::index_type localSize =
      this->m_nMyRows * this->m_nMyCols;
  std::transform(this->vector_, this->vector_ + localSize, this->vector_,
                 std::bind1st(std::multiplies<double>(), coef));

  return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator/=(
    const double coef) {
  return this->operator*=(1.0 / coef);
}

double TlDenseVector_ImplScalapack::operator*(
    const TlDenseVector_ImplScalapack& rhs) const {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  TlDenseVector_ImplScalapack& Xtmp =
      const_cast<TlDenseVector_ImplScalapack&>(*this);
  TlDenseVector_ImplScalapack& Ytmp =
      const_cast<TlDenseVector_ImplScalapack&>(rhs);

  assert(this->getSize() == rhs.getSize());

  const int N = this->getSize();
  double DOT = 0.0;
  const int IX = 1;
  const int JX = 1;
  const int INCX = 1;
  const int IY = 1;
  const int JY = 1;
  const int INCY = 1;

  pddot_(&N, &DOT, Xtmp.vector_, &IX, &JX, Xtmp.m_pDESC, &INCX, Ytmp.vector_,
         &IY, &JY, Ytmp.m_pDESC, &INCY);

  // INCX=1の場合、列のプロセスにのみ値が返る。
  // したがって、全プロセスに値を送信する必要がある。
  {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.broadcast(DOT);
  }

  return DOT;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
void TlDenseVector_ImplScalapack::sortByGreater() {
  // TODO: implement
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::dotInPlace(
    const TlDenseVector_ImplScalapack& rhs) {
  assert(this->getSize() == rhs.getSize());

  const TlDenseVectorObject::index_type localSize =
      this->m_nMyRows * this->m_nMyCols;
  std::transform(this->vector_, this->vector_ + localSize, rhs.vector_,
                 this->vector_, std::multiplies<double>());

  return *this;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseVector_ImplScalapack::initialize() {
  TlScalapackContext::getData(this->m_nContext, this->m_nProc, this->m_nRank,
                              this->m_nProcGridRow, this->m_nProcGridCol);

  // my process position on the process matrix
  Cblacs_gridinfo(this->m_nContext, &(this->m_nProcGridRow),
                  &(this->m_nProcGridCol), &(this->m_nMyProcRow),
                  &(this->m_nMyProcCol));

  // determine sizes of local matrix
  const int nStartRowProc = 0;
  const int nStartColProc = 0;
  this->m_nMyRows = std::max(
      1, numroc_(&(this->m_nRows), &(this->m_nBlockSize), &(this->m_nMyProcRow),
                 &nStartRowProc, &(this->m_nProcGridRow)));
  this->m_nMyCols = std::max(
      1, numroc_(&(this->m_nCols), &(this->m_nBlockSize), &(this->m_nMyProcCol),
                 &nStartColProc, &(this->m_nProcGridCol)));

  // make parameter, desca
  int nInfo = 0;
  descinit_(this->m_pDESC, &(this->m_nRows), &(this->m_nCols),
            &(this->m_nBlockSize), &(this->m_nBlockSize), &nStartRowProc,
            &nStartColProc, &(this->m_nContext), &(this->m_nMyRows), &nInfo);
  assert(nInfo == 0);

  // データ用バッファメモリの確保
  if (this->vector_ != NULL) {
    delete[] this->vector_;
    this->vector_ = NULL;
  }
  this->vector_ = new double[this->m_nMyRows * this->m_nMyCols];
  std::fill(this->vector_, this->vector_ + (this->m_nMyRows * this->m_nMyCols),
            0.0);

  // 行方向のglobal_index v.s. local_indexのリストを作成
  {
    const int nMyRows = this->m_nMyRows;
    const int nBlockSize = this->m_nBlockSize;
    const int nBlockIndex =
        this->m_nMyProcRow * nBlockSize;  // 各ローカル行列の最初のインデックス
    const int nIncrementBlockIndex =
        this->m_nProcGridRow * nBlockSize;  // ブロック最初のインデックスの増分
    this->m_RowIndexTable.resize(nMyRows);
    for (int r = 0; r < nMyRows; ++r) {
      const div_t d = std::div(r, nBlockSize);
      this->m_RowIndexTable[r] =
          nBlockIndex + (nIncrementBlockIndex * d.quot) + d.rem;
    }
  }

  // 列方向のglobal_index v.s. local_indexのリストを作成
  {
    const int nMyCols = this->m_nMyCols;
    const int nBlockSize = this->m_nBlockSize;
    const int nBlockIndex =
        this->m_nMyProcCol *
        this->m_nBlockSize;  // 各ローカル行列の最初のインデックス
    const int nIncrementBlockIndex =
        this->m_nProcGridCol *
        this->m_nBlockSize;  // ブロック最初のインデックスの増分
    this->m_ColIndexTable.resize(nMyCols);
    for (int c = 0; c < nMyCols; ++c) {
      const div_t d = std::div(c, nBlockSize);
      this->m_ColIndexTable[c] =
          nBlockIndex + (nIncrementBlockIndex * d.quot) + d.rem;
    }
  }
}

int TlDenseVector_ImplScalapack::getIndex(int nGlobalRow,
                                          int nGlobalCol) const {
  int nAnswer = -1;

  // +1 for fortran code
  ++nGlobalRow;
  ++nGlobalCol;

  int nLocalRow = 0;
  int nLocalCol = 0;
  int nLocalProcRow = 0;
  int nLocalProcCol = 0;
  infog2l_(&nGlobalRow, &nGlobalCol, this->m_pDESC, &(this->m_nProcGridRow),
           &(this->m_nProcGridCol), &(this->m_nMyProcRow),
           &(this->m_nMyProcCol), &nLocalRow, &nLocalCol, &nLocalProcRow,
           &nLocalProcCol);

  if ((this->m_nMyProcRow == nLocalProcRow) &&
      (this->m_nMyProcCol == nLocalProcCol)) {
    nAnswer = (nLocalRow - 1) + (nLocalCol - 1) * (this->m_nMyRows);
  }

  return nAnswer;
}
