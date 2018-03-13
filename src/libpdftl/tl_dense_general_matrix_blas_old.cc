// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

#include "TlLogging.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_matrix_utils.h"
#include "tl_dense_vector_blas.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

const TlDenseGeneralMatrix_BLAS_old::size_type TlDenseGeneralMatrix_BLAS_old::MAX_LOOP =
    std::numeric_limits<int>::max();

TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(const index_type row,
                                                     const index_type col)
    : TlDenseGeneralMatrixAbstract(), data_(NULL) {
  assert((0 < row) && (0 < col));
  this->row_ = row;
  this->col_ = col;
  this->initialize();
}

// for sub-class
TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(const index_type row,
                                                     const index_type col,
                                                     double* pData)
    : TlDenseGeneralMatrixAbstract(), data_(pData) {
  this->row_ = row;
  this->col_ = col;
}

TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(
    const TlDenseGeneralMatrix_BLAS_old& rhs)
    : TlDenseGeneralMatrixAbstract(), data_(NULL) {
  this->row_ = rhs.getNumOfRows();
  this->col_ = rhs.getNumOfCols();
  this->initialize(false);
  const std::size_t size = this->getNumOfRows() * this->getNumOfCols();
  std::copy(rhs.data_, rhs.data_ + size, this->data_);
}

TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(
    const TlDenseSymmetricMatrix_BLAS_Old& rhs)
    : TlDenseGeneralMatrixAbstract(), data_(NULL) {
  this->row_ = rhs.getNumOfRows();
  this->col_ = rhs.getNumOfCols();
  this->initialize(false);

  const index_type numOfRows = this->getNumOfRows();

  // FX1用富士通コンパイラでは
  // OpenMPをコンストラクタで使ってはいけない。
  for (int row = 0; row < numOfRows; ++row) {
    for (int col = 0; col < row; ++col) {
      const double tmp = rhs.get(row, col);

      this->set(row, col, tmp);
      this->set(col, row, tmp);
    }

    // case : row == col
    this->set(row, row, rhs.get(row, row));
  }
}

TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(
    const TlSerializeData& data)
    : TlDenseGeneralMatrixAbstract(), data_(NULL) {
  this->row_ = std::max(data["row"].getInt(), 1);
  this->col_ = std::max(data["col"].getInt(), 1);
  this->initialize(false);

  size_type index = 0;
  const index_type maxRow = this->getNumOfRows();
  const index_type maxCol = this->getNumOfCols();
  for (index_type row = 0; row < maxRow; ++row) {
    for (index_type col = 0; col < maxCol; ++col) {
      this->set(row, col, data["data"].getAt(index).getDouble());
      ++index;
    }
  }
}

TlDenseGeneralMatrix_BLAS_old::TlDenseGeneralMatrix_BLAS_old(const TlVector_BLAS& vct,
                                                     const index_type row,
                                                     const index_type col)
    : TlDenseGeneralMatrixAbstract(), data_(NULL) {
  assert((0 < row) && (0 < col));
  this->row_ = row;
  this->col_ = col;

  const size_type size = this->getNumOfElements();
  assert(vct.getSize() == size);

  this->initialize(false);

#ifndef __FUJITSU  // cannot use OpenMP in constructor
#pragma omp parallel for
#endif  // __FUJITSU
  for (size_type index = 0; index < size; ++index) {
    this->data_[index] = vct[index];
  }
}

TlDenseGeneralMatrix_BLAS_old::~TlDenseGeneralMatrix_BLAS_old() { this->clear(); }

void TlDenseGeneralMatrix_BLAS_old::initialize(bool isZeroClear) {
  const TlMatrixObject::index_type row = this->getNumOfRows();
  const TlMatrixObject::index_type col = this->getNumOfCols();
  const std::size_t size = this->getNumOfElements();
  try {
    this->data_ = new double[size];
  } catch (std::bad_alloc& ba) {
    this->log_.critical(
        TlUtils::format("bad_alloc caught: %s: row=%d, col=%d, size=%ld",
                        ba.what(), row, col, size));
    throw;
  } catch (...) {
    this->log_.critical("unknown error.");
    throw;
  }
  assert(this->data_ != NULL);

  if (isZeroClear == true) {
    std::fill(this->data_, this->data_ + size, 0.0);
  }
}

void TlDenseGeneralMatrix_BLAS_old::clear() {
  if (this->data_ != NULL) {
    delete[] this->data_;
    this->data_ = NULL;
  }
}

void TlDenseGeneralMatrix_BLAS_old::resize(const TlMatrixObject::index_type nRow,
                                       const TlMatrixObject::index_type nCol) {
  assert((nRow > 0) && (nCol > 0));

  TlDenseGeneralMatrix_BLAS_old oldMatrix(*this);

  this->clear();
  this->row_ = nRow;
  this->col_ = nCol;
  this->initialize(true);

  const int nMaxRowsForCopy = std::min(oldMatrix.getNumOfRows(), nRow);
  const int nMaxColsForCopy = std::min(oldMatrix.getNumOfCols(), nCol);
#pragma omp parallel for
  for (int c = 0; c < nMaxColsForCopy; ++c) {
    for (int r = 0; r < nMaxRowsForCopy; ++r) {
      this->set(r, c, oldMatrix.get(r, c));
    }
  }
}

TlMatrixObject::size_type TlDenseGeneralMatrix_BLAS_old::index(
    index_type row, index_type col) const {
  // Column-major order
  // see. https://en.wikipedia.org/wiki/Row-_and_column-major_order
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

  const size_type index = row + this->getNumOfRows() * col;
  return index;
}

double TlDenseGeneralMatrix_BLAS_old::get(const index_type row,
                                      const index_type col) const {
  return this->data_[this->index(row, col)];
}

void TlDenseGeneralMatrix_BLAS_old::set(const index_type row, const index_type col,
                                    const double value) {
  const size_type index = this->index(row, col);

#pragma omp critical(TlDenseGeneralMatrix_BLAS_old__set)
  { this->data_[index] = value; }
}

void TlDenseGeneralMatrix_BLAS_old::add(const index_type row, const index_type col,
                                    const double value) {
  const size_type index = this->index(row, col);

#pragma omp atomic
  this->data_[index] += value;
}

TlVector_BLAS TlDenseGeneralMatrix_BLAS_old::getVector() const {
  return TlVector_BLAS(this->data_, this->getNumOfElements());
}

TlVector_BLAS TlDenseGeneralMatrix_BLAS_old::getRowVector(
    const TlMatrixObject::index_type row) const {
  assert((0 <= row) && (row < this->getNumOfRows()));

  const TlMatrixObject::index_type numOfCols = this->getNumOfCols();
  TlVector_BLAS answer(numOfCols);

  for (TlMatrixObject::index_type i = 0; i < numOfCols; ++i) {
    answer[i] = this->get(row, i);
  }

  return answer;
}

TlVector_BLAS TlDenseGeneralMatrix_BLAS_old::getColVector(
    const TlMatrixObject::index_type col) const {
  assert((0 <= col) && (col < this->getNumOfCols()));

  const TlMatrixObject::index_type numOfRows = this->getNumOfRows();
  TlVector_BLAS answer(numOfRows);

  for (TlMatrixObject::index_type i = 0; i < numOfRows; ++i) {
    answer[i] = this->get(i, col);
  }

  return answer;
}

TlDenseGeneralMatrix_BLAS_old TlDenseGeneralMatrix_BLAS_old::getBlockMatrix(
    const int nRow, const int nCol, const int nRowDistance,
    const int nColDistance) const {
  assert((0 <= nRow) && (nRow < this->getNumOfRows()));
  assert((0 <= nCol) && (nCol < this->getNumOfCols()));
  assert(0 < nRowDistance);
  assert(0 < nColDistance);

  // for debug
  //   {
  //     std::cout << TlUtils::format("matrix size = (%d, %d)",
  //     this->getNumOfRows(), this->getNumOfCols())
  //          << std::endl;
  //     std::cout << TlUtils::format("(%d, %d) -> (%d, %d)", row, col, row +
  //     row_distance, col + col_distance)
  //          << std::endl;
  //   }

  assert(0 <= (nRow + nRowDistance) &&
         (nRow + nRowDistance) <= this->getNumOfRows());
  assert(0 <= (nCol + nColDistance) &&
         (nCol + nColDistance) <= this->getNumOfCols());

  TlDenseGeneralMatrix_BLAS_old answer(nRowDistance, nColDistance);
#pragma omp parallel for
  for (int dr = 0; dr < nRowDistance; ++dr) {
    const int r = nRow + dr;
    for (int dc = 0; dc < nColDistance; ++dc) {
      const int c = nCol + dc;

      answer.set(dr, dc, this->get(r, c));
    }
  }

  return answer;
}

void TlDenseGeneralMatrix_BLAS_old::setBlockMatrix(
    const index_type row, const index_type col,
    const TlDenseGeneralMatrix_BLAS_old& matrix) {
  TlLogging& log = TlLogging::getInstance();
  const int row_distance = matrix.getNumOfRows();
  const int col_distance = matrix.getNumOfCols();

  if (!((0 <= row && row < this->getNumOfRows()) &&
        (0 <= col && col < this->getNumOfCols()) &&
        (0 < (row + row_distance) &&
         (row + row_distance) <= this->getNumOfRows()) &&
        (0 < (col + col_distance) &&
         (col + col_distance) <= this->getNumOfCols()))) {
    log.critical(TlUtils::format(
        "setBlockMatrix() start(%d, %d) mat(%d, %d) -> (%d, %d)", row, col,
        matrix.getNumOfRows(), matrix.getNumOfCols(), this->getNumOfRows(),
        this->getNumOfCols()));
  }
  assert(0 <= row && row < this->getNumOfRows());
  assert(0 <= col && col < this->getNumOfCols());
  assert(0 < (row + row_distance) &&
         (row + row_distance) <= this->getNumOfRows());
  assert(0 < (col + col_distance) &&
         (col + col_distance) <= this->getNumOfCols());

#pragma omp parallel for
  for (int dr = 0; dr < row_distance; ++dr) {
    const int r = row + dr;
    for (int dc = 0; dc < col_distance; ++dc) {
      const int c = col + dc;

      this->set(r, c, matrix.get(dr, dc));
    }
  }
}

void TlDenseGeneralMatrix_BLAS_old::addBlockMatrix(
    const int row, const int col, const TlDenseGeneralMatrix_BLAS_old& matrix) {
  const int row_distance = matrix.getNumOfRows();
  const int col_distance = matrix.getNumOfCols();

  assert(0 <= row && row < this->getNumOfRows());
  assert(0 <= col && col < this->getNumOfCols());
  assert(0 <= (row + row_distance) &&
         (row + row_distance) < this->getNumOfRows());
  assert(0 <= (col + col_distance) &&
         (col + col_distance) < this->getNumOfCols());

#pragma omp parallel for
  for (int dr = 0; dr < row_distance; ++dr) {
    const int r = row + dr;
    for (int dc = 0; dc < col_distance; ++dc) {
      const int c = col + dc;

      this->add(r, c, matrix.get(dr, dc));
    }
  }
}

double TlDenseGeneralMatrix_BLAS_old::trace() const {
  const index_type dim =
      std::min<index_type>(this->getNumOfRows(), this->getNumOfCols());
  double answer = 0.0;
#pragma omp parallel for reduction(+ : answer)
  for (index_type i = 0; i < dim; ++i) {
    answer += this->get(i, i);
  }

  return answer;
}

double TlDenseGeneralMatrix_BLAS_old::sum() const {
  // double answer = this->m_aMatrix.sum();
  return std::accumulate(this->data_, this->data_ + this->getNumOfElements(),
                         0.0);
}

// operator
TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator=(
    const TlDenseGeneralMatrix_BLAS_old& rhs) {
  if (this != &rhs) {
    this->clear();
    this->row_ = rhs.getNumOfRows();
    this->col_ = rhs.getNumOfCols();
    this->initialize(false);
    const std::size_t size = this->getNumOfElements();
    std::copy(rhs.data_, rhs.data_ + size, this->data_);
  }

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator=(
    const TlDenseSymmetricMatrix_BLAS_Old& rhs) {
  this->clear();
  this->row_ = rhs.getNumOfRows();
  this->col_ = rhs.getNumOfCols();
  this->initialize(false);

  const index_type numOfRows = this->getNumOfRows();
  for (index_type row = 0; row < numOfRows; ++row) {
    // case: row != col
    for (index_type col = 0; col < row; ++col) {
      const double tmp = rhs.get(row, col);
      this->set(row, col, tmp);
      this->set(col, row, tmp);
    }

    // case : row == col
    this->set(row, row, rhs.get(row, row));
  }

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator+=(
    const TlDenseGeneralMatrix_BLAS_old& rhs) {
  assert(this->getNumOfRows() == rhs.getNumOfRows());
  assert(this->getNumOfCols() == rhs.getNumOfCols());

  const size_type size = this->getNumOfElements();

#ifdef _OPENMP
#pragma omp critical(TlDenseGeneralMatrix_BLAS_old__operator_addplus)
  {
    // use OpenMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
      for (size_type block = 0; block < quot; ++block) {
        const size_type index_base = block * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < MAX_LOOP; ++i) {
          const size_type index = index_base + i;
          this->data_[index] += rhs.data_[index];
        }
      }

      const size_type index_base = quot * MAX_LOOP;
#pragma omp for
      for (int i = 0; i < rem; ++i) {
        const size_type index = index_base + i;
        this->data_[index] += rhs.data_[index];
      }
    }
  }

#else
  // not use OpenMP
  for (size_type index = 0; index < size; ++index) {
    this->data_[index] += rhs.data_[index];
  }
#endif  // _OPENMP

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator+=(
    const TlDenseSymmetricMatrix_BLAS_Old& rhs) {
  assert(this->getNumOfRows() == rhs.getNumOfRows());
  assert(this->getNumOfCols() == rhs.getNumOfCols());

  const int numOfRows = rhs.getNumOfRows();

#pragma omp parallel for
  for (int row = 0; row < numOfRows; ++row) {
    // row != col
    for (int col = 0; col < row; ++col) {
      const double tmp = rhs.get(row, col);
#pragma omp critical(TlDenseGeneralMatrix_BLAS_old_operator_plus_equal)
      {
        this->add(row, col, tmp);
        this->add(col, row, tmp);
      }
    }

    // row == col
    this->add(row, row, rhs.get(row, row));
  }

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator-=(
    const TlDenseGeneralMatrix_BLAS_old& rhs) {
  assert(this->getNumOfRows() == rhs.getNumOfRows());
  assert(this->getNumOfCols() == rhs.getNumOfCols());

  const size_type size = this->getNumOfElements();

#ifdef _OPENMP
  // use OpenMP
  const size_type quot = size / MAX_LOOP;
  const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
  {
    for (size_type block = 0; block < quot; ++block) {
      const size_type index_base = block * MAX_LOOP;
#pragma omp for
      for (int i = 0; i < MAX_LOOP; ++i) {
        const size_type index = index_base + i;
        this->data_[index] -= rhs.data_[index];
      }
    }

    const size_type index_base = quot * MAX_LOOP;
#pragma omp for
    for (int i = 0; i < rem; ++i) {
      const size_type index = index_base + i;
      this->data_[index] -= rhs.data_[index];
    }
  }
#else
  for (size_type index = 0; index < size; ++index) {
    this->data_[index] -= rhs.data_[index];
  }
#endif  // _OPENMP

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator-=(
    const TlDenseSymmetricMatrix_BLAS_Old& rhs) {
  assert(this->getNumOfRows() == rhs.getNumOfRows());
  assert(this->getNumOfCols() == rhs.getNumOfCols());

  const int numOfRows = rhs.getNumOfRows();
#pragma omp parallel for
  for (int row = 0; row < numOfRows; ++row) {
    // row != col
    for (int col = 0; col < row; ++col) {
      const double tmp = -rhs.get(row, col);
#pragma omp critical(TlDenseGeneralMatrix_BLAS_old_operator_minus_equal)
      {
        this->add(row, col, tmp);
        this->add(col, row, tmp);
      }
    }

    // row == col
    this->add(row, row, rhs.get(row, row));
  }

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator*=(
    const TlDenseGeneralMatrix_BLAS_old& rhs) {
  TlDenseGeneralMatrix_BLAS_old tmp = (*this) * rhs;
  (*this) = tmp;

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator*=(
    const TlDenseSymmetricMatrix_BLAS_Old& rhs) {
  return this->operator*=(TlDenseGeneralMatrix_BLAS_old(rhs));
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator*=(
    const double& rhs) {
  const size_type size = this->getNumOfElements();

#ifdef _OPENMP
  // use OpenMP
  const size_type quot = size / MAX_LOOP;
  const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
  {
    for (size_type block = 0; block < quot; ++block) {
      const size_type index_base = block * MAX_LOOP;
#pragma omp for
      for (int i = 0; i < MAX_LOOP; ++i) {
        const size_type index = index_base + i;
        this->data_[index] *= rhs;
      }
    }

    const size_type index_base = quot * MAX_LOOP;
#pragma omp for
    for (int i = 0; i < rem; ++i) {
      const size_type index = index_base + i;
      this->data_[index] *= rhs;
    }
  }
#else
  // not use OpenMP
  for (size_type index = 0; index < size; ++index) {
    this->data_[index] *= rhs;
  }
#endif  // _OPENMP

  return (*this);
}

TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::operator/=(
    const double& rhs) {
  return (this->operator*=(1.0 / rhs));
}

TlDenseGeneralMatrix_BLAS_old operator+(const TlDenseGeneralMatrix_BLAS_old& X,
                                    const TlDenseGeneralMatrix_BLAS_old& Y) {
  assert(X.getNumOfRows() == Y.getNumOfRows());
  assert(X.getNumOfCols() == Y.getNumOfCols());

  TlDenseGeneralMatrix_BLAS_old answer(X);
  answer += Y;

  return answer;
}

TlDenseGeneralMatrix_BLAS_old operator-(const TlDenseGeneralMatrix_BLAS_old& X,
                                    const TlDenseGeneralMatrix_BLAS_old& Y) {
  assert(X.getNumOfRows() == Y.getNumOfRows());
  assert(X.getNumOfCols() == Y.getNumOfCols());

  TlDenseGeneralMatrix_BLAS_old answer(X);
  answer -= Y;

  return answer;
}

// X x Y
TlDenseGeneralMatrix_BLAS_old operator*(const TlDenseGeneralMatrix_BLAS_old& X,
                                    const TlDenseGeneralMatrix_BLAS_old& Y) {
#ifdef HAVE_LAPACK
  return multiplicationByLapack(X, Y);
#else
  assert(X.getNumOfCols() == Y.getNumOfRows());
  TlDenseGeneralMatrix_BLAS_old answer(X.getNumOfRows(), Y.getNumOfCols());
  for (index_type row = 0; row < X.getNumOfRows(); ++row) {
    for (index_type col = 0; col < Y.getNumOfCols(); ++col) {
      for (index_type t = 0; t < X.getNumOfCols(); ++t) {
        answer(row, col) += X(row, t) * Y(t, col);
      }
    }
  }

  return answer;
#endif  // HAVE_LAPACK
}

TlDenseGeneralMatrix_BLAS_old operator*(const TlDenseGeneralMatrix_BLAS_old& X,
                                    double Y) {
  TlDenseGeneralMatrix_BLAS_old answer(X);
  answer *= Y;

  return answer;
}

// 行列x縦ベクトル
TlVector_BLAS operator*(const TlDenseGeneralMatrix_BLAS_old& A,
                        const TlVector_BLAS& X) {
#ifdef HAVE_LAPACK
  return multiplicationByLapack(A, X);
#else
  {
    assert(A.getNumOfCols() == X.getSize());
    TlVector_BLAS answer(A.getNumOfRows());
    for (int row = 0; row < X.getNumOfRows(); ++row) {
      double tmp = 0.0;
      for (int col = 0; col < X.getNumOfCols(); ++col) {
        tmp += A.get(row, col) * X.get(col);
      }
      answer.set(row, tmp);
    }
    return answer;
  }
#endif  // HAVE_LAPACK
}

// 横ベクトルx行列
TlVector_BLAS operator*(const TlVector_BLAS& X,
                        const TlDenseGeneralMatrix_BLAS_old& A) {
#ifdef HAVE_LAPACK
  return multiplicationByLapack(X, A);
#else
  {
    assert(A.getNumOfRows() == X.getSize());
    TlVector_BLAS answer(A.getNumOfCol());
    for (int col = 0; col < X.getNumOfCols(); ++col) {
      double tmp = 0.0;
      for (int row = 0; row < X.getNumOfRows(); ++row) {
        tmp += A.get(row, col) * X.get(row);
      }
      answer.set(col, tmp);
    }
    return answer;
  }
#endif  // HAVE_LAPACK
}

// ddot?
const TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::dotInPlace(
    const TlDenseGeneralMatrix_BLAS_old& X) {
  assert(this->getNumOfRows() == X.getNumOfRows());
  assert(this->getNumOfCols() == X.getNumOfCols());

  const size_type size = this->getNumOfElements();

#ifdef _OPENMP
  // use OpenMP
  const size_type quot = size / MAX_LOOP;
  const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
  {
    for (size_type block = 0; block < quot; ++block) {
      const size_type index_base = block * MAX_LOOP;
#pragma omp for
      for (int i = 0; i < MAX_LOOP; ++i) {
        const size_type index = index_base + i;
        this->data_[index] *= X.data_[index];
      }
    }

    const size_type index_base = quot * MAX_LOOP;
#pragma omp for
    for (int i = 0; i < rem; ++i) {
      const size_type index = index_base + i;
      this->data_[index] *= X.data_[index];
    }
  }
#else
  // not use OpenMP
  for (size_type index = 0; index < size; ++index) {
    this->data_[index] *= X.data_[index];
  }
#endif  // _OPENMP

  return (*this);
}

const TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::dotInPlace(
    const TlDenseSymmetricMatrix_BLAS_Old& X) {
  this->dotInPlace(TlDenseGeneralMatrix_BLAS_old(X));

  return *this;
}

TlDenseGeneralMatrix_BLAS_old dot(const TlDenseGeneralMatrix_BLAS_old& X,
                              const TlDenseGeneralMatrix_BLAS_old& Y) {
  assert(X.getNumOfRows() == Y.getNumOfRows());
  assert(X.getNumOfCols() == Y.getNumOfCols());

  TlDenseGeneralMatrix_BLAS_old Z = X;
  Z.dotInPlace(Y);

  return Z;
}

std::size_t TlDenseGeneralMatrix_BLAS_old::getMemSize() const {
  std::size_t answer = sizeof(TlDenseGeneralMatrix_BLAS_old);
  answer += sizeof(double) * this->getNumOfElements();

  return answer;
}

bool TlDenseGeneralMatrix_BLAS_old::load(const std::string& filePath) {
  bool answer = false;
  TlMatrixObject::MatrixType matrixType;
  TlMatrixObject::index_type row, col;

  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &matrixType, &row, &col);
  if (headerSize > 0) {
    this->clear();
    this->row_ = row;
    this->col_ = col;
    this->initialize();

    const TlMatrixObject::size_type length = row * col;
    answer = TlDenseGeneralMatrixAbstract::load(filePath, this->data_, length);
  }

  return answer;
}

bool TlDenseGeneralMatrix_BLAS_old::save(const std::string& filePath) const {
  const TlMatrixObject::size_type length =
      this->getNumOfRows() * this->getNumOfCols();
  return TlDenseGeneralMatrixAbstract::save(filePath, this->data_, length);
}

bool TlDenseGeneralMatrix_BLAS_old::saveText(const std::string& sFilePath) const {
  bool bAnswer = true;

  std::ofstream ofs;
  ofs.open(sFilePath.c_str(), std::ofstream::out);

  if (ofs.good()) {
    bAnswer = this->saveText(ofs);
  } else {
    bAnswer = false;
  }

  ofs.close();

  return bAnswer;
}

bool TlDenseGeneralMatrix_BLAS_old::saveText(std::ostream& os) const {
  bool bAnswer = true;
  const int nRows = this->getNumOfRows();
  const int nCols = this->getNumOfCols();

  os << "TEXT\n";
  os << nRows << "\n";
  os << nCols << "\n";

  // print out LCAO coefficent
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      os << TlUtils::format(" %10.6lf", this->get(i, j));
    }
    os << "\n";
  }
  os << "\n";

  return bAnswer;
}

#ifdef HAVE_HDF5
bool TlDenseGeneralMatrix_BLAS_old::saveHdf5(const std::string& filepath,
                                         const std::string& h5path) const {
  TlHdf5Utils h5(filepath);

  assert(this->matrixType_ == CSFD);
  h5.write(h5path, this->data_, this->getNumOfElements());
  h5.setAttr(h5path, "type", static_cast<int>(this->matrixType_));
  h5.setAttr(h5path, "row", this->getNumOfRows());
  h5.setAttr(h5path, "col", this->getNumOfCols());

  return true;
}

bool TlDenseGeneralMatrix_BLAS_old::loadHdf5(const std::string& filepath,
                                         const std::string& h5path) {
  TlHdf5Utils h5(filepath);

  int matrixType;
  h5.getAttr(h5path, "type", &matrixType);

  index_type row = 0;
  index_type col = 0;
  h5.getAttr(h5path, "row", &row);
  h5.getAttr(h5path, "col", &col);
  this->resize(row, col);

  switch (matrixType) {
    case RSFD: {
      std::vector<double> readBuf(this->getNumOfElements());
      h5.get(h5path, &(readBuf[0]), this->getNumOfElements());
      TlMatrixUtils::RSFD2CSFD(this->getNumOfRows(), this->getNumOfCols(),
                               &(readBuf[0]), this->data_);
    } break;

    case CSFD:
      h5.get(h5path, this->data_, this->getNumOfElements());
      break;

    default:
      this->log_.critical(TlUtils::format(
          "illegal matrix type for TlDenseGeneralMatrix_BLAS_old: %d", matrixType));
      break;
  }

  return true;
}
#endif  // HAVE_HDF5

TlSerializeData TlDenseGeneralMatrix_BLAS_old::getSerialize() const {
  TlSerializeData data;
  data["row"] = this->getNumOfRows();
  data["col"] = this->getNumOfCols();
  data["type"] = "RSFD";

  TlSerializeData tmp;
  const size_type size = this->getNumOfElements();
  for (size_type index = 0; index < size; ++index) {
    tmp.pushBack(this->data_[index]);
  }
  data["data"] = tmp;

  return data;
}

/**
 *  出力ストリーム
 */
std::ostream& operator<<(std::ostream& out,
                         const TlDenseGeneralMatrix_BLAS_old& rhs) {
  rhs.print(out);
  return out;
}

/**
 *  CSV出力
 */
std::string TlDenseGeneralMatrix_BLAS_old::getCsv() const {
  std::ostringstream out;

  for (index_type row = 0; row < this->getNumOfRows(); ++row) {
    for (index_type col = 0; col < this->getNumOfCols(); ++col) {
      out << this->get(row, col) << ", ";
    }
    out << std::endl;
  }

  return out.str();
}

// -----------------------------------------------------------------------------
// computation
// -----------------------------------------------------------------------------
// 転置
// TODO: 高速化できる？
const TlDenseGeneralMatrix_BLAS_old& TlDenseGeneralMatrix_BLAS_old::transposeInPlace() {
  const TlMatrixObject::index_type numOfRows = this->getNumOfRows();
  const TlMatrixObject::index_type numOfCols = this->getNumOfCols();

  if (numOfRows == numOfCols) {
    // 正方行列の場合
    for (TlMatrixObject::index_type row = 0; row < numOfRows; ++row) {
      for (TlMatrixObject::index_type col = 0; col < row; ++col) {
        const double tmp = this->get(row, col);
        this->set(row, col, this->get(col, row));
        this->set(col, row, tmp);
      }
    }
  } else {
    // 正方行列でない場合
    TlDenseGeneralMatrix_BLAS_old tmp(*this);
    this->resize(numOfCols, numOfRows);
    for (TlMatrixObject::index_type row = 0; row < numOfRows; ++row) {
      for (TlMatrixObject::index_type col = 0; col < numOfCols; ++col) {
        this->set(col, row, tmp.get(row, col));
      }
    }
  }

  return (*this);
}

bool TlDenseGeneralMatrix_BLAS_old::inverse() {
#ifdef HAVE_LAPACK
  // using LAPACK
  return inverseByLapack(*this);
#else
  // without LAPACK
  std::cerr << "sorry. this code is not implemented." << std::endl;
  abort();
  return false;
#endif  // HAVE_LAPACK
}

TlDenseGeneralMatrix_BLAS_old
TlDenseGeneralMatrix_BLAS_old::solveLinearLeastSquaresProblem(
    const TlDenseGeneralMatrix_BLAS_old& inB) const {
  TlDenseGeneralMatrix_BLAS_old X;
  const bool answer = solveLinearLeastSquaresProblemByLapack(*this, inB, &X);
  if (answer != true) {
    std::cerr << "fail to solve: " << __FILE__ << ", " << __LINE__ << std::endl;
  }

  return X;
}

TlVector_BLAS TlDenseGeneralMatrix_BLAS_old::getDiagonalElements() const {
  const index_type dim = std::min(this->getNumOfRows(), this->getNumOfCols());
  TlVector_BLAS answer(dim);
  for (index_type i = 0; i < dim; ++i) {
    const double value = this->get(i, i);
    answer.set(i, value);
  }

  return answer;
}

// 要素の絶対値の最大値を返す
// int* outRow, outCol にはその要素がある行と列を返す
double TlDenseGeneralMatrix_BLAS_old::getMaxAbsoluteElement(
    TlMatrixObject::index_type* outRow,
    TlMatrixObject::index_type* outCol) const {
  double dAnswer = -1.0;
  int nMaxRow = 0;
  int nMaxCol = 0;

  const int numOfRows = this->getNumOfRows();
  const int numOfCols = this->getNumOfCols();
  for (int row = 0; row < numOfRows; ++row) {
    for (int col = 0; col < numOfCols; ++col) {
      const double val = std::fabs(this->get(row, col));
      if (dAnswer < val) {
        dAnswer = val;
        nMaxRow = row;
        nMaxCol = col;
      }
    }
  }

  if (outRow != NULL) {
    *outRow = nMaxRow;
  }
  if (outCol != NULL) {
    *outCol = nMaxCol;
  }

  return dAnswer;
}

double TlDenseGeneralMatrix_BLAS_old::getRMS() const {
  double sum2 = 0.0;
  const std::size_t numOfDataSize = this->getNumOfElements();
  for (std::size_t i = 0; i < numOfDataSize; ++i) {
    const double tmp = this->data_[i];
    sum2 += tmp * tmp;
  }

  const double elements = this->getNumOfRows() * this->getNumOfCols();

  const double rms = std::sqrt(sum2 / elements);
  return rms;
}

#ifdef HAVE_LAPACK
extern "C" {
void dgemv_(const char* TRANS, const int* M, const int* N, const double* ALPHA,
            const double* A, const int* LDA, const double* X, const int* INCX,
            const double* BETA, double* Y, const int* INCY);

void dgemm_(const char*, const char*, const int*, const int*, const int*,
            const double*, const double*, int*, const double*, const int*,
            const double*, double*, int*);

void dgetrf_(const int* M, const int* N, double* A, const int* LDA, int* IPIV,
             int* INFO);
void dgetri_(const int* N, double* A, const int* LDA, int* IPIV, double* WORK,
             int* LWORK, int* INFO);

void dgelss_(const int* M, const int* N, const int* NRHS, double* A,
             const int* LDA, double* B, const int* LDB, double* S,
             const double* RCOND, int* RANK, double* WORK, const int* LWORK,
             int* INFO);
}

TlVector_BLAS multiplicationByLapack(const TlDenseGeneralMatrix_BLAS_old& A,
                                     const TlVector_BLAS& X) {
  if (A.getNumOfCols() != X.getSize()) {
    std::cerr
        << TlUtils::format(
               "multiplicationByLapack(): parameter error at file:%s, line:%d",
               __FILE__, __LINE__)
        << std::endl;
    std::abort();
  }

  TlVector_BLAS Y(A.getNumOfRows());
  const char TRANS = 'N';
  const int M = A.getNumOfRows();
  const int N = A.getNumOfCols();
  const double alpha = 1.0;
  const int LDA = std::max(1, M);
  const int INCX = 1;
  const double beta = 1.0;
  const int INCY = 1;

  // *  DGEMV  performs one of the matrix-vector operations
  // *
  // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  // *
  // *  where alpha and beta are scalars, x and y are vectors and A is an
  // *  m by n matrix.
  dgemv_(&TRANS, &M, &N, &alpha,
         const_cast<TlDenseGeneralMatrix_BLAS_old&>(A).data_, &LDA,
         const_cast<TlVector_BLAS&>(X).data_, &INCX, &beta, Y.data_, &INCY);

  return Y;
}

TlVector_BLAS multiplicationByLapack(const TlVector_BLAS& X,
                                     const TlDenseGeneralMatrix_BLAS_old& A) {
  if (X.getSize() != A.getNumOfRows()) {
    std::cerr
        << TlUtils::format(
               "multiplicationByLapack(): parameter error at file:%s, line:%d",
               __FILE__, __LINE__)
        << std::endl;
    std::abort();
  }

  TlVector_BLAS Y(A.getNumOfCols());
  const char TRANS = 'T';
  const int M = A.getNumOfRows();
  const int N = A.getNumOfCols();
  const double alpha = 1.0;
  const int LDA = std::max(1, M);
  const int INCX = 1;
  const double beta = 1.0;
  const int INCY = 1;

  // *  DGEMV  performs one of the matrix-vector operations
  // *
  // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  // *
  // *  where alpha and beta are scalars, x and y are vectors and A is an
  // *  m by n matrix.
  dgemv_(&TRANS, &M, &N, &alpha,
         const_cast<TlDenseGeneralMatrix_BLAS_old&>(A).data_, &LDA,
         const_cast<TlVector_BLAS&>(X).data_, &INCX, &beta, Y.data_, &INCY);

  return Y;
}

TlDenseGeneralMatrix_BLAS_old multiplicationByLapack(
    const TlDenseGeneralMatrix_BLAS_old& X, const TlDenseGeneralMatrix_BLAS_old& Y) {
  // std::cerr << "use dgemm() of LAPACK" << std::endl;
  if (X.getNumOfCols() != Y.getNumOfRows()) {
    std::cerr
        << TlUtils::format(
               "multiplicationByLapack(): parameter error at file:%s, line:%d",
               __FILE__, __LINE__)
        << std::endl;
    std::abort();
  }

  TlDenseGeneralMatrix_BLAS_old answer(X.getNumOfRows(), Y.getNumOfCols());

  char TRANSA = 'N';
  char TRANSB = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int M = X.getNumOfRows();
  int N = Y.getNumOfCols();
  int K = X.getNumOfCols();
  int LDA = M;
  int LDB = K;
  int LDC = M;

  // DGEMM  performs one of the matrix-matrix operations
  // C := alpha*op( A )*op( B ) + beta*C,
  // where  op( X ) is one of op( X ) = X   or   op( X ) = X',
  // alpha and beta are scalars, and A, B and C are matrices, with op( A )
  // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha,
         const_cast<TlDenseGeneralMatrix_BLAS_old&>(X).data_, &LDA,
         const_cast<TlDenseGeneralMatrix_BLAS_old&>(Y).data_, &LDB, &beta,
         answer.data_, &LDC);

  return answer;
}

bool inverseByLapack(TlDenseGeneralMatrix_BLAS_old& X) {
  bool bAnswer = false;

  const int M = X.getNumOfRows();
  const int N = X.getNumOfCols();

  const int LDA = std::max(1, M);
  int* IPIV = new int[std::min(M, N)];
  int INFO = 0;

  double* A = X.data_;
  int LWORK = std::max(1, N);
  double* WORK = new double[LWORK];

  dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
  if (INFO == 0) {
    dgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
    if (INFO == 0) {
      bAnswer = true;
    } else {
      std::cerr << "inverseByLapack() failed.: dgetri() return code = " << INFO
                << std::endl;
    }
  } else {
    std::cerr << "inverseByLapack() failed.: dgetrf() return code = " << INFO
              << std::endl;
  }

  std::swap(X.row_, X.col_);

  delete[] WORK;
  WORK = NULL;

  delete[] IPIV;
  IPIV = NULL;

  return bAnswer;
}

bool solveLinearLeastSquaresProblemByLapack(
    const TlDenseGeneralMatrix_BLAS_old& inA, const TlDenseGeneralMatrix_BLAS_old& inB,
    TlDenseGeneralMatrix_BLAS_old* pX) {
  bool answer = true;
  TlLogging& log = TlLogging::getInstance();

  const int M = inA.getNumOfRows();
  if (M != inB.getNumOfRows()) {
    log.critical(
        TlUtils::format("the numbers of A and B are not consistent: %d != %d",
                        M, inB.getNumOfRows()));
    answer = false;
    return answer;
  }
  const int N = inA.getNumOfCols();
  const int NRHS = inB.getNumOfCols();
  TlDenseGeneralMatrix_BLAS_old A = inA;
  const int LDA = std::max(1, M);
  const int LDB = std::max(1, std::max(M, N));
  TlDenseGeneralMatrix_BLAS_old B = inB;
  B.resize(LDB, NRHS);  // size extended.
  double* pS = new double[std::min(M, N)];
  const double RCOND =
      -1.0;  // If RCOND < 0, machine precision is used instead.
  int RANK = 0;
  const int LWORK =
      3 * std::min(M, N) +
      std::max(std::max(2 * std::min(M, N), std::max(M, N)), NRHS);
  double* pWORK = new double[std::max(1, LWORK)];
  int INFO = 0;

  dgelss_(&M, &N, &NRHS, A.data_, &LDA, B.data_, &LDB, pS, &RCOND, &RANK, pWORK,
          &LWORK, &INFO);
  if (INFO != 0) {
    answer = false;
    if (INFO > 0) {
      log.critical(
          TlUtils::format("%d-th argument had an illegal value.", -INFO));
    } else {
      log.critical(
          TlUtils::format("%d-th off-diagonal elements of an intermediate "
                          "bidiagonal form did not converge to zero.",
                          INFO));
    }
  }

  delete[] pS;
  pS = NULL;

  *pX = B;
  pX->resize(N, NRHS);

  return answer;
}

#endif  // HAVE_LAPACK
