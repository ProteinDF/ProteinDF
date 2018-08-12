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

#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include "TlLogging.h"
#include "TlMath.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "scalapack.h"
#include "tl_dense_symmetric_matrix_blacs.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_io.h"
#include "tl_dense_vector_blas.h"
#include "tl_matrix_utils.h"

TlDenseSymmetricMatrix_blacs::TlDenseSymmetricMatrix_blacs(const index_type dim)
    : TlDenseGeneralMatrix_blacs(dim, dim) {}

TlDenseSymmetricMatrix_blacs::TlDenseSymmetricMatrix_blacs(
    const TlDenseSymmetricMatrix_blacs& rhs)
    : TlDenseGeneralMatrix_blacs((TlDenseGeneralMatrix_blacs&)rhs) {}

TlDenseSymmetricMatrix_blacs::TlDenseSymmetricMatrix_blacs(
    const TlDenseGeneralMatrix_blacs& rhs)
    : TlDenseGeneralMatrix_blacs(rhs) {
  assert(rhs.getNumOfRows() == rhs.getNumOfCols());
  // コピーされたバッファの下半分しか使わない
}

TlDenseSymmetricMatrix_blacs::TlDenseSymmetricMatrix_blacs(
    const TlDistributedVector& rhs, const index_type dim)
    : TlDenseGeneralMatrix_blacs(rhs, dim, dim) {}

TlDenseSymmetricMatrix_blacs::~TlDenseSymmetricMatrix_blacs() {}

TlMatrixObject::size_type TlDenseSymmetricMatrix_blacs::getNumOfElements()
    const {
  return TlMatrixObject::getNumOfElements_RLHD();
}

void TlDenseSymmetricMatrix_blacs::resize(const index_type size) {
  assert(size > 0);
  if (size == this->getNumOfRows()) {
    // no need this operation
    return;
  }

  // backup old object
  TlDenseSymmetricMatrix_blacs tmp(*this);

  // new size and zero clear
  this->m_nRows = size;
  this->m_nCols = size;
  this->initialize();

  // copy
  const index_type copyMaxRows = std::min(size, tmp.getNumOfRows());
  const index_type copyMaxCols = std::min(size, tmp.getNumOfCols());
  assert(this->getBlockSize() == tmp.getBlockSize());
  {
    // ブロックサイズが同じでプロセスグリッドが同一なら
    // ローカルコピーで十分。
    std::vector<index_type>::const_iterator rBegin =
        tmp.m_RowIndexTable.begin();
    std::vector<index_type>::const_iterator cBegin =
        tmp.m_ColIndexTable.begin();
    std::vector<index_type>::const_iterator rEnd = std::lower_bound(
        tmp.m_RowIndexTable.begin(), tmp.m_RowIndexTable.end(), copyMaxRows);
    std::vector<index_type>::const_iterator cEnd = std::lower_bound(
        tmp.m_ColIndexTable.begin(), tmp.m_ColIndexTable.end(), copyMaxCols);
    for (std::vector<index_type>::const_iterator r = rBegin; r != rEnd; ++r) {
      const index_type row = *r;
      for (std::vector<index_type>::const_iterator c = cBegin; c != cEnd; ++c) {
        const index_type col = *c;
        if (row < col) {
          break;
        }

        const index_type oldIndex =
            std::distance(rBegin, r) + std::distance(cBegin, c) * tmp.m_nMyRows;
        assert(0 <= oldIndex);
        assert(oldIndex < tmp.getNumOfMyElements());

        const double value = tmp.pData_[oldIndex];
        this->set(row, col, value);
      }
    }
  }
}

TlDenseSymmetricMatrix_blacs& TlDenseSymmetricMatrix_blacs::operator=(
    const TlDenseSymmetricMatrix_blacs& rhs) {
  if (&rhs != this) {
    TlDenseGeneralMatrix_blacs::operator=(rhs);
  }

  return (*this);
}

TlDenseSymmetricMatrix_blacs& TlDenseSymmetricMatrix_blacs::operator=(
    const TlDenseGeneralMatrix_blacs& rhs) {
  if (&rhs != this) {
    TlDenseSymmetricMatrix_blacs tmp(rhs);
    *this = tmp;
  }

  return (*this);
}

int TlDenseSymmetricMatrix_blacs::getProcIdForIndex(
    index_type globalRow, index_type globalCol) const {
  if (globalRow < globalCol) {
    std::swap(globalRow, globalCol);
  }

  return TlDenseGeneralMatrix_blacs::getProcIdForIndex(globalRow, globalCol);
}

TlMatrixObject::size_type TlDenseSymmetricMatrix_blacs::getIndex(
    index_type globalRow, index_type globalCol) const {
  if (globalRow < globalCol) {
    std::swap(globalRow, globalCol);
  }

  return TlDenseGeneralMatrix_blacs::getIndex(globalRow, globalCol);
}

double TlDenseSymmetricMatrix_blacs::get(index_type row, index_type col) const {
  if (row < col) {
    std::swap(row, col);
  }
  return TlDenseGeneralMatrix_blacs::get(row, col);
}

void TlDenseSymmetricMatrix_blacs::set(index_type row, index_type col,
                                       const double value) {
  if (row < col) {
    std::swap(row, col);
  }
  TlDenseGeneralMatrix_blacs::set(row, col, value);
}

void TlDenseSymmetricMatrix_blacs::add(index_type row, index_type col,
                                       const double value) {
  if (row < col) {
    std::swap(row, col);
  }
  TlDenseGeneralMatrix_blacs::add(row, col, value);
}

double TlDenseSymmetricMatrix_blacs::getLocal(index_type row,
                                              index_type col) const {
  if (row < col) {
    std::swap(row, col);
  }
  return TlDenseGeneralMatrix_blacs::getLocal(row, col);
}

TlDenseVector_Lapack TlDenseSymmetricMatrix_blacs::getRowVector(
    const index_type inRow) const {
  assert((0 <= inRow) && (inRow < this->getNumOfRows()));
  const index_type numOfGlobalRows = this->getNumOfRows();
  const index_type numOfGlobalCols = this->getNumOfCols();
  std::vector<double> v(numOfGlobalCols, 0.0);

  if (std::binary_search(this->m_RowIndexTable.begin(),
                         this->m_RowIndexTable.end(), inRow) == true) {
    const index_type row = inRow;
    std::vector<index_type>::const_iterator pEnd = std::upper_bound(
        this->m_ColIndexTable.begin(), this->m_ColIndexTable.end(), row);
    for (std::vector<index_type>::const_iterator p =
             this->m_ColIndexTable.begin();
         p != pEnd; ++p) {
      const index_type col = *p;
      assert(row >= col);
      if (col < numOfGlobalCols) {
        const index_type local_index = this->getIndex(row, col);
        assert(local_index != -1);
        v[col] = this->pData_[local_index];
      }
    }
  }

  // row-col 交換
  if (std::binary_search(this->m_ColIndexTable.begin(),
                         this->m_ColIndexTable.end(), inRow) == true) {
    const index_type col = inRow;
    std::vector<index_type>::const_iterator pBegin = std::lower_bound(
        this->m_RowIndexTable.begin(), this->m_RowIndexTable.end(), col);
    std::vector<index_type>::const_iterator pEnd = this->m_RowIndexTable.end();
    for (std::vector<index_type>::const_iterator p = pBegin; p != pEnd; ++p) {
      const index_type row = *p;
      assert(row >= col);
      if (row < numOfGlobalRows) {
        const index_type local_index = this->getIndex(row, col);
        assert(local_index != -1);
        v[row] = this->pData_[local_index];
      }
    }
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(&(v[0]), numOfGlobalCols);

  return TlDenseVector_Lapack(v);
}

TlDenseVector_Lapack TlDenseSymmetricMatrix_blacs::getColumnVector(
    const index_type col) const {
  return this->getRowVector(col);
}

TlDenseVector_Lapack TlDenseSymmetricMatrix_blacs::getPartialMatrix(
    int* pStartRow, int* pEndRow, int* pStartCol, int* pEndCol) const {
  assert(pStartRow != NULL);
  assert(pEndRow != NULL);
  assert(pStartCol != NULL);
  assert(pEndCol != NULL);

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int proc = rComm.getNumOfProc();
  const int rank = rComm.getRank();

  const int globalSize = this->getNumOfRows();
  assert(this->getNumOfRows() == this->getNumOfCols());
  const int range = globalSize * (globalSize + 1) / 2;

  const int interval = (range + (proc - 1)) / proc;
  const int localStart = interval * rank;
  const int localEnd = localStart + interval - 1;

  TlDenseVector_Lapack P(interval);
  for (int i = 0; i < proc; ++i) {
    index_type numOfLocalRows = 0;
    index_type numOfLocalCols = 0;
    std::vector<int> rowIndexTable;
    std::vector<int> colIndexTable;
    double* pBuf = NULL;

    if (i == rank) {
      numOfLocalRows = this->m_nMyRows;
      numOfLocalCols = this->m_nMyCols;
      rowIndexTable = this->m_RowIndexTable;
      colIndexTable = this->m_ColIndexTable;
      const std::size_t bufSize = numOfLocalRows * numOfLocalCols;
      pBuf = new double[bufSize];
      std::copy(this->pData_, this->pData_ + bufSize, pBuf);
    }
    rComm.broadcast(numOfLocalRows, i);
    rComm.broadcast(numOfLocalCols, i);
    rComm.broadcast(rowIndexTable, i);
    rComm.broadcast(colIndexTable, i);
    const std::size_t bufSize = numOfLocalRows * numOfLocalCols;
    if (i != rank) {
      pBuf = new double[bufSize];
    }
    rComm.broadcast(pBuf, bufSize, i);

    for (int c = 0; c < numOfLocalCols; ++c) {
      const int globalColIndex = colIndexTable[c];

      for (int r = 0; r < numOfLocalRows; ++r) {
        const int globalRowIndex = rowIndexTable[r];

        const int globalIndex =
            globalRowIndex +
            (2 * globalSize - globalColIndex) * (globalColIndex - 1) / 2;
        if (localStart == globalIndex) {
          *pStartRow = globalRowIndex;
          *pStartCol = globalColIndex;
        }
        if (localEnd == globalIndex) {
          *pEndRow = globalRowIndex;
          *pEndCol = globalColIndex;
        }

        if ((localStart <= globalIndex) && (globalIndex <= localEnd)) {
          const std::size_t index = r + numOfLocalRows * c;  // row-major
          const double value = pBuf[index];
          P.set((globalIndex - localStart), value);
        }
      }
    }

    if (pBuf != NULL) {
      delete[] pBuf;
      pBuf = NULL;
    }
  }

  return P;
}

// 複数スレッドから同時に呼び出さない
bool TlDenseSymmetricMatrix_blacs::getSparseMatrixX(
    TlSparseSymmetricMatrix* pMatrix, bool isFinalize) const {
  return TlDenseGeneralMatrix_blacs::getSparseMatrixX(pMatrix, isFinalize);
}

void TlDenseSymmetricMatrix_blacs::mergeSparseMatrix(
    const TlSparseSymmetricMatrix& M) {
  assert(M.getNumOfRows() == this->getNumOfRows());
  assert(M.getNumOfCols() == this->getNumOfCols());

  // 送信すべきインデックスリストの作成
  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int numOfProcs = rComm.getNumOfProc();
  std::vector<std::vector<index_type> > indexArrays(numOfProcs);
  std::vector<std::vector<double> > values(numOfProcs);

  TlSparseMatrix::const_iterator itEnd = M.end();
  for (TlSparseMatrix::const_iterator it = M.begin(); it != itEnd; ++it) {
    index_type globalRow = it->first.row;
    index_type globalCol = it->first.col;
    if (globalRow < globalCol) {
      std::swap(globalRow, globalCol);
    }
    const double value = it->second;

    const int index = this->getIndex(globalRow, globalCol);
    if (index != -1) {
      this->pData_[index] += value;
    } else {
      const int targetProc = this->getProcIdForIndex(globalRow, globalCol);
      indexArrays[targetProc].push_back(globalRow);
      indexArrays[targetProc].push_back(globalCol);
      values[targetProc].push_back(value);
    }
  }

  this->mergeMatrix_common(indexArrays, values);
}

void TlDenseSymmetricMatrix_blacs::mergeSparseMatrixAsync(
    const TlSparseMatrix* pMatrix, bool isFinalize) {
  if (pMatrix != NULL) {
    assert(pMatrix->getNumOfRows() == this->getNumOfRows());
    assert(pMatrix->getNumOfCols() == this->getNumOfCols());

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();

    // 送信すべきインデックスリストの作成
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    TlSparseSymmetricMatrix::const_iterator itEnd = pMatrix->end();
    for (TlSparseSymmetricMatrix::const_iterator it = pMatrix->begin();
         it != itEnd; ++it) {
      index_type globalRow = it->first.row;
      index_type globalCol = it->first.col;
      const double value = it->second;

      if (globalRow < globalCol) {
        std::swap(globalRow, globalCol);
      }

      const int index = this->getIndex(globalRow, globalCol);
      if (index != -1) {
        this->pData_[index] += value;
      } else {
        const int targetProc = this->getProcIdForIndex(globalRow, globalCol);
        indexArrays[targetProc].push_back(globalRow);
        indexArrays[targetProc].push_back(globalCol);
        values[targetProc].push_back(value);
      }
    }

    this->mergeMatrixAsync_send(indexArrays, values);
  }

  this->mergeMatrixAsync_recv(isFinalize);
}

double TlDenseSymmetricMatrix_blacs::getMaxAbsoluteElement(int* pOutRow,
                                                           int* pOutCol) const {
  return TlDenseGeneralMatrix_blacs::getMaxAbsoluteElement(pOutRow, pOutCol);
}

double TlDenseSymmetricMatrix_blacs::getRMS() const {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  double sum2 = 0.0;

  // 対角項と下半分非対角項の二乗和
  const int numOfLocalRows = this->m_nMyRows;
  const int numOfLocalCols = this->m_nMyCols;
  for (int localRowIndex = 0; localRowIndex < numOfLocalRows; ++localRowIndex) {
    const index_type row = this->m_RowIndexTable[localRowIndex];
    for (int localColIndex = 0; localColIndex < numOfLocalCols;
         ++localColIndex) {
      const index_type col = this->m_ColIndexTable[localColIndex];
      if (row > col) {
        const double tmp = this->getLocal(row, col);
        sum2 += 2.0 * tmp * tmp;
      } else if (row == col) {
        const double tmp = this->getLocal(row, col);
        sum2 += tmp * tmp;
      }
    }
  }

  rComm.allReduce_SUM(sum2);

  const double elements = this->getNumOfRows() * this->getNumOfCols();

  const double rms = std::sqrt(sum2 / elements);
  return rms;
}

double TlDenseSymmetricMatrix_blacs::sum() const {
  double answer = 0.0;

  // 対角項と下半分非対角項の和
  const int numOfLocalRows = this->m_nMyRows;
  const int numOfLocalCols = this->m_nMyCols;
  for (int localRowIndex = 0; localRowIndex < numOfLocalRows; ++localRowIndex) {
    const index_type row = this->m_RowIndexTable[localRowIndex];
    for (int localColIndex = 0; localColIndex < numOfLocalCols;
         ++localColIndex) {
      const index_type col = this->m_ColIndexTable[localColIndex];
      if (row > col) {
        answer += this->getLocal(row, col) * 2.0;
      } else if (row == col) {
        answer += this->getLocal(row, col);
      }
    }
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(answer);

  return answer;
}

// =============================================================================
// matrix I/O
// =============================================================================
// bool TlDenseSymmetricMatrix_blacs::isLoadable(std::ifstream& ifs) {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   bool answer = false;
//
//   if (rComm.isMaster() == true) {
//     answer = TlDenseSymmetricMatrix_BLAS_Old::isLoadable(ifs);
//   }
//   rComm.broadcast(answer);
//   return answer;
// }
//
// bool TlDenseSymmetricMatrix_blacs::isLoadable(const std::string& rFilePath) {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   bool answer = false;
//
//   if (rComm.isMaster() == true) {
//     answer = TlDenseSymmetricMatrix_BLAS_Old::isLoadable(rFilePath);
//   }
//   rComm.broadcast(answer);
//   return answer;
// }

bool TlDenseSymmetricMatrix_blacs::load(const std::string& sFilePath) {
  if (TlDenseGeneralMatrix_blacs::isUsingPartialIO == true) {
    return this->loadLocal(sFilePath);
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();

  std::fstream fs;
  if (rComm.isMaster() == true) {
    fs.open(sFilePath.c_str(), std::ios::binary | std::ios::in);
  }

  bool bIsFail = false;
  if (rComm.isMaster() == true) {
    bIsFail = fs.fail();
  }
  rComm.broadcast(bIsFail);
  if (bIsFail) {
#ifdef DEBUG
    std::cerr
        << "[error] TlDenseSymmetricMatrix_blacs::load(): could not open file. "
        << sFilePath << std::endl;
#endif  // DEBUG
    return false;
  }

  bool bAnswer = this->load(fs);

  if (bAnswer != true) {
    if (rComm.isMaster() == true) {
      std::cerr << "TlDenseSymmetricMatrix_blacs::load() is not supported: "
                << sFilePath << std::endl;
      std::abort();
    }
  }

  if (rComm.isMaster() == true) {
    fs.close();
  }

  rComm.broadcast(bAnswer);
  return bAnswer;
}

bool TlDenseSymmetricMatrix_blacs::load(std::fstream& fs) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());

  const int numOfProcs = rComm.getNumOfProc();

  // read header
  bool answer = true;
  TlMatrixObject::MatrixType matType;
  index_type dim;
  if (rComm.isMaster() == true) {
    TlMatrixUtils::getHeaderInfo(fs, &matType, &dim);
    assert(matType == TlMatrixObject::RLHD);
  }
  // rComm.broadcast(matType);
  rComm.broadcast(dim);

  // switch (matType) {
  //   case RLHD:
  //     break;
  //
  //   case CLHD:
  //     break;
  //
  //   default:
  //     if (rComm.isMaster() == true) {
  //       this->log_.critical("this matrix type is not supported. stop.");
  //     }
  //     answer = false;
  //     break;
  // }

  this->m_nRows = dim;
  this->m_nCols = dim;
  assert(this->m_nRows == this->m_nCols);
  this->initialize();

  if (rComm.isMaster() == true) {
    static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
    std::vector<double> buf(bufferCount, 0.0);

    index_type row = 0;
    index_type col = 0;
    size_type count = 0;
    const size_type maxCount = this->getNumOfElements();
    bool isFinished = false;

    std::vector<index_type> sizeLists(numOfProcs, 0);
    std::vector<std::vector<TlMatrixObject::MatrixElement> > elementsBuf(
        numOfProcs);
    std::vector<bool> isSendData(numOfProcs, false);

    while (isFinished == false) {
      // buffer分を一度に読み込み
      fs.read((char*)(&buf[0]), sizeof(double) * bufferCount);

      // 各プロセスのバッファに振り分ける
      std::vector<std::vector<TlMatrixObject::MatrixElement> > elements(
          numOfProcs);
      for (std::size_t i = 0; i < bufferCount; ++i) {
        const int proc = this->getProcIdForIndex(row, col);
        elements[proc].push_back(
            TlMatrixObject::MatrixElement(row, col, buf[i]));

        {
          ++col;
          if (col > row) {
            col = 0;
            ++row;
          }
        }

        ++count;
        if (count >= maxCount) {
          isFinished = true;
          break;
        }
      }

      // データを送信
      for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0は送信しない
        if (isSendData[proc] == true) {
          rComm.wait(sizeLists[proc]);
          if (sizeLists[proc] > 0) {
            rComm.wait(&(elementsBuf[proc][0]));
          }
          isSendData[proc] = false;
        }

        sizeLists[proc] = elements[proc].size();
        elementsBuf[proc] = elements[proc];

        rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
        if (sizeLists[proc] > 0) {
          rComm.iSendDataX(&(elementsBuf[proc][0]), sizeLists[proc], proc,
                           TAG_LOAD_VALUES);
        }
        isSendData[proc] = true;
      }

      // proc=0分データの書き込み
      {
        const index_type sizeList = elements[0].size();
        for (index_type i = 0; i < sizeList; ++i) {
          const index_type r = elements[0][i].row;
          const index_type c = elements[0][i].col;
          const size_type index = this->getIndex(r, c);
          assert(index != -1);
          this->pData_[index] = elements[0][i].value;
        }
      }
    }  // end while

    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      if (isSendData[proc] == true) {
        rComm.wait(sizeLists[proc]);
        if (sizeLists[proc] > 0) {
          rComm.wait(&(elementsBuf[proc][0]));
        }
        isSendData[proc] = false;
      }
    }

    assert(rComm.checkNonBlockingCommunications());

    // 終了メッセージを全ノードに送る
    std::vector<int> endMsg(numOfProcs, 0);
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.iSendData(endMsg[proc], proc, TAG_LOAD_END);
    }
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.wait(endMsg[proc]);
    }
  } else {
    // slave
    const int root = 0;
    index_type sizeList = 0;
    std::vector<TlMatrixObject::MatrixElement> elements;
    int endMsg = 0;

    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
    rComm.iReceiveData(endMsg, root, TAG_LOAD_END);
    bool isBreakLoop = false;
    while (isBreakLoop == false) {
      if (rComm.test(sizeList) == true) {
        rComm.wait(sizeList);
        if (sizeList > 0) {
          elements.resize(sizeList);
          rComm.receiveDataX(&(elements[0]), sizeList, root, TAG_LOAD_VALUES);

          for (index_type i = 0; i < sizeList; ++i) {
            const index_type row = elements[i].row;
            const index_type col = elements[i].col;
            const size_type index = this->getIndex(row, col);
            assert(index != -1);
            this->pData_[index] = elements[i].value;
          }
        }
        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
      }

      if (rComm.test(endMsg) == true) {
        rComm.wait(endMsg);
        rComm.cancel(sizeList);
        isBreakLoop = true;
      }
    }
  }

  assert(rComm.checkNonBlockingCommunications());
  return answer;
}

bool TlDenseSymmetricMatrix_blacs::save(const std::string& sFilePath) const {
  if (TlDenseGeneralMatrix_blacs::isUsingPartialIO == true) {
    return this->saveLocal(sFilePath);
  }

  bool answer = true;

  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());

  if (rComm.isMaster() == true) {
    // master
    assert(this->m_nRows == this->m_nCols);
    const int globalSize = this->m_nRows;
    TlFileSymmetricMatrix fm(sFilePath, globalSize);

    // store local matrix
    {
      const std::vector<TlMatrixObject::MatrixElement> buf =
          this->getMatrixElementsInLocal();
      TlDenseGeneralMatrix_blacs::saveElements(&fm, buf);
    }

    // recive submatrix & write
    const int numOfSlaves = rComm.getNumOfProc() - 1;
    for (int i = 0; i < numOfSlaves; ++i) {
      int src = 0;
      int tag = TAG_SAVE_HANDSHAKE;
      size_type bufSize = 0;
      rComm.iReceiveDataFromAnySource(bufSize, tag);
      rComm.wait(&bufSize, &src);

      std::vector<TlMatrixObject::MatrixElement> buf(bufSize);
      rComm.receiveDataX(&(buf[0]), bufSize, src, TAG_SAVE_DATA);

      TlDenseGeneralMatrix_blacs::saveElements(&fm, buf);
    }
  } else {
    // slave: send submatrix
    const int root = 0;
    const std::vector<TlMatrixObject::MatrixElement> buf =
        this->getMatrixElementsInLocal();
    const size_type bufSize = buf.size();
    rComm.sendData(bufSize, root, TAG_SAVE_HANDSHAKE);
    rComm.iSendDataX(&(buf[0]), buf.size(), root, TAG_SAVE_DATA);
    rComm.wait(&(buf[0]));
  }

  rComm.broadcast(answer);
  assert(rComm.checkNonBlockingCommunications());
  return answer;
}

// HDF5 ------------------------------------------------------------------------
#ifdef HAVE_HDF5
bool TlDenseSymmetricMatrix_blacs::saveHdf5(const std::string& filepath,
                                            const std::string& h5path) const {
  return TlDenseGeneralMatrix_blacs::saveHdf5(filepath, h5path,
                                              TlMatrixObject::RLHD);
}

TlMatrixObject::size_type TlDenseSymmetricMatrix_blacs::getArrayIndex(
    const index_type row, const index_type col) const {
  assert(row >= col);
  const size_type index = row * (row + 1) / 2 + col;
  assert(index < this->getNumOfElements());

  return index;
}

bool TlDenseSymmetricMatrix_blacs::loadHdf5(const std::string& filepath,
                                            const std::string& h5path) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());

  const int numOfProcs = rComm.getNumOfProc();

  // read header
  bool answer = true;
  int matType = 0;
  std::vector<index_type> rowcol(2);
  if (rComm.isMaster() == true) {
    TlHdf5Utils h5(filepath);
    h5.getAttr(h5path, "type", &matType);
    h5.getAttr(h5path, "row", &(rowcol[0]));
    h5.getAttr(h5path, "col", &(rowcol[1]));
  }
  rComm.broadcast(matType);
  rComm.broadcast(&(rowcol[0]), 2, 0);

  switch (matType) {
    case RLHD:
      break;

    default:
      this->log_.critical(TlUtils::format(
          "this matrix type is not supported. stop. type: %d", matType));
      answer = false;
      break;
  }

  this->m_nRows = rowcol[0];
  this->m_nCols = rowcol[1];
  if (this->m_nRows != this->m_nCols) {
    this->log_.critical(
        TlUtils::format("The dims of symmetric matrix are mismatch: (%d, %d)",
                        this->m_nRows, this->m_nCols));
  }
  this->initialize();

  if (rComm.isMaster() == true) {
    TlHdf5Utils h5(filepath);

    static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
    std::vector<double> buf(bufferCount, 0.0);

    index_type row = 0;
    index_type col = 0;
    size_type count = 0;
    const size_type maxCount = this->getNumOfElements();
    bool isFinished = false;

    std::vector<index_type> sizeLists(numOfProcs, 0);
    std::vector<std::vector<TlMatrixObject::MatrixElement> > elementsBuf(
        numOfProcs);
    std::vector<bool> isSendData(numOfProcs, false);

    int cycle = 0;
    while (isFinished == false) {
      // buffer分を一度に読み込み
      const size_type numOfLoadItems =
          std::min<size_type>((cycle + 1) * bufferCount, maxCount);
      std::vector<TlHdf5Utils::size_type> coord(numOfLoadItems);
      for (size_type i = 0; i < numOfLoadItems; ++i) {
        coord[i] = cycle * bufferCount + i;
      }
      h5.getSelectedElements(h5path, coord, &buf);
      ++cycle;

      // 各プロセスのバッファに振り分ける
      std::vector<std::vector<TlMatrixObject::MatrixElement> > elements(
          numOfProcs);
      for (size_type i = 0; i < numOfLoadItems; ++i) {
        const int proc = this->getProcIdForIndex(row, col);
        elements[proc].push_back(
            TlMatrixObject::MatrixElement(row, col, buf[i]));

        {
          ++col;
          if (col > row) {
            col = 0;
            ++row;
          }
        }

        ++count;
        if (count >= maxCount) {
          isFinished = true;
          break;
        }
      }

      // データを送信
      for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
        if (isSendData[proc] == true) {
          rComm.wait(sizeLists[proc]);
          if (sizeLists[proc] > 0) {
            rComm.wait(&(elementsBuf[proc][0]));
          }
          isSendData[proc] = false;
        }

        sizeLists[proc] = elements[proc].size();
        elementsBuf[proc] = elements[proc];

        rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
        if (sizeLists[proc] > 0) {
          rComm.iSendDataX(&(elementsBuf[proc][0]), sizeLists[proc], proc,
                           TAG_LOAD_VALUES);
        }
        isSendData[proc] = true;
      }

      // proc=0分データの書き込み
      {
        const index_type sizeList = elements[0].size();
        for (index_type i = 0; i < sizeList; ++i) {
          const index_type r = elements[0][i].row;
          const index_type c = elements[0][i].col;
          const size_type index = this->getIndex(r, c);
          assert(index != -1);
          this->pData_[index] = elements[0][i].value;
        }
      }
    }  // end while

    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      if (isSendData[proc] == true) {
        rComm.wait(sizeLists[proc]);
        if (sizeLists[proc] > 0) {
          rComm.wait(&(elementsBuf[proc][0]));
        }
        isSendData[proc] = false;
      }
    }

    assert(rComm.checkNonBlockingCommunications());

    // 終了メッセージを全ノードに送る
    std::vector<int> endMsg(numOfProcs, 0);
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.iSendData(endMsg[proc], proc, TAG_LOAD_END);
    }
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.wait(endMsg[proc]);
    }
  } else {
    // slave
    const int root = 0;
    index_type sizeList = 0;
    std::vector<TlMatrixObject::MatrixElement> elements;
    int endMsg = 0;

    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
    rComm.iReceiveData(endMsg, root, TAG_LOAD_END);
    bool isLoopBreak = false;
    while (isLoopBreak == false) {
      if (rComm.test(sizeList) == true) {
        rComm.wait(sizeList);
        if (sizeList > 0) {
          elements.resize(sizeList);
          rComm.receiveDataX(&(elements[0]), sizeList, root, TAG_LOAD_VALUES);

          for (index_type i = 0; i < sizeList; ++i) {
            const index_type row = elements[i].row;
            const index_type col = elements[i].col;
            const size_type index = this->getIndex(row, col);
            assert(index != -1);
            this->pData_[index] = elements[i].value;
          }
        }
        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
      }

      if (rComm.test(endMsg) == true) {
        rComm.wait(endMsg);
        rComm.cancel(sizeList);
        isLoopBreak = true;
      }
    }
  }

  assert(rComm.checkNonBlockingCommunications());
  return answer;
}

#endif  // HAVE_HDF5

// local I/O -------------------------------------------------------------------
bool TlDenseSymmetricMatrix_blacs::saveLocal(
    const std::string& filePath) const {
  // std::cerr << "TlDenseSymmetricMatrix_blacs::saveLocal() file=" << filePath
  // << std::endl;

  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  // DataType& data_tmp = const_cast<DataType&>(this->data_);

  bool answer = true;
  std::ofstream ofs;
  ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

  const int type = 18;  // means Distributed(16) + RLHD(2)
  const index_type globalRow = this->m_nRows;
  const index_type globalCol = this->m_nCols;
  const index_type myRows = this->m_nMyRows;
  const index_type myCols = this->m_nCols;
  ofs.write(reinterpret_cast<const char*>(&type), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&globalRow), sizeof(index_type));
  ofs.write(reinterpret_cast<const char*>(&globalCol), sizeof(index_type));
  ofs.write(reinterpret_cast<const char*>(&myRows), sizeof(index_type));
  ofs.write(reinterpret_cast<const char*>(&myCols), sizeof(index_type));
  ofs.write(reinterpret_cast<const char*>(&this->m_nRank), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nProc), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nProcGridRow), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nProcGridCol), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nMyProcRow), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nMyProcCol), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(&this->m_nBlockSize), sizeof(int));
  ofs.write(reinterpret_cast<const char*>(this->pData_),
            sizeof(double) * this->getNumOfMyElements());

  ofs.close();

  // std::cerr << "save: rank=" << this->m_nRank << std::endl;
  // std::cerr << "save: proc=" << this->m_nProc << std::endl;
  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.barrier();

  return answer;
}

bool TlDenseSymmetricMatrix_blacs::loadLocal(const std::string& filePath) {
  // std::cerr << "TlDenseSymmetricMatrix_blacs::loadLocal() file=" << filePath
  // << std::endl;
  bool answer = false;

  std::ifstream ifs;
  ifs.open(filePath.c_str());
  if (ifs.fail() != true) {
    int type = 0;
    index_type globalRow = 0;
    index_type globalCol = 0;
    index_type myRows = 0;
    index_type myCols = 0;
    int rank = 0;
    int proc = 0;
    int procGridRow = 0;
    int procGridCol = 0;
    int myProcRow = 0;
    int myProcCol = 0;
    int blockSize = 0;
    ifs.read((char*)&type, sizeof(int));
    ifs.read((char*)&globalRow, sizeof(index_type));
    ifs.read((char*)&globalCol, sizeof(index_type));
    ifs.read((char*)&myRows, sizeof(index_type));
    ifs.read((char*)&myCols, sizeof(index_type));
    ifs.read((char*)&rank, sizeof(int));
    ifs.read((char*)&proc, sizeof(int));
    ifs.read((char*)&procGridRow, sizeof(int));
    ifs.read((char*)&procGridCol, sizeof(int));
    ifs.read((char*)&myProcRow, sizeof(int));
    ifs.read((char*)&myProcCol, sizeof(int));
    ifs.read((char*)&blockSize, sizeof(int));

    if (type != 18) {
      std::cerr << TlUtils::format("file type mismatch(%d). ", type) << __FILE__
                << "," << __LINE__ << std::endl;
    }
    this->m_nRows = globalRow;
    this->m_nCols = globalCol;
    // std::cerr << "rank = " << rank << std::endl;
    // std::cerr << "proc = " << proc << std::endl;
    assert(rank == this->m_nRank);
    assert(proc == this->m_nProc);
    assert(procGridRow == this->m_nProcGridRow);
    assert(procGridCol == this->m_nProcGridCol);
    assert(myProcRow == this->m_nMyProcRow);
    assert(myProcCol == this->m_nMyProcCol);
    assert(blockSize == this->m_nBlockSize);
    this->initialize();

    ifs.read((char*)this->pData_, sizeof(double) * this->getNumOfMyElements());

    answer = true;
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.barrier();

  return answer;
}

std::vector<TlMatrixObject::MatrixElement>
TlDenseSymmetricMatrix_blacs::getMatrixElementsInLocal() const {
  const std::size_t numOfMyElements = this->getNumOfMyElements();
  std::vector<TlMatrixObject::MatrixElement> answer;
  answer.reserve(numOfMyElements);

  const index_type numOfRowIndeces = this->m_RowIndexTable.size();
  const index_type numOfColIndeces = this->m_ColIndexTable.size();
  for (index_type i = 0; i < numOfRowIndeces; ++i) {
    const index_type globalRow = this->m_RowIndexTable[i];

    for (index_type j = 0; j < numOfColIndeces; ++j) {
      const index_type globalCol = this->m_ColIndexTable[j];

      if (globalRow >= globalCol) {
        answer.push_back(TlMatrixObject::MatrixElement(
            globalRow, globalCol, this->getLocal(globalRow, globalCol)));
      }
    }
  }

  return answer;
}

// =============================================================================
// matrix operation
// =============================================================================
bool TlDenseSymmetricMatrix_blacs::eig(
    TlDenseVector_Lapack* pEigVal, TlDenseGeneralMatrix_blacs* pEigVec,
    TlDenseSymmetricMatrix_blacs::DIAGONAL_METHOD method) const {
#ifdef HAVE_SCALAPACK
  if (method == DIVIDE_AND_CONQUER) {
    return diagonalByScaLapack_DC(*this, pEigVal, pEigVec);
  } else {
    return diagonalByScaLapack_QR(*this, pEigVal, pEigVec);
  }
#else
  std::cerr << "sorry. this code is not implemented." << std::endl;
  abort();
#endif  // HAVE_SCALAPACK
}

TlDenseSymmetricMatrix_blacs TlDenseSymmetricMatrix_blacs::inverse() const {
#ifdef HAVE_SCALAPACK
  // using SCALAPACK
  return inverseByScaLapack(*this);
#else
// without SCALAPACK
#error NOT found algebra package: need SCALAPACK library
#endif  // HAVE_SCALAPACK
}

// !!!
TlDenseGeneralMatrix_blacs multiplicationByScaLapack(
    const TlDenseSymmetricMatrix_blacs& X,
    const TlDenseGeneralMatrix_blacs& Y) {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  // TlDenseSymmetricMatrix_blacs& Xtmp =
  // const_cast<TlDenseSymmetricMatrix_blacs&>(X);  TlDenseGeneralMatrix_blacs&
  // Ytmp =
  // const_cast<TlDenseGeneralMatrix_blacs&>(Y);

  assert(X.getNumOfCols() == Y.getNumOfRows());

  // X.symmetrize();
  TlDenseGeneralMatrix_blacs Z(X.m_nRows, Y.m_nCols);

  // use SCALAPACK
  const char SIDE = 'L';  // L means "C := alpha*A*B + beta*C",
  // R means "C := alpha*B*A + beta*C"
  const char UPLO =
      'L';  // L means the lower triangular part of the symmetric matrix
  // U means the upper triangular part of the symmetric matrix
  const int M = Z.getNumOfRows();  // the number of rows of the matrix  C
  const int N = Z.getNumOfCols();  // the number of columns of the matrix C
  const double ALPHA = 1.0;        // ALPHA specifies the scalar alpha

  const double* A = X.pData_;

  const int IA = 1;
  const int JA = 1;

  // const double* B =
  // &(const_cast<TlDenseGeneralMatrix_blacs&>(Y).m_aMatrix[0]); //
  // DIMENSION (LDB, n)
  const double* B = Y.pData_;
  const int IB = 1;
  const int JB = 1;

  const double BETA = 0.0;  // BETA  specifies the scalar  beta
  double* C = Z.pData_;     // DIMENSION (LDC, n)
  const int IC = 1;
  const int JC = 1;

  pdsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &IA, &JA, X.m_pDESC, B, &IB, &JB,
          Y.m_pDESC, &BETA, C, &IC, &JC, Z.m_pDESC);

  return Z;
}

TlDenseGeneralMatrix_blacs multiplicationByScaLapack(
    const TlDenseGeneralMatrix_blacs& X,
    const TlDenseSymmetricMatrix_blacs& Y) {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  // TlDenseGeneralMatrix_blacs& Xtmp =
  // const_cast<TlDenseGeneralMatrix_blacs&>(X);
  // TlDenseSymmetricMatrix_blacs& Ytmp =
  // const_cast<TlDenseSymmetricMatrix_blacs&>(Y);

  assert(X.getNumOfCols() == Y.getNumOfRows());

  // Y.symmetrize();
  TlDenseGeneralMatrix_blacs Z(X.m_nRows, Y.m_nCols);

  // use SCALAPACK
  const char SIDE = 'R';  // L means "C := alpha*A*B + beta*C",
                          // R means "C := alpha*B*A + beta*C"
  const char UPLO =
      'L';  // L means the lower triangular part of the symmetric matrix
            // U means the upper triangular part of the symmetric matrix
  const int M = Z.getNumOfRows();  // the number of rows of the matrix  C
  const int N = Z.getNumOfCols();  // the number of columns of the matrix C
  const double ALPHA = 1.0;        // ALPHA specifies the scalar alpha

  // const double* A =
  // &(const_cast<TlDenseGeneralMatrix_blacs&>(X).m_aMatrix[0]); //
  // DIMENSION (LDA, ka)
  const double* A = Y.pData_;
  const int IA = 1;
  const int JA = 1;

  // const double* B =
  // &(const_cast<TlDenseSymmetricMatrix_blacs&>(Y).m_aMatrix[0]);    //
  // DIMENSION (LDB, n)
  const double* B = X.pData_;  // DIMENSION (LDB, n)
  const int IB = 1;
  const int JB = 1;

  const double BETA = 0.0;  // BETA  specifies the scalar  beta
  double* C = Z.pData_;     // DIMENSION (LDC, n)
  const int IC = 1;
  const int JC = 1;

  pdsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &IA, &JA, Y.m_pDESC, B, &IB, &JB,
          X.m_pDESC, &BETA, C, &IC, &JC, Z.m_pDESC);

  return Z;
}

bool diagonalByScaLapack_QR(const TlDenseSymmetricMatrix_blacs& inMatrix,
                            TlDenseVector_Lapack* outEigVal,
                            TlDenseGeneralMatrix_blacs* outEigVec) {
  assert(outEigVal != NULL);
  assert(outEigVec != NULL);

  const char* JOBZ = "V";
  const char* UPLO = "L";

  assert(inMatrix.m_nRows == inMatrix.m_nCols);
  const int N = inMatrix.m_nRows;

  double* pBufA = inMatrix.pData_;

  const int IA = 1;
  const int JA = 1;
  const int* DESCA = inMatrix.m_pDESC;

  outEigVal->resize(N);
  double* W = outEigVal->data();

  outEigVec->resize(N, N);
  double* Z = outEigVec->pData_;

  const int IZ = 1;
  const int JZ = 1;
  const int* DESCZ = outEigVec->m_pDESC;

  int LWORK = -1;
  double* pWorkSizeCheck = new double[1];
  int info = 0;

  pdsyev_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ,
          pWorkSizeCheck, &LWORK, &info);

  if (info != 0) {
    std::cout << "pdsyev_ error @1st call: " << info << std::endl;
    return false;
  }

  LWORK = (int)pWorkSizeCheck[0];
  double* pWork = new double[LWORK];

  delete[] pWorkSizeCheck;
  pWorkSizeCheck = NULL;

  pdsyev_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ, pWork,
          &LWORK, &info);

  delete[] pWork;
  pWork = NULL;

  if (info != 0) {
    std::cout << "pdsyev_ error @2nd call: " << info << std::endl;
  }
  return ((info == 0) ? true : false);
}

// Divide-and-Conquer Algorithm
bool diagonalByScaLapack_DC(const TlDenseSymmetricMatrix_blacs& inMatrix,
                            TlDenseVector_Lapack* outEigVal,
                            TlDenseGeneralMatrix_blacs* outEigVec) {
  assert(outEigVal != NULL);
  assert(outEigVec != NULL);

  const char* JOBZ = "V";
  const char* UPLO = "L";

  assert(inMatrix.m_nRows == inMatrix.m_nCols);
  const int N = inMatrix.m_nRows;

  double* pBufA = inMatrix.pData_;

  const int IA = 1;
  const int JA = 1;
  const int* DESCA = inMatrix.m_pDESC;

  outEigVal->resize(N);
  double* W = outEigVal->data();

  outEigVec->resize(N, N);
  double* Z = outEigVec->pData_;

  const int IZ = 1;
  const int JZ = 1;
  const int* DESCZ = outEigVec->m_pDESC;

  int LWORK = -1;
  double* pWorkSizeCheck = new double[1];

  const int NPCOL = inMatrix.m_nProcGridCol;
  const int LIWORK = 7 * N + 8 * NPCOL + 2;
  int* IWORK = new int[LIWORK];
  int info = 0;

  pdsyevd_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ,
           pWorkSizeCheck, &LWORK, IWORK, &LIWORK, &info);

  if (info != 0) {
    std::cout << "pdsyev_ error @1st call: " << info << std::endl;
    return false;
  }

  LWORK = (int)pWorkSizeCheck[0];
  double* pWork = new double[LWORK];

  delete[] pWorkSizeCheck;
  pWorkSizeCheck = NULL;

  pdsyevd_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ, pWork,
           &LWORK, IWORK, &LIWORK, &info);

  delete[] IWORK;
  IWORK = NULL;

  delete[] pWork;
  pWork = NULL;

  if (info != 0) {
    std::cout << "pdsyev_ error @2nd call: " << info << std::endl;
  }
  return ((info == 0) ? true : false);
}

TlDenseSymmetricMatrix_blacs inverseByScaLapack(const TlDenseSymmetricMatrix_blacs& X) {
  TlDenseGeneralMatrix_blacs Y = X;
  TlDenseSymmetricMatrix_blacs answer = Y.inverse();
  return answer;
  // 以下のルーチンは正の対称行列でなければならない
  //   const int N = X.getNumOfRows();
  //   const int IA = 1;
  //   const int JA = 1;
  //   int info = 0;

  //   pdpotrf_("U", &N, &(X.data_[0]), &IA, &JA, X.m_pDESC, &info);
  //   if (info != 0){
  //     std::cout << "pdpotrf_ returns " << info << std::endl;
  //     return false;
  //   }

  //   pdpotri_("U", &N, &(X.data_[0]), &IA, &JA, X.m_pDESC, &info);
  //   if (info != 0){
  //     std::cout << "pdpotri_ returns " << info << std::endl;
  //     return false;
  //   }
}

// void TlDenseSymmetricMatrix_blacs::symmetrize() const {
//   const int nDim = this->m_nRows;
//   for (int row = 0; row < nDim; ++row){
//     const int nGlobalRow = row +1;
//     for (int col = 0; col < row; ++col){
//       assert(row > col);
//       const double dValue = this->get(row, col);

//       // copy (row, col) to (col, row)
//       const int nGlobalCol = col +1;
//       pdelset_(&(this->data_[0]), &nGlobalCol, &nGlobalRow, this->m_pDESC,
//       &dValue);
//     }
//   }
// }

TlDenseGeneralMatrix_blacs operator*(const TlDenseSymmetricMatrix_blacs& X,
                                     const TlDenseSymmetricMatrix_blacs& Y) {
  TlDenseGeneralMatrix_blacs Z = Y;
  return multiplicationByScaLapack(X, Z);
}

TlDenseGeneralMatrix_blacs operator*(const TlDenseSymmetricMatrix_blacs& X,
                                     const TlDenseGeneralMatrix_blacs& Y) {
  return multiplicationByScaLapack(X, Y);
}

TlDenseGeneralMatrix_blacs operator*(const TlDenseGeneralMatrix_blacs& X,
                                     const TlDenseSymmetricMatrix_blacs& Y) {
  return multiplicationByScaLapack(X, Y);
}

TlDenseSymmetricMatrix_blacs operator*(const double X,
                                       const TlDenseSymmetricMatrix_blacs& Y) {
  TlDenseSymmetricMatrix_blacs ans = Y;
  ans *= X;

  return ans;
}

TlDistributedVector operator*(const TlDenseSymmetricMatrix_blacs& A,
                              const TlDistributedVector& X) {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  TlDenseSymmetricMatrix_blacs& Atmp =
      const_cast<TlDenseSymmetricMatrix_blacs&>(A);
  TlDistributedVector& Xtmp = const_cast<TlDistributedVector&>(X);

  const int N = A.getNumOfCols();
  assert(N == X.getSize());

  TlDistributedVector Y(N);

  const char UPLO = 'L';
  const double alpha = 1.0;
  const double beta = 0.0;
  const int IA = 1;
  const int JA = 1;
  const int IX = 1;
  const int JX = 1;
  const int INCX = 1;
  const int IY = 1;
  const int JY = 1;
  const int INCY = 1;

  pdsymv_(&UPLO, &N, &alpha, Atmp.pData_, &IA, &JA, A.m_pDESC, &(Xtmp.data_[0]),
          &IX, &JX, X.m_pDESC, &INCX, &beta, &(Y.data_[0]), &IY, &JY, Y.m_pDESC,
          &INCY);

  return Y;
}

TlDenseSymmetricMatrix_blacs operator*(const TlDenseSymmetricMatrix_blacs& X,
                                       const double Y) {
  TlDenseSymmetricMatrix_blacs ans = X;
  ans *= Y;

  return ans;
}

// // Harbrecht, Peter, Schneider, 2011
// TlDenseGeneralMatrix_blacs
// TlDenseSymmetricMatrix_blacs::choleskyFactorization(
//     const double threshold) const {
//   // TlLogging& log = TlLogging::getInstance();
//
//   const index_type N = this->getNumOfRows();
//   TlDenseVector_Lapack d = this->getDiagonalElements();
//   double error = d.sum();
//   std::vector<TlDenseVector_Lapack::size_type> pivot(N);
//   for (index_type i = 0; i < N; ++i) {
//     pivot[i] = i;
//   }
//
//   TlDenseGeneralMatrix_blacs L(N, N);
//   index_type m = 0;
//   // double sum_ll = 0.0;
//
//   while (error > threshold) {
//     std::vector<TlDenseVector_Lapack::size_type>::const_iterator it =
//         d.argmax(pivot.begin() + m, pivot.end());
//     const index_type i = it - pivot.begin();
//     std::swap(pivot[m], pivot[i]);
//
//     const double l_m_pm = std::sqrt(d[pivot[m]]);
//     L.set(m, pivot[m], l_m_pm);
//
//     const double inv_l_m_pm = 1.0 / l_m_pm;
//
//     for (index_type i = m + 1; i < N; ++i) {
//       double sum_ll = 0.0;
//       for (index_type j = 0; j < m; ++j) {
//         sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
//       }
//       const double l_m_pi =
//           (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
//       L.set(m, pivot[i], l_m_pi);
//
//       d[pivot[i]] -= l_m_pi * l_m_pi;
//     }
//
//     error = 0.0;
//     for (index_type i = m + 1; i < N; ++i) {
//       error += d[pivot[i]];
//     }
//
//     // log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
//     ++m;
//   }
//
//   L.transposeInPlace();
//   L.resize(N, m);
//
//   return L;
// }
//
// TlDenseGeneralMatrix_blacs
// TlDenseSymmetricMatrix_blacs::choleskyFactorization_mod(
//     const double threshold) const {
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   TlLogging& log = TlLogging::getInstance();
//
//   const index_type N = this->getNumOfRows();
//   TlDenseVector_Lapack d = this->getDiagonalElements();
//   double error = d.sum();
//   std::vector<TlDenseVector_Lapack::size_type> pivot(N);
//   for (index_type i = 0; i < N; ++i) {
//     pivot[i] = i;
//   }
//
//   TlDenseGeneralMatrix_blacs L(N, N);
//   index_type m = 0;
//   // double sum_ll = 0.0;
//
//   while (error > threshold) {
//     std::vector<TlDenseVector_Lapack::size_type>::const_iterator it =
//         d.argmax(pivot.begin() + m, pivot.end());
//     const index_type i = it - pivot.begin();
//     std::swap(pivot[m], pivot[i]);
//
//     const double l_m_pm = std::sqrt(d[pivot[m]]);
//     L.set(m, pivot[m], l_m_pm);
//
//     const double inv_l_m_pm = 1.0 / l_m_pm;
//
//     // L(0:m, pivot[m])を一時保管
//     const TlDenseVector_Lapack L_x_pm =
//         L.getColVector(pivot[m]);  // TODO: バッファの取り過ぎ
//
//     for (index_type i = m + 1; i < N; ++i) {
//       double sum_ll = 0.0;
//       for (index_type j = 0; j < m; ++j) {
//         // sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
//         sum_ll += L_x_pm[j] * L.getLocal(j, pivot[i]);
//       }
//       rComm.allReduce_SUM(sum_ll);
//
//       double l_m_pi = 0.0;
//       {
//         // double G_pm_pi_exact = this->get(pivot[m], pivot[i]);
//         // double l_m_pi_exact = (this->get(pivot[m], pivot[i]) - sum_ll) *
//         // inv_l_m_pm; double G_pm_pi_local = this->getLocal(pivot[m],
//         // pivot[i]);
//
//         // for (int p = 0; p < rComm.getNumOfProcs(); ++p) {
//         //     if (p == rComm.getRank()) {
//         //         std::cerr << TlUtils::format("[%d] G(%d, %d)=% f(% f),
//         //         holder=%d, index=%d",
//         //                                      rComm.getRank(),
//         //                                      pivot[m], pivot[i],
//         //                                      G_pm_pi_exact, G_pm_pi_local,
//         // this->getProcIdForIndex(pivot[m],
//         //                                      pivot[i]),
//         //                                      this->getIndex(pivot[m],
//         //                                      pivot[i]))
//         //                   << std::endl;
//         //     }
//         //     rComm.barrier();
//         // }
//
//         int holder = this->getProcIdForIndex(pivot[m], pivot[i]);
//         if (holder == rComm.getRank()) {
//           assert(this->getIndex(pivot[m], pivot[i]) != -1);
//           // assert(std::fabs(this->getLocal(pivot[m], pivot[i]) -
//           // G_pm_pi_exact) < 1.0E-5);
//           l_m_pi = (this->getLocal(pivot[m], pivot[i]) - sum_ll) *
//           inv_l_m_pm;
//         }
//         rComm.broadcast(l_m_pi, holder);
//       }
//       // assert(std::fabs(l_m_pi - l_m_pi_exact) < 1.0E-5);
//       L.set(m, pivot[i], l_m_pi);
//
//       d[pivot[i]] -= l_m_pi * l_m_pi;
//     }
//
//     error = 0.0;
//     for (index_type i = m + 1; i < N; ++i) {
//       error += d[pivot[i]];
//     }
//
//     log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
//     ++m;
//   }
//
//   L.transposeInPlace();
//   L.resize(N, m);
//
//   return L;
// }
//
TlDenseGeneralMatrix_blacs
TlDenseSymmetricMatrix_blacs::choleskyFactorization_mod2(
    const double threshold) const {
  // timing data
  TlTime CD_all_time;
  TlTime CD_rowvec_time;
  TlTime CD_bcast_time;
  TlTime CD_calc_time;
  TlTime CD_allreduce_time;

  CD_all_time.start();
  this->log_.info("TlDenseSymmetricMatrix_blacs::choleskyFactorization_mod2()");
  const bool isEnableMmap = false;  // mmapをつかうかどうか
  TlCommunicate& rComm = TlCommunicate::getInstance();

  // prepare variables
  const index_type N = this->getNumOfRows();
  assert(N == this->getNumOfCols());
  TlDenseVector_Lapack global_diagonals =
      this->getDiagonalElements();  // 対角成分
  TlDenseGeneralMatrix_arrays_RowOriented L(
      N, 1, rComm.getNumOfProcs(), rComm.getRank(),
      isEnableMmap);  // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
  const index_type local_N = L.getNumOfLocalVectors();
  TlDenseVector_Lapack L_pm(N);
  std::vector<int> global_pivot(N);   //
  std::vector<int> reverse_pivot(N);  // global_pivotの逆引き
  std::vector<int> local_pivot(local_N);
  std::vector<double> local_diagonals(local_N);
  const int myRank = rComm.getRank();

  double error = 0.0;
  index_type error_global_loc = 0;
  index_type error_local_loc = 0;
  {
    int local_i = 0;
    for (int global_i = 0; global_i < N; ++global_i) {
      global_pivot[global_i] = global_i;
      reverse_pivot[global_i] = global_i;
      if (L.getSubunitID(global_i) == myRank) {
        local_pivot[local_i] = global_i;
        local_diagonals[local_i] = global_diagonals.get(global_i);
        if (error < local_diagonals[local_i]) {
          error_local_loc = local_i;
          error_global_loc = local_pivot[local_i];
          error = local_diagonals[local_i];
        }
        ++local_i;
        // this->log_.warn(TlUtils::format("%d/%d", local_i, local_N));
      }
    }
    assert(local_i == local_N);
  }
  rComm.allReduce_MAXLOC(&error, &error_global_loc);

  index_type m = 0;
  index_type local_m = 0;
  if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
    std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
    ++local_m;
  }

  int progress = 0;
  index_type division = index_type(N * 0.1);
  while (error > threshold) {
    // progress
    if (m >= progress * division) {
      this->log_.info(
          TlUtils::format("cd progress %8d/%8d, error=%f", m, N, error));
      ++progress;
      L.reserveColSize(progress * division);  // メモリの確保
    }
    L.resize(N, m + 1);

    // pivot
    std::swap(global_pivot[m], global_pivot[error_global_loc]);
    reverse_pivot[global_pivot[m]] = m;
    reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

    const double l_m_pm = std::sqrt(error);
    const index_type pivot_m = global_pivot[m];
    L.set(pivot_m, m, l_m_pm);  // 通信発生せず。関係無いPEは値を捨てる。
    const double inv_l_m_pm = 1.0 / l_m_pm;

    // get row vector
    CD_rowvec_time.start();
    const index_type numOf_G_cols = N - (m + 1);
    std::vector<double> G_pm(N - (m + 1));
    // const index_type numOf_G_cols = N -(m+1);
    {
      const TlDenseVector_Lapack Gm = this->getRowVector(pivot_m);
      assert(Gm.getSize() == N);
      for (index_type i = 0; i < numOf_G_cols; ++i) {
        const index_type pivot_i = global_pivot[m + 1 + i];  // from (m+1) to N
        assert(0 <= pivot_i);
        assert(pivot_i < N);
        G_pm[i] = Gm.get(pivot_i);
      }
    }
    CD_rowvec_time.stop();

    CD_bcast_time.start();
    {
      // 全PEに分配
      const int PEinCharge = L.getSubunitID(pivot_m);
      if (PEinCharge == rComm.getRank()) {
        // const index_type copySize = L.getRowVector(pivot_m, &(L_pm[0]), m
        // +1);
        L_pm = L.getVector(pivot_m);
        assert(L_pm.getSize() == m + 1);
      }
      // rComm.broadcast(&(L_pm[0]), m +1, PEinCharge);
      rComm.broadcast(&L_pm, PEinCharge);
    }
    CD_bcast_time.stop();

    CD_calc_time.start();
    error = 0.0;
#pragma omp parallel
    {
      TlDenseVector_Lapack L_pi(m + 1);
      double my_error = 0.0;
      int my_error_global_loc = 0;
      int my_error_local_loc = 0;

#pragma omp for schedule(runtime)
      for (int i = local_m; i < local_N; ++i) {
        const int pivot_i = local_pivot[i];
        // const index_type copySize = L.getRowVector(pivot_i, &(L_pi[0]), m
        // +1);
        L_pi = L.getVector(pivot_i);
        assert(L_pi.getSize() == m + 1);
        double sum_ll = 0.0;
        for (index_type j = 0; j < m; ++j) {
          sum_ll += L_pm.get(j) * L_pi.get(j);
        }

        const int G_pm_index = reverse_pivot[pivot_i] - (m + 1);
        const double l_m_pi = (G_pm[G_pm_index] - sum_ll) * inv_l_m_pm;
        const double ll = l_m_pi * l_m_pi;
#pragma omp critical(DfCD_Parallel__calcCholeskyVectors_onTheFly)
        {
          L.set(pivot_i, m, l_m_pi);
          global_diagonals.add(pivot_i, -ll);
        }

        if (global_diagonals.get(pivot_i) > my_error) {
          my_error = global_diagonals.get(pivot_i);
          my_error_global_loc = reverse_pivot[pivot_i];  // == m +1 + i
          my_error_local_loc = i;
        }
      }

#ifdef _OPENMP
      const int numOfThreads = omp_get_num_threads();
      const int myThreadID = omp_get_thread_num();
      for (int thread = 0; thread < numOfThreads; ++thread) {
        if (thread == myThreadID) {
          if (error < my_error) {
            error = my_error;
            error_global_loc = my_error_global_loc;
            error_local_loc = my_error_local_loc;
          }
        }
#pragma omp flush(error, error_global_loc)
      }
#else
      error = my_error;
      error_global_loc = my_error_global_loc;
      error_local_loc = my_error_local_loc;
#endif  // _OPENMP
    }
    CD_calc_time.stop();

    CD_allreduce_time.start();
    rComm.allReduce_MAXLOC(&error, &error_global_loc);
    global_diagonals.set(global_pivot[error_global_loc], error);

    ++m;
    if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
      std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
      ++local_m;
    }
    CD_allreduce_time.stop();
  }
  this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));
  CD_all_time.stop();

  // timing data
  this->log_.info(TlUtils::format("CD all:       %10.1f sec.",
                                  CD_all_time.getElapseTime()));
  this->log_.info(TlUtils::format("CD rowvec:    %10.1f sec.",
                                  CD_rowvec_time.getElapseTime()));
  this->log_.info(TlUtils::format("CD bacst:     %10.1f sec.",
                                  CD_bcast_time.getElapseTime()));
  this->log_.info(TlUtils::format("CD calc:      %10.1f sec.",
                                  CD_calc_time.getElapseTime()));
  this->log_.info(TlUtils::format("CD allreduce: %10.1f sec.",
                                  CD_allreduce_time.getElapseTime()));

  return TlDenseGeneralMatrix_blacs(L);
}
