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

#include "tl_dense_vector_blacs.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include "TlLogging.h"
#include "scalapack.h"
#include "tl_dense_general_matrix_blacs.h"
#include "tl_dense_vector_lapack.h"
#include "tl_scalapack_context.h"

int TlDistributedVector::systemBlockSize_ = 64;
const std::size_t TlDistributedVector::FILE_BUFFER_SIZE =
    100 * 1024 * 1024;  // 100 MB

TlDistributedVector::TlDistributedVector(const TlVectorAbstract::index_type nSize)
    : m_nContext(0),
      m_nRows(nSize),
      m_nCols(1),
      m_nBlockSize(TlDistributedVector::systemBlockSize_) {
  this->initialize();
}

TlDistributedVector::TlDistributedVector(const TlDistributedVector& rhs)
    : m_nContext(0),
      m_nRows(rhs.m_nRows),
      m_nCols(rhs.m_nCols),
      m_nBlockSize(rhs.m_nBlockSize) {
  this->initialize();
  this->data_ = rhs.data_;
}

TlDistributedVector::TlDistributedVector(
    const std::vector<double>& rhs, const TlVectorAbstract::size_type globalSize)
    : m_nContext(0),
      m_nRows(globalSize),
      m_nCols(1),
      m_nBlockSize(TlDistributedVector::systemBlockSize_) {
  this->initialize();
  const size_t copySize = std::min(rhs.size(), this->data_.size());
  std::memcpy(&(this->data_[0]), &(rhs[0]), copySize);
}

TlDistributedVector::TlDistributedVector(const TlDenseVector_Lapack& rhs)
    : m_nContext(0),
      m_nRows(rhs.getSize()),
      m_nCols(1),
      m_nBlockSize(TlDistributedVector::systemBlockSize_) {
  this->initialize();

  size_type size = rhs.getSize();
  for (size_type i = 0; i < size; ++i) {
    const int index = this->getIndex(i, 0);
    if (index != -1) {
      this->data_[index] = rhs.get(i);
    }
  }
}

TlDistributedVector::~TlDistributedVector() {}

void TlDistributedVector::setSystemBlockSize(int blockSize) {
  TlDistributedVector::systemBlockSize_ = blockSize;
}

void TlDistributedVector::initialize() {
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

  this->data_.resize((this->m_nMyRows * this->m_nMyCols), 0.0);

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

void TlDistributedVector::resize(const TlVectorAbstract::index_type size) {
  assert(size > 0);

  TlDistributedVector tmp = *this;
  this->m_nRows = size;
  this->m_nCols = 1;
  this->initialize();

  const TlVectorAbstract::index_type copySize = std::min(size, tmp.getSize());
  for (TlVectorAbstract::index_type i = 0; i < copySize; ++i) {
    this->set(i, tmp.get(i));
  }
}

TlDistributedVector& TlDistributedVector::operator=(
    const TlDistributedVector& rhs) {
  if (&rhs != this) {
    this->m_nRows = rhs.m_nRows;
    this->m_nCols = rhs.m_nCols;
    this->m_nBlockSize = rhs.m_nBlockSize;
    this->initialize();
    this->data_ = rhs.data_;
  }

  return *this;
}

TlDistributedVector& TlDistributedVector::operator=(const TlDenseVector_Lapack& rhs) {
  this->m_nRows = rhs.getSize();
  this->m_nCols = 1;
  this->initialize();

  size_type size = rhs.getSize();
  for (size_type i = 0; i < size; ++i) {
    const int index = this->getIndex(i, 0);
    if (index != -1) {
      this->data_[index] = rhs.get(i);
    }
  }

  return *this;
}

TlDistributedVector::size_type TlDistributedVector::getSize() const {
  return this->m_nRows;
}

int TlDistributedVector::getIndex(int nGlobalRow, int nGlobalCol) const {
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

void TlDistributedVector::set(const size_type i, const double dValue) {
  const int index = this->getIndex(i, 0);
  if (index != -1) {
    this->data_[index] = dValue;
  }
}

double TlDistributedVector::get(const size_type index) const {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  DataType& dataTmp = const_cast<DataType&>(this->data_);

  double dAnswer = 0.0;
  const int nGlobalRow = index + 1;
  const int nGlobalCol = 0 + 1;
  pdelget_("A", " ", &dAnswer, &(dataTmp[0]), &nGlobalRow, &nGlobalCol,
           this->m_pDESC);

  return dAnswer;
}

void TlDistributedVector::add(const size_type index, const double value) {
  const int localIndex = this->getIndex(index, 0);
  if (localIndex != -1) {
    this->data_[localIndex] = value;
  }
}

double TlDistributedVector::operator[](const size_type index) const {
  return this->get(index);
}

double& TlDistributedVector::operator[](const size_type index) {
  const int nGlobalRow = index + 1;
  const int nGlobalCol = 0 + 1;

  int nLocalRow = 0;
  int nLocalCol = 0;
  int nLocalProcRow = 0;
  int nLocalProcCol = 0;
  infog2l_(&nGlobalRow, &nGlobalCol, this->m_pDESC, &(this->m_nProcGridRow),
           &(this->m_nProcGridCol), &(this->m_nMyProcRow),
           &(this->m_nMyProcCol), &nLocalRow, &nLocalCol, &nLocalProcRow,
           &nLocalProcCol);

  this->m_dTempVar = this->get(index);

  if ((this->m_nMyProcRow == nLocalProcRow) &&
      (this->m_nMyProcCol == nLocalProcCol)) {
    // const size_t index = (nLocalRow -1) +(nLocalCol -1)*(this->m_nMyCols);
    const size_t index = (nLocalRow - 1) + (nLocalCol - 1) * (this->m_nMyRows);
    return this->data_[index];
  } else {
    return this->m_dTempVar;  // 使い捨て
  }
}

int TlDistributedVector::getProcIdForIndex(const index_type globalIndex) const {
  const int rowBlockIndex = globalIndex / this->m_nBlockSize;
  const int rowProcID = rowBlockIndex % this->m_nProcGridRow;
  // const int colBlockIndex = 0;
  const int colProcID = 0;

  // "Row-major" for Cblacs_gridinit
  const int procMatrixIndex = rowProcID * this->m_nProcGridCol + colProcID;
  // const int procID = this->processMatrix_[procMatrixIndex];
  int procID = procMatrixIndex;

  return procID;
}

TlDistributedVector& TlDistributedVector::operator+=(
    const TlDistributedVector& rhs) {
  assert(this->m_nRows == rhs.m_nRows);
  assert(this->m_nCols == rhs.m_nCols);
  assert(this->m_nBlockSize == rhs.m_nBlockSize);

  this->data_ += rhs.data_;

  return (*this);
}

TlDistributedVector& TlDistributedVector::operator-=(
    const TlDistributedVector& rhs) {
  assert(this->m_nRows == rhs.m_nRows);
  assert(this->m_nCols == rhs.m_nCols);
  assert(this->m_nBlockSize == rhs.m_nBlockSize);

  this->data_ -= rhs.data_;

  return (*this);
}

TlDistributedVector& TlDistributedVector::operator*=(const double dCoef) {
  this->data_ *= dCoef;

  return (*this);
}

TlDistributedVector& TlDistributedVector::operator/=(const double dCoef) {
  this->data_ /= dCoef;

  return (*this);
}

TlDistributedVector operator+(const TlDistributedVector& X,
                              const TlDistributedVector& Y) {
  TlDistributedVector Z(X);
  Z += Y;

  return Z;
}

TlDistributedVector operator-(const TlDistributedVector& X,
                              const TlDistributedVector& Y) {
  TlDistributedVector Z(X);
  Z -= Y;

  return Z;
}

bool TlDistributedVector::load(const std::string& sFilePath) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  bool bAnswer = false;

  std::ifstream ifs;
  if (rComm.isMaster() == true) {
    ifs.open(sFilePath.c_str());
  }

  bool bIsFail = false;
  if (rComm.isMaster() == true) {
    bIsFail = ifs.fail();
  }
  rComm.broadcast(bIsFail);

  if (bIsFail) {
    std::cerr << "[error] TlDistributedVector::load(): could not open file. "
              << sFilePath << std::endl;
    abort();
  }

  bAnswer = this->load(ifs);

  if (bAnswer != true) {
    std::cerr << "TlDenseGeneralMatrix_blacs::load() is not supported: "
              << sFilePath << std::endl;
  }

  if (rComm.isMaster() == true) {
    ifs.close();
  }

  return bAnswer;
}

bool TlDistributedVector::load(std::ifstream& ifs) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  TlLogging& log = TlLogging::getInstance();

  const int numOfProcs = rComm.getNumOfProc();

  bool bAnswer = true;

  // binary mode
  // read header
  int nSize = 0;
  if (rComm.isMaster() == true) {
    ifs.read((char*)&(nSize), sizeof(int));
  }
  rComm.broadcast(nSize);
  this->resize(nSize);  // initialize?

  if (rComm.isMaster() == true) {
    static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
    std::vector<double> buf(bufferCount, 0.0);

    index_type currentIndex = 0;
    bool isFinished = false;

    std::vector<std::size_t> sizeLists(numOfProcs, 0);
    std::vector<std::vector<index_type> > indexLists(numOfProcs,
                                                     std::vector<index_type>());
    std::vector<std::vector<double> > valueLists(numOfProcs,
                                                 std::vector<double>());
    std::vector<bool> isSendData(numOfProcs, false);

    while (isFinished == false) {
      // buffer分を一度に読み込み
      ifs.read((char*)(&buf[0]), sizeof(double) * bufferCount);

      // 各プロセスのバッファに振り分ける
      std::map<int, std::vector<index_type> > tmpIndexLists;
      std::map<int, std::vector<double> > tmpValueLists;
      for (std::size_t i = 0; i < bufferCount; ++i) {
        const int proc = this->getProcIdForIndex(currentIndex);

        if (proc == 0) {
          // masterが持っている
          const size_type index = this->getIndex(currentIndex, 0);
          assert(index != -1);
          this->data_[index] = buf[i];
        } else {
          log.debug(TlUtils::format("SEND [%d] (%4d)", proc, currentIndex));

          tmpIndexLists[proc].push_back(currentIndex);
          tmpValueLists[proc].push_back(buf[i]);
        }

        // count up
        ++currentIndex;
        if (currentIndex >= this->getSize()) {
          isFinished = true;
          break;
        }
      }

      // データを送信
      std::map<int, std::vector<index_type> >::const_iterator itEnd =
          tmpIndexLists.end();
      for (std::map<int, std::vector<index_type> >::const_iterator it =
               tmpIndexLists.begin();
           it != itEnd; ++it) {
        const int proc = it->first;
        const std::size_t numOfContents = it->second.size();
        assert(numOfContents == tmpValueLists[proc].size());

        if (numOfContents == 0) {
          // 送る要素がない場合は送らない
          continue;
        }

        if (isSendData[proc] == true) {
          rComm.wait(sizeLists[proc]);
          rComm.wait(&(indexLists[proc][0]));
          rComm.wait(&(valueLists[proc][0]));
          isSendData[proc] = false;
        }

        sizeLists[proc] = numOfContents;
        indexLists[proc] = tmpIndexLists[proc];
        valueLists[proc] = tmpValueLists[proc];

        rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
        rComm.iSendDataX(&(indexLists[proc][0]), sizeLists[proc], proc,
                         TAG_LOAD_ROWCOLS);
        rComm.iSendDataX(&(valueLists[proc][0]), sizeLists[proc], proc,
                         TAG_LOAD_VALUES);
        isSendData[proc] = true;
      }
    }  // end while

    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      if (isSendData[proc] == true) {
        rComm.wait(sizeLists[proc]);
        rComm.wait(&(indexLists[proc][0]));
        rComm.wait(&(valueLists[proc][0]));
        isSendData[proc] = false;
      }
    }

    // 終了メッセージを全ノードに送る
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      if (isSendData[proc] == true) {
        rComm.wait(sizeLists[proc]);
        rComm.wait(&(indexLists[proc][0]));
        rComm.wait(&(valueLists[proc][0]));
      }
    }
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
    int sizeList = 0;
    std::vector<index_type> indexList;
    std::vector<double> valueList;
    int endMsg = 0;

    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
    rComm.iReceiveData(endMsg, root, TAG_LOAD_END);
    bool isLoopBreak = false;
    while (isLoopBreak == false) {
      if (rComm.test(sizeList) == true) {
        rComm.wait(sizeList);
        log.debug(
            TlUtils::format("RECV [%d] size=%d", rComm.getRank(), sizeList));

        assert(sizeList != -1);
        indexList.resize(sizeList);
        valueList.resize(sizeList);
        rComm.iReceiveDataX(&(indexList[0]), sizeList, root, TAG_LOAD_ROWCOLS);
        rComm.iReceiveDataX(&(valueList[0]), sizeList, root, TAG_LOAD_VALUES);
        rComm.wait(&(indexList[0]));
        rComm.wait(&(valueList[0]));
        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);

        for (int i = 0; i < sizeList; ++i) {
          const index_type row = indexList[i];
          const size_type index = this->getIndex(row, 0);

          log.debug(
              TlUtils::format("RECV [%d] (%4d, %4d)", rComm.getRank(), row));

          assert(index != -1);
          this->data_[index] = valueList[i];
        }
      }

      if (rComm.test(endMsg) == true) {
        rComm.wait(endMsg);
        rComm.cancel(sizeList);
        log.debug(TlUtils::format("RECV [%d] END", rComm.getRank()));
        isLoopBreak = true;
      }
    }
  }

  rComm.barrier(true);
  assert(rComm.checkNonBlockingCommunications());
  return bAnswer;
}

bool TlDistributedVector::save(const std::string& sFilePath) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  // const double dEndTime = rComm.getTime();

  // this->update();
  bool bAnswer = true;

  std::ofstream ofs;
  if (rComm.isMaster() == true) {
    ofs.open(sFilePath.c_str(), std::ofstream::out | std::ofstream::binary);
  }

  bAnswer = this->save(ofs);

  if (rComm.isMaster() == true) {
    ofs.close();
  }

  return bAnswer;
}

bool TlDistributedVector::save(std::ofstream& ofs) const {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  DataType& dataTmp = const_cast<DataType&>(this->data_);

  bool bAnswer = true;

  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int nSize = this->getSize();
  if (rComm.isMaster() == true) {
    ofs.write(reinterpret_cast<const char*>(&nSize), sizeof(int));
  }
  const int offset = sizeof(int);

  const int nNumOfProc = rComm.getNumOfProc();
  const int nRank = rComm.getRank();
  const int nGlobalRows = this->m_nRows;
  const int nGlobalCols = this->m_nCols;
  ;

  if (rComm.isMaster() == true) {
    const int nRows = this->m_nMyRows;
    const int nCols = this->m_nMyCols;
    // const int nBlockSize = this->m_nBlockSize;
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
        const int nGlobalIndex =
            nGlobalRowIndex * nGlobalCols + nGlobalColIndex;
        ofs.seekp(sizeof(double) * nGlobalIndex + offset, std::ios_base::beg);

        const int index = r + nRows * c;  // row-major

        ofs.write(reinterpret_cast<const char*>(&(dataTmp[index])),
                  sizeof(double));
      }
    }
  }

  for (int i = 1; i < nNumOfProc; ++i) {
    if (i == nRank) {
      // slave
      rComm.sendData(this->m_nMyRows);
      rComm.sendData(this->m_nMyCols);
      rComm.sendData(this->m_RowIndexTable);
      rComm.sendData(this->m_ColIndexTable);
      rComm.sendData(this->data_);
    } else if (nRank == 0) {
      // master
      int nRows;
      rComm.receiveData(nRows, i);
      int nCols;
      rComm.receiveData(nCols, i);
      std::vector<int> rowIndexTable;
      rComm.receiveData(rowIndexTable, i);
      std::vector<int> colIndexTable;
      rComm.receiveData(colIndexTable, i);
      std::vector<double> buf;
      rComm.receiveData(buf, i);

      for (int c = 0; c < nCols; ++c) {
        const int nGlobalColIndex = colIndexTable[c];
        if (nGlobalColIndex >= nGlobalCols) {
          continue;
        }

        for (int r = 0; r < nRows; ++r) {
          const int nGlobalRowIndex = rowIndexTable[r];
          if (nGlobalRowIndex >= nGlobalRows) {
            continue;
          }

          const int nGlobalIndex =
              nGlobalRowIndex * nGlobalCols + nGlobalColIndex;
          ofs.seekp(sizeof(double) * nGlobalIndex + offset, std::ios_base::beg);

          const int index = r + nRows * c;  // row-major
          ofs.write(reinterpret_cast<const char*>(&(buf[index])),
                    sizeof(double));
        }
      }
    }

    rComm.barrier();
  }

  // const double dEndTime = rComm.getTime();
  // if (rComm.isMaster() == true){
  //  std::cout << "TlDenseGeneralMatrix_blacs::save() time:" << (dEndTime -
  //  dStartTime)
  //  << std::endl;
  //}

  rComm.barrier();
  assert(rComm.checkNonBlockingCommunications());
  return bAnswer;
}

std::vector<double> TlDistributedVector::getVector() const {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  DataType& dataTmp = const_cast<DataType&>(this->data_);

  std::vector<double> ans(this->getSize(), 0.0);
  const int nRows = this->m_nMyRows;
  const int nCols = this->m_nMyCols;
  const int nGlobalRows = this->m_nRows;
  const int nGlobalCols = this->m_nCols;
  ;
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
      ans[nGlobalIndex] = dataTmp[index];
    }
  }

  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(ans);

  return ans;
}

double operator*(const TlDistributedVector& X, const TlDistributedVector& Y) {
  // this const_cast is requiered for PGI compiler
  // "error: expression must be an lvalue or a function designator"
  TlDistributedVector& Xtmp = const_cast<TlDistributedVector&>(X);
  TlDistributedVector& Ytmp = const_cast<TlDistributedVector&>(Y);

  assert(X.getSize() == Y.getSize());

  const int N = X.getSize();
  double DOT = 0.0;
  const int IX = 1;
  const int JX = 1;
  const int INCX = 1;
  const int IY = 1;
  const int JY = 1;
  const int INCY = 1;

  pddot_(&N, &DOT, &(Xtmp.data_[0]), &IX, &JX, X.m_pDESC, &INCX,
         &(Ytmp.data_[0]), &IY, &JY, Y.m_pDESC, &INCY);

  // INCX=1の場合、列のプロセスにのみ値が返る。
  // したがって、全プロセスに値を送信する必要がある。
  {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.broadcast(DOT);
  }

  return DOT;
}

TlDistributedVector operator*(const TlDistributedVector& X, const double& Y) {
  TlDistributedVector Z(X);
  Z *= Y;

  return Z;
}

TlDistributedVector operator*(const double& Y, const TlDistributedVector& X) {
  return (X * Y);
}
