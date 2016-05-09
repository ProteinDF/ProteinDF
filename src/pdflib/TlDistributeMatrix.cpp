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
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ios>
#include <fstream>
#include <limits>
#include <algorithm>
#include <numeric>

#include "scalapack.h"
#include "TlMath.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlFileMatrix.h"
#include "TlCommunicate.h"

//#define DEBUGOUT_GET_PARTIAL_MATRIX

/// @file
/// OpenMPに関する注意事項
/// #pragma omp critical (TlDistributeMatrix_gpmServerTasks_) は
/// this->gpmServerTasks_ メンバ変数の書き換え処理を行う箇所に設置している。
/// 同様に #pragma omp critical (TlDistributeMatrix_gpmClientTasks_) は
/// this->gpmClientTasks_ メンバ変数の書き換え処理を行う箇所に設置している。


int TlDistributeMatrix::systemBlockSize_ = 64;
const std::size_t TlDistributeMatrix::FILE_BUFFER_SIZE = 100*1024*1024; // 100 MB
const TlDistributeMatrix::size_type TlDistributeMatrix::MAX_LOOP = std::numeric_limits<int>::max();
bool TlDistributeMatrix::isUsingPartialIO = false;

TlScalapackContext* TlScalapackContext::m_pTlScalapackContextInstance = NULL;
int TlScalapackContext::m_nContext = 0;
int TlScalapackContext::m_nProc = 0;
int TlScalapackContext::m_nRank = 0;
int TlScalapackContext::m_nProcGridRow = 0;
int TlScalapackContext::m_nProcGridCol = 0;

TlScalapackContext::TlScalapackContext()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // initialize
    Cblacs_pinfo(&(this->m_nRank), &(this->m_nProc));
    if (this->m_nProc < 1) {
        if (this->m_nRank == 0) {
            this->m_nProc = rComm.getNumOfProc();
        }
        Cblacs_setup(&(this->m_nRank), &(this->m_nProc));
    }
    Cblacs_get(-1, 0, &(this->m_nContext));

    // // process grid
    TlScalapackContext::m_nProc = rComm.getNumOfProc();
    {
        std::vector<int> f = TlMath::factor(TlScalapackContext::m_nProc);
        int r = 1;
        int c = 1;
        int f_size = f.size();
        for (int i = 0; i < f_size; ++i) {
            if ((i % 2) == 0) {
                c *= f[i];
            } else {
                r *= f[i];
            }
        }
        TlScalapackContext::m_nProcGridRow = r;
        TlScalapackContext::m_nProcGridCol = c;
    }

    //   Cblacs_gridinit(&(TlScalapackContext::m_nContext), "Col",
    //          TlScalapackContext::m_nProcGridRow, TlScalapackContext::m_nProcGridCol);
    Cblacs_gridinit(&(TlScalapackContext::m_nContext), "Row-major",
                    TlScalapackContext::m_nProcGridRow, TlScalapackContext::m_nProcGridCol);
}

TlScalapackContext::~TlScalapackContext()
{
    Cblacs_gridexit(TlScalapackContext::m_nContext);
    TlScalapackContext::m_nContext = 0;
    Cblacs_exit(1);
}

void TlScalapackContext::getData(int& rContext, int& rProc, int& rRank,
                                 int& rProcGridRow, int& rProcGridCol)
{
    if (TlScalapackContext::m_pTlScalapackContextInstance == NULL) {
        TlScalapackContext::m_pTlScalapackContextInstance = new TlScalapackContext();
    }
    assert(TlScalapackContext::m_pTlScalapackContextInstance != NULL);

    rContext =  TlScalapackContext::m_nContext;
    rProc = TlScalapackContext::m_nProc;
    rRank = TlScalapackContext::m_nRank;
    rProcGridRow = TlScalapackContext::m_nProcGridRow;
    rProcGridCol = TlScalapackContext::m_nProcGridCol;
}

void TlScalapackContext::finalize()
{
    if (TlScalapackContext::m_pTlScalapackContextInstance != NULL) {
        delete TlScalapackContext::m_pTlScalapackContextInstance;
        TlScalapackContext::m_pTlScalapackContextInstance = NULL;
    }
}

////////////////////////////////////////////////////////////////////////
void TlDistributeMatrix::LocalMatrixHeader::load(std::ifstream* pIs)
{
    pIs->read((char*)&(this->type), sizeof(int));
    pIs->read((char*)&(this->globalRow), sizeof(index_type));
    pIs->read((char*)&(this->globalCol), sizeof(index_type));
    pIs->read((char*)&(this->myRows), sizeof(index_type));
    pIs->read((char*)&(this->myCols), sizeof(index_type));
    pIs->read((char*)&(this->rank), sizeof(int));
    pIs->read((char*)&(this->numOfProcs), sizeof(int));
    pIs->read((char*)&(this->procGridRow), sizeof(int));
    pIs->read((char*)&(this->procGridCol), sizeof(int));
    pIs->read((char*)&(this->myProcRow), sizeof(int));
    pIs->read((char*)&(this->myProcCol), sizeof(int));
    pIs->read((char*)&(this->blockSize), sizeof(int));
}


void TlDistributeMatrix::LocalMatrixHeader::save(std::ofstream* pOs)
{
    pOs->write(reinterpret_cast<const char*>(&(this->type)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->globalRow)), sizeof(index_type));
    pOs->write(reinterpret_cast<const char*>(&(this->globalCol)), sizeof(index_type));
    pOs->write(reinterpret_cast<const char*>(&(this->myRows)), sizeof(index_type));
    pOs->write(reinterpret_cast<const char*>(&(this->myCols)), sizeof(index_type));
    pOs->write(reinterpret_cast<const char*>(&(this->rank)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->numOfProcs)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->procGridRow)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->procGridCol)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->myProcRow)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->myProcCol)), sizeof(int));
    pOs->write(reinterpret_cast<const char*>(&(this->blockSize)), sizeof(int));
}


////////////////////////////////////////////////////////////////////////
TlDistributeMatrix::TlDistributeMatrix(const index_type row, const index_type col)
    : log_(TlLogging::getInstance()),
      m_nContext(0), m_nRows(row), m_nCols(col), m_nBlockSize(TlDistributeMatrix::systemBlockSize_),
      pData_(NULL)
{
    this->initialize();
}


TlDistributeMatrix::TlDistributeMatrix(const TlDistributeMatrix& rhs)
    : log_(TlLogging::getInstance()),
      m_nContext(0), m_nRows(rhs.m_nRows), m_nCols(rhs.m_nCols), m_nBlockSize(rhs.m_nBlockSize),
      pData_(NULL)
{
    this->initialize();
    std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(), this->pData_);
}


TlDistributeMatrix::TlDistributeMatrix(const TlDistributeSymmetricMatrix& rhs)
    : log_(TlLogging::getInstance()),
      m_nContext(0), m_nRows(rhs.m_nRows), m_nCols(rhs.m_nCols), m_nBlockSize(rhs.m_nBlockSize),
      pData_(NULL)
{
    this->initialize();
    std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(), this->pData_);
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());
    const int numOfProcs = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    std::vector<std::vector<index_type> > indexArray(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    // 送信準備
    const index_type numOfLocalRows = rhs.m_nMyRows;
    const index_type numOfLocalCols = rhs.m_nMyCols;
    for (index_type c = 0; c < numOfLocalCols; ++c) { 
        const index_type globalCol = rhs.m_ColIndexTable[c];
        assert(globalCol < this->getNumOfCols());
        for (index_type r = 0; r < numOfLocalRows; ++r) {
            const index_type globalRow = rhs.m_RowIndexTable[r];
            assert(globalCol < this->getNumOfRows());

            if (globalRow >= globalCol) {
                // (globalRow, globalCol)のデータを(globalCol, globalRow)にコピー
                const double value = rhs.pData_[numOfLocalRows * c + r];
                const size_type localIndex = this->getIndex(globalCol, globalRow);
                if (localIndex != -1) {
                    assert((std::size_t)localIndex < this->getNumOfMyElements());
                    this->pData_[localIndex] = value;
                } else {
                    const int proc = this->getProcIdForIndex(globalCol, globalRow);
                    assert(proc != rank);
                    indexArray[proc].push_back(globalCol);
                    indexArray[proc].push_back(globalRow);
                    values[proc].push_back(value);
                }
            }
        }
    }
    
    std::vector<std::size_t> sizeArray(numOfProcs);
    for (int i = 0; i < numOfProcs; ++i) {
        const std::size_t size = indexArray[i].size() / 2;
        assert(size == values[i].size());
        sizeArray[i] = size;
    }

    // 送信
    for (int i = 0; i < numOfProcs; ++i) {
        if (i == rank) {
            continue;
        }

        rComm.iSendData(sizeArray[i], i, TAG_CONSTRUCTOR_SIZE);
        if (sizeArray[i] > 0) {
            rComm.iSendDataX(&(indexArray[i][0]), sizeArray[i] * 2, i, TAG_CONSTRUCTOR_INDEX);
            rComm.iSendDataX(&(values[i][0]), sizeArray[i], i, TAG_CONSTRUCTOR_VALUE);
        }
    }

    // 受信
    enum {
        WAIT_SIZE = 1,
        RECV_SIZE = 2,
        WAIT_INDEX = 4,
        RECV_INDEX = 8,
        WAIT_VALUES = 16,
        RECV_VALUES = 32,
        FINISHED = 64
    };
    bool isFinished = false;
    std::vector<unsigned int> state(numOfProcs, 0);
    std::vector<std::size_t> recvSize(numOfProcs);
    std::vector<index_type*> recvIndex(numOfProcs, (index_type*)NULL);
    std::vector<double*> recvValues(numOfProcs, (double*)NULL);
    while (isFinished == false) {
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rank) {
                continue;
            }

            if ((state[i] & FINISHED) != 0) {
                continue;
            }

            if ((state[i] & WAIT_SIZE) == 0) {
                rComm.iReceiveData(recvSize[i], i, TAG_CONSTRUCTOR_SIZE);
                state[i] |= WAIT_SIZE;
            } 
                
            if (((state[i] & WAIT_SIZE) == WAIT_SIZE) &&
                ((state[i] & RECV_SIZE) == 0) &&
                (rComm.test(recvSize[i]) == true)) {
                rComm.wait(recvSize[i]);
                state[i] |= RECV_SIZE;

                if (recvSize[i] == 0) {
                    state[i] |= WAIT_INDEX;
                    state[i] |= RECV_INDEX;
                    state[i] |= WAIT_VALUES;
                    state[i] |= RECV_VALUES;
                    state[i] |= FINISHED;
                }
                continue;
            }

            if ((state[i] & RECV_SIZE) != 0) {
                if ((state[i] & WAIT_INDEX) == 0) {
                    const std::size_t bufSize = recvSize[i] * 2;
                    recvIndex[i] = new index_type[bufSize];
                    rComm.iReceiveDataX(recvIndex[i], bufSize, i, TAG_CONSTRUCTOR_INDEX);
                    state[i] |= WAIT_INDEX;
                }

                if ((state[i] & WAIT_VALUES) == 0) {
                    const std::size_t bufSize = recvSize[i];
                    recvValues[i] = new double[bufSize];
                    rComm.iReceiveDataX(recvValues[i], bufSize, i, TAG_CONSTRUCTOR_VALUE);
                    state[i] |= WAIT_VALUES;
                }
            }
            
            if (((state[i] & WAIT_INDEX) != 0) &&
                ((state[i] & RECV_INDEX) == 0) &&
                (rComm.test(recvIndex[i]) == true)) {
                rComm.wait(recvIndex[i]);
                state[i] |= RECV_INDEX;
            }

            if (((state[i] & WAIT_VALUES) != 0) &&
                ((state[i] & RECV_VALUES) == 0) &&
                (rComm.test(recvValues[i]) == true)) {
                rComm.wait(recvValues[i]);
                state[i] |= RECV_VALUES;
            }

            if (((state[i] & RECV_INDEX) == RECV_INDEX) &&
                ((state[i] & RECV_VALUES) == RECV_VALUES)) {
                const std::size_t bufSize = recvSize[i];
                for (std::size_t j = 0; j < bufSize; ++j) {
                    const index_type globalRow = recvIndex[i][j * 2];
                    const index_type globalCol = recvIndex[i][j * 2 + 1];
                    const double value = recvValues[i][j];
                    const size_type localIndex = this->getIndex(globalRow, globalCol);
                    assert(localIndex >= 0);
                    this->pData_[localIndex] = value;
                }
                state[i] |= FINISHED;
            }
        }

        // check if job is finished
        int numOfFinishedProc = 0;
        for (int i = 0; i < numOfProcs; ++i) {
            if ((state[i] & FINISHED) != 0) {
                ++numOfFinishedProc;
            }
        }
        if (numOfFinishedProc == numOfProcs -1) {
            isFinished = true;
        }
    }

    // バッファの後始末
    for (int i = 0; i < numOfProcs; ++i) {
        if (recvIndex[i] != NULL) {
            delete[] recvIndex[i];
            recvIndex[i] = NULL;
        }
        if (recvValues[i] != NULL) {
            delete[] recvValues[i];
            recvValues[i] = NULL;
        }
    }

    // 送信wait
    for (int i = 0; i < numOfProcs; ++i) {
        if (i == rank) {
            continue;
        }

        rComm.wait(sizeArray[i]);
        if (sizeArray[i] > 0) {
            rComm.wait(&(indexArray[i][0]));
            rComm.wait(&(values[i][0]));
        }
    }

    rComm.barrier(true);
    assert(rComm.checkNonBlockingCommunications());
}


TlDistributeMatrix::TlDistributeMatrix(const TlDistributeVector& rhs,
                                       const index_type row, const index_type col)
    : log_(TlLogging::getInstance()),
      m_nContext(0), m_nRows(row), m_nCols(col), m_nBlockSize(TlDistributeMatrix::systemBlockSize_),
      pData_(NULL)
{
    this->initialize();

    TlDistributeVector rhs_tmp = const_cast<TlDistributeVector&>(rhs); // for PGI compiler
    const std::size_t bufSize = this->getNumOfMyElements();
    std::copy(&(rhs_tmp.data_[0]), &(rhs_tmp.data_[0]) + bufSize, this->pData_);
}


TlDistributeMatrix::TlDistributeMatrix(const TlRowVectorMatrix& rhs)
    : log_(TlLogging::getInstance()),
      m_nContext(0), m_nRows(rhs.getNumOfRows()), m_nCols(rhs.getNumOfCols()),
      m_nBlockSize(TlDistributeMatrix::systemBlockSize_),
      pData_(NULL)
{
    this->initialize();

    TlCommunicate& rComm = TlCommunicate::getInstance();
    // const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlVector rowVec(numOfCols);
    for (index_type row = 0; row < numOfRows; ++row) {
        const int charge = rhs.getSubunitID(row);
        if (charge == myRank) {
            rowVec = rhs.getVector(row);
            assert(rowVec.getSize() == numOfCols);
        }
        rComm.broadcast(&(rowVec[0]), numOfCols, charge);
        
        for (index_type col = 0; col < numOfCols; ++col) {
            this->set(row, col, rowVec[col]);
        }
    }
}


TlDistributeMatrix::~TlDistributeMatrix()
{
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
}


void TlDistributeMatrix::setSystemBlockSize(int blockSize)
{
    TlDistributeMatrix::systemBlockSize_ = blockSize;
}


void TlDistributeMatrix::setUsingPartialIO(bool isUsePIO)
{
    TlDistributeMatrix::isUsingPartialIO = isUsePIO;
}


std::size_t TlDistributeMatrix::getMemSize() const
{
//     std::size_t answer = sizeof(TlDistributeMatrix);
//     answer += (sizeof(int) * 2 + sizeof(double)) * this->m_nMyRows * this->m_nMyCols;
//     return answer;
    const int blockSize = this->getBlockSize();
    const int rowCycle = (this->getNumOfRows() / blockSize) / this->m_nProcGridRow +1;
    const int colCycle = (this->getNumOfCols() / blockSize) / this->m_nProcGridCol +1;
    const std::size_t numOfBlocks = blockSize * blockSize * rowCycle * colCycle;
    const std::size_t answer = numOfBlocks * sizeof(double);

    return answer;
}


void TlDistributeMatrix::resize(const index_type row, const index_type col)
{
    assert(row > 0);
    assert(col > 0);
    if ((row == this->getNumOfRows()) &&
        (col == this->getNumOfCols())) {
        // do not need operation.
        return;
    }
    
    // backup old condition
    TlDistributeMatrix tmp(*this);

    // new size and zero clear
    this->m_nRows = row;
    this->m_nCols = col;
    this->initialize();

    // copy
    const int copyMaxRows = std::min(row, tmp.getNumOfRows());
    const int copyMaxCols = std::min(col, tmp.getNumOfCols());
    if (this->getBlockSize() == tmp.getBlockSize()) {
        // ブロックサイズが同じでプロセスグリッドが同一なら
        // ローカルコピーで十分。
        std::vector<int>::const_iterator rEnd = tmp.m_RowIndexTable.end();
        std::vector<int>::const_iterator cEnd = tmp.m_ColIndexTable.end();
        int index = 0;
        for (std::vector<int>::const_iterator c = tmp.m_ColIndexTable.begin(); c != cEnd; ++c) {
            const int col = *c;
            
            for (std::vector<int>::const_iterator r = tmp.m_RowIndexTable.begin(); r != rEnd; ++r) {
                const int row = *r;
                
                if ((row < copyMaxRows) && (col < copyMaxCols)) {
                    const double value = tmp.pData_[index];
                    const int l = this->getIndex(row, col);
                    if (l != -1) {
                        this->pData_[l] = value;
                    }
                }
                ++index;
            }
        }
    } else {
        // 全通信を行う: too slow
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int nProc = rComm.getNumOfProc();
        const int nRank = rComm.getRank();
        for (int i = 0; i < nProc; ++i) {
            index_type numOfLocalRows;
            index_type numOfLocalCols;
            std::vector<index_type> rowIndexTable;
            std::vector<index_type> colIndexTable;
            double* pBuf = NULL;
            
            if (i == nRank) {
                numOfLocalRows = this->m_nMyRows;
                numOfLocalCols = this->m_nMyCols;
                rowIndexTable = tmp.m_RowIndexTable;
                colIndexTable = tmp.m_ColIndexTable;
                const std::size_t bufSize = numOfLocalRows * numOfLocalCols;
                pBuf = new double[bufSize];
                std::copy(this->pData_, this->pData_ + bufSize, pBuf);
            }
            rComm.broadcast(numOfLocalRows, i);
            rComm.broadcast(numOfLocalCols, i);
            rComm.broadcast(rowIndexTable, i);
            rComm.broadcast(colIndexTable, i);
            const std::size_t bufSize = numOfLocalRows * numOfLocalCols;
            if (i != nRank) {
                pBuf = new double[bufSize];
            }
            rComm.broadcast(pBuf, bufSize, i);
            
            std::vector<int>::const_iterator rEnd = rowIndexTable.end();
            std::vector<int>::const_iterator cEnd = colIndexTable.end();
            int index = 0;
            for (std::vector<int>::const_iterator c = colIndexTable.begin(); c != cEnd; ++c) {
                const int col = *c;
                
                for (std::vector<int>::const_iterator r = rowIndexTable.begin(); r != rEnd; ++r) {
                    const int row = *r;
                    
                    if ((row < copyMaxRows) && (col < copyMaxCols)) {
                        const double value = pBuf[index];
                        const int l = this->getIndex(row, col);
                        if (l != -1) {
                            this->pData_[l] = value;
                        }
                    }
                    ++index;
                }
            }

            if (pBuf != NULL) {
                delete[] pBuf;
                pBuf = NULL;
            }
            rComm.barrier();
        }
    }
}


TlDistributeVector TlDistributeMatrix::getVector() const
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //DataType& data_tmp = const_cast<DataType&>(this->data_);

    std::vector<double> answer(this->getNumOfMyElements());
    std::copy(this->pData_, this->pData_ + this->getNumOfMyElements(), answer.begin());

    return TlDistributeVector(answer, this->getNumOfRows() * this->getNumOfCols());
}


TlVector TlDistributeMatrix::getRowVector(const index_type nRow) const
{
    assert((0 <= nRow) && (nRow < this->m_nRows));

    const int nGlobalCols = this->m_nCols;
    TlVector answer(nGlobalCols);

    if (std::binary_search(this->m_RowIndexTable.begin(), this->m_RowIndexTable.end(), nRow) == true) {
        std::vector<int>::const_iterator pEnd = this->m_ColIndexTable.end();
        for (std::vector<int>::const_iterator p = this->m_ColIndexTable.begin(); p != pEnd; ++p) {
            const int nCol = *p;
            if (nCol < nGlobalCols) {
                const int local_index = this->getIndex(nRow, nCol);
                assert(local_index != -1);
                answer[nCol] = this->pData_[local_index];
            }
        }
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);

    return answer;
}


TlVector TlDistributeMatrix::getColumnVector(const index_type nCol) const
{
    assert((0 <= nCol) && (nCol < this->m_nCols));

    const int nGlobalRows = this->m_nRows;
    TlVector answer(nGlobalRows);

    if (std::binary_search(this->m_ColIndexTable.begin(), this->m_ColIndexTable.end(), nCol) == true) {
        std::vector<int>::const_iterator pEnd = this->m_RowIndexTable.end();
        for (std::vector<int>::const_iterator p = this->m_RowIndexTable.begin(); p != pEnd; ++p) {
            const int nRow = *p;
            if (nRow < nGlobalRows) {
                const int local_index = this->getIndex(nRow, nCol);
                assert(local_index != -1);
                answer[nRow] = this->pData_[local_index];
            }
        }
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);

    return answer;
}


double TlDistributeMatrix::getMaxAbsoluteElement(index_type* pOutRow, index_type* pOutCol) const
{
    int nRow = 0;
    int nCol = 0;
    double dMaxValue = this->getLocalMaxAbsoluteElement(&nRow, &nCol);

    // communication
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int nProc = rComm.getNumOfProc();
        const int nRank = rComm.getRank();
        for (int i = 1; i < nProc; ++i) {
            if (i == nRank) {
                rComm.sendData(dMaxValue);
                rComm.sendData(nRow);
                rComm.sendData(nCol);
            } else if (rComm.isMaster() == true) {
                double receivedMaxValue = 0.0;
                int receivedRow = 0;
                int receivedCol = 0;
                rComm.receiveData(receivedMaxValue, i);
                rComm.receiveData(receivedRow, i);
                rComm.receiveData(receivedCol, i);
                if (dMaxValue < receivedMaxValue) {
                    dMaxValue = receivedMaxValue;
                    nRow = receivedRow;
                    nCol = receivedCol;
                }
            }
            //rComm.barrier();
        }

        rComm.broadcast(dMaxValue);
        rComm.broadcast(nRow);
        rComm.broadcast(nCol);
    }

    if (pOutRow != NULL) {
        *pOutRow = nRow;
    }
    if (pOutCol != NULL) {
        *pOutCol = nCol;
    }

    return dMaxValue;
}


double TlDistributeMatrix::getLocalMaxAbsoluteElement(index_type* pOutRow, index_type* pOutCol) const
{
    double dMaxValue = 0.0;
    std::size_t index = 0;
    const size_type nMaxIndex = this->getNumOfMyElements();
    for (size_type i = 0; i < nMaxIndex; ++i) {
        const double v = std::fabs(this->pData_[i]);
        if (dMaxValue < v) {
            dMaxValue = v;
            index = i;
        }
    }

    if ((pOutRow != NULL) || (pOutCol != NULL)) {
        ldiv_t t = std::ldiv(index, this->m_nMyRows);
        index_type myRow = t.rem;
        index_type myCol = t.quot;
        index_type nRow = 0;
        index_type nCol = 0;
        if ((myRow < static_cast<index_type>(this->m_RowIndexTable.size())) &&
            (myCol < static_cast<index_type>(this->m_ColIndexTable.size()))) {
            nRow = this->m_RowIndexTable[myRow];
            nCol = this->m_ColIndexTable[myCol];
        } else {
            dMaxValue = 0.0;
        }

        if (pOutRow != NULL) {
            *pOutRow = nRow;
        }
        if (pOutCol != NULL) {
            *pOutCol = nCol;
        }
    }

    return dMaxValue;
}


double TlDistributeMatrix::getRMS() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double sum2 = 0.0;
    const size_type maxIndex = this->getNumOfMyElements();
    for (size_type i = 0; i < maxIndex; ++i) {
        double tmp = this->pData_[i];
        sum2 += tmp * tmp;
    }
    rComm.allReduce_SUM(sum2);

    const double elements = this->getNumOfRows() * this->getNumOfCols();
    
    const double rms = std::sqrt(sum2 / elements);
    return rms;
}


// TlSparseMatrix TlDistributeMatrix::getPartialMatrix(const double threshold) const
// {
//     const int globalRows = this->m_nRows;
//     const int globalCols = this->m_nCols;

//     // 自分のデータのうち、有意なデータを抽出する
//     TlSparseMatrix m(globalRows, globalCols);
//     {
//         const index_type numOfLocalRows = this->m_nMyRows;
//         const index_type numOfLocalCols = this->m_nMyCols;
//         for (int c = 0; c < numOfLocalCols; ++c) {
//             const int globalColIndex = this->m_ColIndexTable[c];
//             const int index_tmp = numOfLocalRows * c;

//             for (int r = 0; r < numOfLocalRows; ++r) {
//                 const int index = r + index_tmp;
//                 const double value = this->pData_[index];
//                 if (std::fabs(value) > threshold) {
//                     const int globalRowIndex = this->m_ColIndexTable[r];
//                     m.set(globalRowIndex, globalColIndex, value);
//                 }
//             }
//         }
//     }

//     // 集計
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int proc = rComm.getNumOfProc();
//     const int rank = rComm.getRank();

//     int numOfElements = m.getSize();
//     rComm.allReduce_SUM(numOfElements);
//     const int averageNumOfElements = (numOfElements + proc -1) / proc;

//     // 平均よりも多い要素を持っているSparseMatrixは他に渡す
//     TlSparseMatrix diffM(globalRows, globalCols);
//     const int maxCycle = numOfElements - averageNumOfElements;
//     for (int i = 0; i < maxCycle; ++i) {
//         int row = 0;
//         int col = 0;
//         const double value = m.pop(&row, &col);
//         diffM.set(row, col, value);
//     }

//     // diffMの集計
//     int numOfDiffMElements = diffM.getSize();
//     rComm.allReduce_SUM(numOfDiffMElements);
//     const int numOfHolds = (numOfDiffMElements + proc -1) / proc; // この数だけ担当する

//     // diffMの平均化
//     const int myStartCount = numOfHolds * rank;
//     const int myEndCount = numOfHolds * (rank +1);
//     int count = 0;
//     TlSparseMatrix addM(globalRows, globalCols);
//     for (int i = 0; i < proc; ++i) {
//         TlSparseMatrix tmp(globalRows, globalCols);
//         if (i == rank) {
//             tmp = diffM;
//         }
//         rComm.broadcast(diffM, i);

//         TlSparseMatrix::const_iterator pEnd = tmp.end();
//         for (TlSparseMatrix::const_iterator p = tmp.begin(); p != pEnd; ++p) {
//             if ((myStartCount <= count) && (count < myEndCount)) {
//                 addM.set(*p);
//             }
//             ++count;
//         }
//     }

//     // 最終成果物
//     m.merge(addM);

//     return m;
// }


// void TlDistributeMatrix::getPartialMatrix(TlSparseMatrix& ioMatrix) const
// {
//     assert(this->getNumOfRows() == ioMatrix.getNumOfRows());
//     assert(this->getNumOfCols() == ioMatrix.getNumOfCols());

//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const int nProc = rComm.getNumOfProc();
//     const int nRank = rComm.getRank();
//     for (int i = 0; i < nProc; ++i) {
//         TlSparseMatrix tmp;
//         if (i == nRank) {
//             tmp = ioMatrix;
//         }
//         rComm.broadcast(tmp, i);

//         TlSparseMatrix::iterator pEnd = tmp.end();
//         for (TlSparseMatrix::iterator p = tmp.begin(); p != pEnd; ++p) {
//             const unsigned long index = p->first;
//             int nRow = 0;
//             int nCol = 0;
//             tmp.index(index, &nRow, &nCol);
//             const double dValue = this->get(nRow, nCol);
//             p->second = dValue;
//         }

//         if (i == nRank) {
//             ioMatrix = tmp;
//         }
//     }
// }


// 複数スレッドから同時に呼び出さない
bool TlDistributeMatrix::getSparseMatrixX(TlSparseMatrix* pMatrix, bool isFinalize) const
{
    bool answer = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (pMatrix != NULL) {
        // Client task通信処理
#pragma omp critical (TlDistributeMatrix__getSparseMatrixX)
        {
            this->getSparseMatrixX_registerTask(pMatrix);
            answer = this->getPartialMatrix_ClientTasks(pMatrix);
        }
    } else {
        // Server処理
        this->getPartialMatrix_ServerTasks(isFinalize);
    }

    if (isFinalize == true) {
        rComm.checkNonBlockingCommunications();
        assert(this->gpmClientTasks_.size() == 0);

#pragma omp critical (TlDistributeMatrix_gpmClientTasks_)
        {
            this->gpmClientTasks_.clear();
        }
    }
    
    return answer;
}


// 複数スレッドから同時に呼び出さない
bool TlDistributeMatrix::getPartialMatrixX(TlPartialMatrix* pMatrix, bool isFinalize) const
{
    bool answer = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

#pragma omp critical (TlDistributeMatrix__getPartialMatrixX)
    {
        // Client task通信処理
        this->getPartialMatrixX_registerTask(pMatrix);

        answer = this->getPartialMatrix_ClientTasks(pMatrix);

        // Server処理
        this->getPartialMatrix_ServerTasks(isFinalize);
    }

    if (isFinalize == true) {
        rComm.checkNonBlockingCommunications();
        assert(this->gpmClientTasks_.size() == 0);
#pragma omp critical (TlDistributeMatrix_gpmClientTasks_)
        {
            this->gpmClientTasks_.clear();
        }
    }

    return answer;
}


void TlDistributeMatrix::getPartialMatrix_ServerTasks(const bool isFinalize) const
{
#pragma omp critical (TlDistributeMatrix__getPartialMatrix_ServerTasks)
    {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    //const int numOfProcs = rComm.getNumOfProc();

    int requestProc = 0;
    if (this->isWaitingRequestHandShake_ == false) {
        // セッション開始
        rComm.iReceiveDataFromAnySource(this->sessionId_, TAG_GPM_SESSION_ID);
        this->isWaitingRequestHandShake_ = true;

        if (this->isDebugOut_GPM_ == true) {
            std::cerr << TlUtils::format("[%d] SRV [**]: recv session id from any.",
                                         rComm.getRank())
                      << std::endl;
        }
    }

    if ((this->isWaitingRequestHandShake_ == true) &&
        (rComm.test(this->sessionId_, &requestProc) == true)) {
        // セッション番号受信
        rComm.wait(this->sessionId_);
        this->isWaitingRequestHandShake_ = false;

        if (this->isDebugOut_GPM_ == true) {
            std::cerr << TlUtils::format("[%d] SRV [OK]: recv session id from [%d].",
                                         rComm.getRank(),
                                         requestProc)
                      << std::endl;
        }

        // taskの登録
        GPM_ServerTask task;
        task.state = 0;
        task.sessionId = this->sessionId_;
        task.requestProc = requestProc;
        task.numOfComponents = 0;

// #pragma omp critical (TlDistributeMatrix_gpmServerTasks_)
//         {
            this->gpmServerTasks_.push_back(task);
//         }
    }

    GpmServerTasksType::iterator itEnd = this->gpmServerTasks_.end();
    for (GpmServerTasksType::iterator it = this->gpmServerTasks_.begin(); it != itEnd; ++it) {
        GPM_ServerTask& task = *it;
        const int proc = task.requestProc;
        if ((task.state & GPM_SERVER_RECV_NUM_OF_COMPONENTS) == 0) {
            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] SRV [**]: recv num of components from [%d]",
                                             rComm.getRank(), proc)
                          << std::endl;
            }

            rComm.iReceiveData(task.numOfComponents,
                               proc,
                               TAG_GPM_NUM_OF_COMPONENTS + task.sessionId);
            task.state |= GPM_SERVER_RECV_NUM_OF_COMPONENTS;
            continue;
        }

        if (((task.state & GPM_SERVER_RECV_NUM_OF_COMPONENTS) != 0) &&
            ((task.state & GPM_SERVER_WAIT_NUM_OF_COMPONENTS) == 0)) {
            if (rComm.test(task.numOfComponents) == true) {
                rComm.wait(task.numOfComponents);
                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] SRV [OK]: recv num of components from [%d]. num = %d",
                                                 rComm.getRank(), proc, task.numOfComponents)
                              << std::endl;
                }
                assert(task.numOfComponents > 0);
                task.state |= GPM_SERVER_WAIT_NUM_OF_COMPONENTS;
                continue;
            }
        }
        
        
        if (((task.state & GPM_SERVER_WAIT_NUM_OF_COMPONENTS) != 0) &&
            ((task.state & GPM_SERVER_RECV_COMPONENTS) == 0)) {
            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] SRV [**]: recv components from [%d]",
                                             rComm.getRank(), proc)
                          << std::endl;
            }

            assert(task.numOfComponents > 0);
            task.components.resize(task.numOfComponents);
            rComm.iReceiveDataX((int*)&(task.components[0]),
                                task.numOfComponents,
                                proc,
                                TAG_GPM_COMPONENTS + task.sessionId);
            task.state |= GPM_SERVER_RECV_COMPONENTS;
            continue;
        }

        if (((task.state & GPM_SERVER_RECV_COMPONENTS) != 0) &&
            ((task.state & GPM_SERVER_WAIT_COMPONENTS) == 0)) {
            if (rComm.test((int*)&(task.components[0])) == true) {
                rComm.wait((int*)&(task.components[0]));
                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] SRV [OK]: recv components from [%d]",
                                                 rComm.getRank(), proc)
                              << std::endl;
                }

                task.state |= GPM_SERVER_WAIT_COMPONENTS;
                continue;
            }
        }

        if (((task.state & GPM_SERVER_WAIT_COMPONENTS) != 0) &&
            ((task.state & GPM_SERVER_SEND_ELEMENT_VALUES) == 0)) {
            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] SRV [**]: send elements to [%d].",
                                             rComm.getRank(), proc)
                          << std::endl;
            }
            
            const int size = task.numOfComponents / 2;
            task.elementValues.resize(size);
            for (int i = 0; i < size; ++i) {
                const index_type row = task.components[i*2];
                const index_type col = task.components[i*2 +1];
                const int localIndex = this->getIndex(row, col);
                assert(localIndex != -1);
                const double value = this->pData_[localIndex];
                task.elementValues[i] = value;
            }
            rComm.iSendDataX((double*)&(task.elementValues[0]),
                             size,
                             proc,
                             TAG_GPM_ELEMENT_VALUES + task.sessionId);
            task.state |= GPM_SERVER_SEND_ELEMENT_VALUES;
            continue;
        }

        if (((task.state & GPM_SERVER_SEND_ELEMENT_VALUES) != 0) &&
            ((task.state & GPM_SERVER_WAIT_ELEMENT_VALUES) == 0)) {
            if (rComm.test((double*)&(task.elementValues[0])) == true) {
                rComm.wait((double*)&(task.elementValues[0]));
                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] SRV [OK]: send elements to [%d].",
                                                 rComm.getRank(), proc)
                              << std::endl;
                }
                task.state |= GPM_SERVER_WAIT_ELEMENT_VALUES;
                continue;
            }
        }
    }

    // 終了したtaskをリストから除く
    //GpmServerTasksType::iterator itEnd = this->gpmServerTasks_.end();
    for (GpmServerTasksType::iterator it = this->gpmServerTasks_.begin(); it != itEnd; ) {
        if ((it->state & GPM_SERVER_WAIT_ELEMENT_VALUES) != 0) {
// #pragma omp critical (TlDistributeMatrix_gpmServerTasks_)
//             {
                it = this->gpmServerTasks_.erase(it);
//             }
            continue;
        }
        ++it;
    }
    
    
    // finalize 処理
    if (isFinalize == true) {
        if (this->isWaitingRequestHandShake_ == true) {
            rComm.cancel(this->sessionId_);
        }
    }
    }
}


void TlDistributeMatrix::getSparseMatrixX_registerTask(TlSparseMatrix* pMatrix) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    
    // Client task登録処理
    if (pMatrix != NULL) {
        // 受信処理
        assert(this->getNumOfRows() == pMatrix->getNumOfRows());
        assert(this->getNumOfCols() == pMatrix->getNumOfCols());
        if (this->gpmClientTasks_.find(pMatrix) == this->gpmClientTasks_.end()) {
            // タスクを初期化
            GPM_ClientTask task;
            task.state.resize(numOfProcs, 0);
            task.sessionIds.resize(numOfProcs, 0);
            for (int proc = 0; proc < numOfProcs; ++proc) {
                task.sessionIds[proc] = this->getPartialMatrix_getSessionId(proc);
            }
            task.components.resize(numOfProcs, std::vector<index_type>());
            task.numOfComponents.resize(numOfProcs);
            task.elementValues.resize(numOfProcs);
            task.isFinished = false;
            
            // 要求された要素をどのプロセスが持っているかをチェック
            for (TlSparseMatrix::iterator it = pMatrix->begin(); it != pMatrix->end(); ++it) {
                index_type row = it->first.row;
                index_type col = it->first.col;
                
                const int localIndex = this->getIndex(row, col);
                if (localIndex != -1) {
                    // 自分が持っていれば格納すべき行列に値を代入する
                    pMatrix->set(row, col, this->getLocal(row, col));
                } else {
                    // 他プロセスが持っていれば、そのプロセスに要求するためのリストを作成する
                    const int proc = this->getProcIdForIndex(row, col);
//                     std::cerr << TlUtils::format("[%d] check elements([%d] has (%d, %d)",
//                                                   rComm.getRank(), proc, row, col)
//                                << std::endl;
                    task.components[proc].push_back(row);
                    task.components[proc].push_back(col);
                }
            }
#pragma omp critical (TlDistributeMatrix_gpmClientTasks_)
            {
                this->gpmClientTasks_[pMatrix] = task;
            }
        }
    }
}


void TlDistributeMatrix::getPartialMatrixX_registerTask(TlPartialMatrix* pMatrix) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();

    // Client task登録処理
    if (pMatrix != NULL) {
        // 受信処理
        assert(this->getNumOfRows() == pMatrix->getNumOfRows());
        assert(this->getNumOfCols() == pMatrix->getNumOfCols());
        if (this->gpmClientTasks_.find(pMatrix) == this->gpmClientTasks_.end()) {
            // タスクを初期化
            GPM_ClientTask task;
            task.state.resize(numOfProcs, 0);
            task.sessionIds.resize(numOfProcs, 0);
            for (int proc = 0; proc < numOfProcs; ++proc) {
                task.sessionIds[proc] = this->getPartialMatrix_getSessionId(proc);
            }
            task.components.resize(numOfProcs, std::vector<index_type>());
            task.numOfComponents.resize(numOfProcs);
            task.elementValues.resize(numOfProcs);
            task.isFinished = false;
            
            // 要求された要素をどのプロセスが持っているかをチェック
            const index_type startRow = pMatrix->getStartRow();
            const index_type startCol = pMatrix->getStartCol();
            const index_type endRow = startRow + pMatrix->getRowRange();
            const index_type endCol = startCol + pMatrix->getColRange();
            for (index_type row = startRow; row < endRow; ++row) {
                for (index_type col = startCol; col < endCol; ++col) {
                    
                    const int localIndex = this->getIndex(row, col);
                    if (localIndex != -1) {
                        // 自分が持っていれば格納すべき行列に値を代入する
                        pMatrix->set(row, col, this->getLocal(row, col));
                    } else {
                        // 他プロセスが持っていれば、そのプロセスに要求するためのリストを作成する
                        const int proc = this->getProcIdForIndex(row, col);
                        task.components[proc].push_back(row);
                        task.components[proc].push_back(col);
                    }
                }
            }
            
#pragma omp critical (TlDistributeMatrix_gpmClientTasks_)
            {
                this->gpmClientTasks_[pMatrix] = task;
            }
        }
    }
}


int TlDistributeMatrix::getPartialMatrix_getSessionId(const int proc) const
{
    int answer = -1;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.getRank() != proc) {
        for (int i = 0; i < MAX_SESSION_ID; ++i) {
            if (this->sessionTable_[proc].test(i) == false) {
                answer = i;
                break;
            }
        }
        assert(answer != -1);
    } else {
        answer = 0;
    }
    
    return answer;
}


bool TlDistributeMatrix::getPartialMatrix_ClientTasks(TlMatrixObject* pMatrix) const
{
    bool answer = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
   
    GpmClientTasks::iterator itEnd = this->gpmClientTasks_.end();
    for (GpmClientTasks::iterator it = this->gpmClientTasks_.begin(); it != itEnd; ++it) {
        TlMatrixObject* pMatrix = it->first;
        GPM_ClientTask& task = it->second;

        assert(task.sessionIds.size() == std::size_t(numOfProcs));

        for (int proc = 0; proc < numOfProcs; ++proc) {
            if ((task.state[proc] & GPM_CLIENT_COUNT_COMPONENTS) == 0) {
                // 要素番号情報の総数を数える
                const std::size_t size = task.components[proc].size();
                task.numOfComponents[proc] = size;
                task.state[proc] |= GPM_CLIENT_COUNT_COMPONENTS;

                if (size == 0) {
                    task.state[proc] = (GPM_CLIENT_COUNT_COMPONENTS |
                                        GPM_CLIENT_SEND_SESSION_ID |
                                        GPM_CLIENT_WAIT_SESSION_ID |
                                        GPM_CLIENT_SEND_NUM_OF_COMPONENTS |
                                        GPM_CLIENT_WAIT_NUM_OF_COMPONENTS |
                                        GPM_CLIENT_SEND_COMPONENTS |
                                        GPM_CLIENT_WAIT_COMPONENTS |
                                        GPM_CLIENT_RECV_ELEMENT_VALUES |
                                        GPM_CLIENT_WAIT_ELEMENT_VALUES);
                }
                continue;
            }
            
            if ((task.state[proc] & GPM_CLIENT_SEND_SESSION_ID) == 0) {
                if (this->trafficControl_[proc] == NULL) {
                    // 相手にセッション番号を送信する
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: send session id to [%d]. num = %d",
                                                     rComm.getRank(), proc, task.numOfComponents[proc])
                                  << std::endl;
                    }
                    
                    rComm.iSendData(task.sessionIds[proc], proc,
                                    TAG_GPM_SESSION_ID);
                    task.state[proc] |= GPM_CLIENT_SEND_SESSION_ID;
                    this->trafficControl_[proc] = pMatrix;
                    continue;
                }
            }

            if (((task.state[proc] & GPM_CLIENT_SEND_SESSION_ID) != 0) &&
                ((task.state[proc] & GPM_CLIENT_WAIT_SESSION_ID) == 0)) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (rComm.test(task.sessionIds[proc]) == true) {
                        rComm.wait(task.sessionIds[proc]);
#ifndef NDEBUG
                        if (this->isDebugOut_GPM_ == true) {
                            std::cerr << TlUtils::format("[%d] CLI [OK]: send session id to [%d].",
                                                         rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif // NDEBUG
                        
                        task.state[proc] |= GPM_CLIENT_WAIT_SESSION_ID;
                        continue;
                    }
                }
            }

            if ((task.state[proc] & GPM_CLIENT_SEND_NUM_OF_COMPONENTS) == 0) {
                if (this->trafficControl_[proc] == pMatrix) {
                    // 相手に要求要素組の総数を送信する
                    assert(task.numOfComponents[proc] != 0);
#ifndef NDEBUG
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: send num of components to [%d]. num = %d",
                                                     rComm.getRank(), proc, task.numOfComponents[proc])
                                  << std::endl;
                    }
#endif // NDEBUG

                    rComm.iSendData(task.numOfComponents[proc], proc,
                                    TAG_GPM_NUM_OF_COMPONENTS + task.sessionIds[proc]);
                    task.state[proc] |= GPM_CLIENT_SEND_NUM_OF_COMPONENTS;
                    this->trafficControl_[proc] = pMatrix;
                    continue;
                }
            }

            if (((task.state[proc] & GPM_CLIENT_SEND_NUM_OF_COMPONENTS) != 0) &&
                ((task.state[proc] & GPM_CLIENT_WAIT_NUM_OF_COMPONENTS) == 0)) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (rComm.test(task.numOfComponents[proc]) == true) {
                        rComm.wait(task.numOfComponents[proc]);
#ifndef NDEBUG
                        if (this->isDebugOut_GPM_ == true) {
                            std::cerr << TlUtils::format("[%d] CLI [OK]: send num of components to [%d].",
                                                         rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif // NDEBUG
                        
                        task.state[proc] |= GPM_CLIENT_WAIT_NUM_OF_COMPONENTS;
                        continue;
                    }
                }
            }
        
            if ((task.state[proc] & GPM_CLIENT_SEND_COMPONENTS) == 0) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: send components to [%d].",
                                                     rComm.getRank(), proc)
                                  << std::endl;
                    }

                    // 相手に要求要素組を送信する
                    rComm.iSendDataX((int*)&(task.components[proc][0]),
                                     task.numOfComponents[proc],
                                     proc,
                                     TAG_GPM_COMPONENTS + task.sessionIds[proc]);
                    task.state[proc] |= GPM_CLIENT_SEND_COMPONENTS;
                    continue;
                }
            }

            if (((task.state[proc] & GPM_CLIENT_SEND_COMPONENTS) != 0) &&
                ((task.state[proc] & GPM_CLIENT_WAIT_COMPONENTS) == 0)) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (rComm.test((int*)&(task.components[proc][0])) == true) {
                        rComm.wait((int*)&(task.components[proc][0]));
#ifndef NDEBUG
                        if (this->isDebugOut_GPM_ == true) {
                            std::cerr << TlUtils::format("[%d] CLI [OK]: send components to [%d].",
                                                         rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif // NDEBUG

                        task.state[proc] |= GPM_CLIENT_WAIT_COMPONENTS;
                        continue;
                    }
                }
            }
            
            if ((task.state[proc] & GPM_CLIENT_RECV_ELEMENT_VALUES) == 0) {
                if (this->trafficControl_[proc] == pMatrix) {
#ifndef NDEBUG
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: recv elements from [%d].",
                                                     rComm.getRank(), proc)
                                  << std::endl;
                    }
#endif // NDEBUG

                    // 行列データを待ち受ける
                    const std::size_t size = task.numOfComponents[proc] / 2;
                    task.elementValues[proc].resize(size);
                    rComm.iReceiveDataX((double*)&(task.elementValues[proc][0]),
                                        size,
                                        proc,
                                        TAG_GPM_ELEMENT_VALUES + task.sessionIds[proc]);
                    task.state[proc] |= GPM_CLIENT_RECV_ELEMENT_VALUES;
                    continue;
                }
            }

            if (((task.state[proc] & GPM_CLIENT_RECV_ELEMENT_VALUES) != 0) &&
                ((task.state[proc] & GPM_CLIENT_WAIT_ELEMENT_VALUES) == 0)) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (rComm.test((double*)&(task.elementValues[proc][0])) == true) {
                        rComm.wait((double*)&(task.elementValues[proc][0]));
                        if (this->isDebugOut_GPM_ == true) {
                            std::cerr << TlUtils::format("[%d] CLI [OK]: recv elements from [%d].",
                                                         rComm.getRank(), proc)
                                      << std::endl;
                        }

                        this->trafficControl_[proc] = NULL;
                        task.state[proc] |= GPM_CLIENT_WAIT_ELEMENT_VALUES;
                        continue;
                    }
                }
            }
        }

        // 受け取りが完了したかどうかチェック
        bool isFinish = true;
        for (int proc = 0; proc < numOfProcs; ++proc) {
            if ((task.state[proc] & GPM_CLIENT_WAIT_ELEMENT_VALUES) == 0) {
                isFinish = false;
                break;
            }
        }

        if (isFinish == true) {
            for (int proc = 0; proc < numOfProcs; ++proc) {
                const std::size_t numOfComponents = task.numOfComponents[proc] / 2;
                for (std::size_t i = 0; i < numOfComponents; ++i) {
                    const index_type row = task.components[proc][i*2   ];
                    const index_type col = task.components[proc][i*2 +1];
                    const double value = task.elementValues[proc][i];

                    pMatrix->set(row, col, value);
//                     if (this->isDebugOut_GPM_ == true) {
//                         std::cerr << TlUtils::format("[%d] CLI [OK]: set (%d, %d) = %f.",
//                                                      rComm.getRank(), row, col, value)
//                                   << std::endl;
//                     }
                }
            }

#ifndef NDEBUG
            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] CLI [##]: task finished.",
                                             rComm.getRank())
                          << std::endl;
            }
#endif // NDEBUG
            
            task.isFinished = true;
        }
    }

    // Client task終了処理
    if (pMatrix != NULL) {
        GpmClientTasks::iterator it = this->gpmClientTasks_.find(pMatrix);
        if (it != this->gpmClientTasks_.end()) {
            if (it->second.isFinished == true) {
                for (int proc = 0; proc < numOfProcs; ++proc) {
                    this->sessionTable_[proc].reset(it->second.sessionIds[proc]);
                }
#pragma omp critical (TlDistributeMatrix_gpmClientTasks_)
                {
                    this->gpmClientTasks_.erase(it);
                }

                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] CLI [##]: task deleted.",
                                                 rComm.getRank())
                              << std::endl;
                }
                
                answer = true;
            }
        }        
    }

    return answer;
}


void TlDistributeMatrix::mergeSparseMatrix(const TlSparseMatrix& M)
{
    assert(M.getNumOfRows() == this->getNumOfRows());
    assert(M.getNumOfCols() == this->getNumOfCols());

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    rComm.checkNonBlockingCommunications();

    // 送信すべきインデックスリストの作成
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    TlSparseMatrix::const_iterator itEnd = M.end();
    for (TlSparseMatrix::const_iterator it = M.begin(); it != itEnd; ++it) {
        index_type globalRow = it->first.row;
        index_type globalCol = it->first.col;
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
    rComm.checkNonBlockingCommunications();
}


void TlDistributeMatrix::mergePartialMatrix(const TlPartialMatrix& M)
{
    assert(M.getNumOfRows() == this->getNumOfRows());
    assert(M.getNumOfCols() == this->getNumOfCols());
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.checkNonBlockingCommunications();
    const int numOfProcs = rComm.getNumOfProc();

    // 送信すべきインデックスリストの作成
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    const index_type startRow = M.getStartRow();
    const index_type startCol = M.getStartCol();
    const index_type endRow = startRow + M.getRowRange();
    const index_type endCol = startCol + M.getColRange();

    for (index_type globalRow = startRow; globalRow < endRow; ++globalRow) {
        for (index_type globalCol = startCol; globalCol < endCol; ++globalCol) {
            const double value = M.get(globalRow, globalCol);
            const size_type index = this->getIndex(globalRow, globalCol);
            if (index != -1) {
                this->pData_[index] += value;
            } else {
                const int targetProc = this->getProcIdForIndex(globalRow, globalCol);
                indexArrays[targetProc].push_back(globalRow);
                indexArrays[targetProc].push_back(globalCol);
                values[targetProc].push_back(value);
            }
        }
    }

    this->mergeMatrix_common(indexArrays, values);
    rComm.checkNonBlockingCommunications();
}


void TlDistributeMatrix::mergeMatrix_common(const std::vector<std::vector<index_type> >& indexArrays,
                                            const std::vector<std::vector<double> >& values)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int numOfProcs = rComm.getNumOfProc();
    const int rank = rComm.getRank();
    
    // インデックスリストを数える
    std::vector<std::size_t> numOfIndeces(numOfProcs);
    for (int proc = 0; proc < numOfProcs; ++proc) {
        numOfIndeces[proc] = indexArrays[proc].size();
    }    

    // インデックスリスト数、インデックスリスト、値を送信する
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc == rank) {
            continue;
        }

        const std::size_t size = numOfIndeces[proc];
        rComm.iSendData(numOfIndeces[proc], proc, TAG_TLDISTRIBUTE_MATRIX__MERGE__NUM_OF_INDECES);
        if (size > 0) {
            assert(values[proc].size() == size / 2);
            rComm.iSendDataX((index_type*)&(indexArrays[proc][0]),
                             size,
                             proc,
                             TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES);
            rComm.iSendDataX((double*)&(values[proc][0]),
                             size / 2,
                             proc,
                             TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES);
        }
    }

    // 行列データを送信する処理 ------------------------------------------------
    // 送られるリストの処理
    std::vector<std::size_t> numOfInputIndeces(numOfProcs, 0);
    std::vector<std::vector<index_type> > inputIndexList(numOfProcs);
    std::vector<std::vector<double> > inputValues(numOfProcs);
    enum {
        RECV_NUM_OF_INDECES = 1,
        WAIT_NUM_OF_INDECES = 2,
        RECV_INDECES = 4,
        WAIT_INDECES = 8,
        RECV_VALUES  = 16,
        WAIT_VALUES  = 32,
        FINISHED = 64
    };
    std::vector<unsigned int> state(numOfProcs, 0);
    state[rank] = (RECV_NUM_OF_INDECES | WAIT_NUM_OF_INDECES |
                   RECV_INDECES | WAIT_INDECES |
                   RECV_VALUES | WAIT_VALUES |
                   FINISHED);

    // input
    int terminateProc = 0;
    while (true) {
        for (int proc = 0; proc < numOfProcs; ++proc) {
            if (rank == proc) {
                continue;
            }
            if ((state[proc] & FINISHED) == FINISHED) {
                // this rank-job has been finished.
                continue;
            }
            
            if ((state[proc] & RECV_NUM_OF_INDECES) != RECV_NUM_OF_INDECES) {
                rComm.iReceiveData(numOfInputIndeces[proc],
                                   proc,
                                   TAG_TLDISTRIBUTE_MATRIX__MERGE__NUM_OF_INDECES);
                state[proc] |= RECV_NUM_OF_INDECES;
            }

            if (((state[proc] & RECV_NUM_OF_INDECES) == RECV_NUM_OF_INDECES) &&
                ((state[proc] & WAIT_NUM_OF_INDECES) != WAIT_NUM_OF_INDECES)) {
                if (rComm.test(numOfInputIndeces[proc]) == true) {
                    rComm.wait(numOfInputIndeces[proc]);
                    state[proc] |= WAIT_NUM_OF_INDECES;
                }
            }

            if (((state[proc] & WAIT_NUM_OF_INDECES) == WAIT_NUM_OF_INDECES) &&
                ((state[proc] & RECV_INDECES) != RECV_INDECES)) {
                const std::size_t size = numOfInputIndeces[proc];
                if (size > 0) {
                    inputIndexList[proc].resize(size);
                    rComm.iReceiveDataX((index_type*)&(inputIndexList[proc][0]),
                                        size,
                                        proc,
                                        TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES);
                    state[proc] |= RECV_INDECES;
                } else {
                    assert(size == 0);
                    state[proc] |= (RECV_INDECES | WAIT_INDECES |
                                    RECV_VALUES  | WAIT_VALUES  |
                                    FINISHED);
                }
            }

            if (((state[proc] & RECV_INDECES) == RECV_INDECES) &&
                ((state[proc] & WAIT_INDECES) != WAIT_INDECES)) {
                if (rComm.test(&(inputIndexList[proc][0])) == true) {
                    rComm.wait(&(inputIndexList[proc][0]));
                    state[proc] |= WAIT_INDECES;
                }
            }
            
            if (((state[proc] & WAIT_NUM_OF_INDECES) == WAIT_NUM_OF_INDECES) &&
                ((state[proc] & RECV_VALUES) != RECV_VALUES)) {
                const std::size_t size = numOfInputIndeces[proc];
                if (size > 0) {
                    inputValues[proc].resize(size / 2);
                    rComm.iReceiveDataX((double*)&(inputValues[proc][0]),
                                        (size / 2),
                                        proc,
                                        TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES);
                    state[proc] |= RECV_VALUES;
                }
            }
            
            if (((state[proc] & RECV_VALUES) == RECV_VALUES) &&
                ((state[proc] & WAIT_VALUES) != WAIT_VALUES)) {
                if (rComm.test(&(inputValues[proc][0])) == true) {
                    rComm.wait(&(inputValues[proc][0]));
                    state[proc] |= WAIT_VALUES;
                }
            }
            
            
            if (((state[proc] & WAIT_VALUES) == WAIT_VALUES) &&
                ((state[proc] & WAIT_INDECES) == WAIT_INDECES)) {
                const std::size_t size = numOfInputIndeces[proc];
                const int max_i = size / 2;
                for (int i = 0; i < max_i; ++i) {
                    const index_type globalRow = inputIndexList[proc][i*2   ];
                    const index_type globalCol = inputIndexList[proc][i*2 +1];
                    const double value = inputValues[proc][i];
                    const int index = this->getIndex(globalRow, globalCol);
                    assert(index != -1);
                    assert((std::size_t)index < this->getNumOfMyElements());
                    this->pData_[index] += value;
                }

                state[proc] |= FINISHED;
                ++terminateProc;
            }
        }

        if (terminateProc >= numOfProcs -1) {
            for (int proc = 0; proc < numOfProcs; ++proc) {
                assert((state[proc] & FINISHED) == FINISHED);
            }
            break;
        }
    }
    
    // wait for my sending
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc == rank) {
            continue;
        }

        rComm.wait(numOfIndeces[proc]);
        const int size = numOfIndeces[proc];
        if (size > 0) {
            rComm.wait(&(indexArrays[proc][0]));
            rComm.wait(&(values[proc][0]));
        }
    }
}


void TlDistributeMatrix::mergeSparseMatrixAsync(const TlSparseMatrix* pMatrix, bool isFinalize)
{
    if (pMatrix != NULL) {
        assert(pMatrix->getNumOfRows() == this->getNumOfRows());
        assert(pMatrix->getNumOfCols() == this->getNumOfCols());
        
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int numOfProcs = rComm.getNumOfProc();
        
        // 送信すべきインデックスリストの作成
        std::vector<std::vector<index_type> > indexArrays(numOfProcs);
        std::vector<std::vector<double> > values(numOfProcs);
        
        TlSparseMatrix::const_iterator itEnd = pMatrix->end();
        for (TlSparseMatrix::const_iterator it = pMatrix->begin(); it != itEnd; ++it) {
            index_type globalRow = it->first.row;
            index_type globalCol = it->first.col;
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

        this->mergeMatrixAsync_send(indexArrays, values);
    }

    this->mergeMatrixAsync_recv(isFinalize);
}


void TlDistributeMatrix::mergeMatrixAsync_send(const std::vector<std::vector<index_type> >& indexArrays,
                                               const std::vector<std::vector<double> >& values)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    std::vector<int> sessionIds(numOfProcs);
    std::vector<std::size_t> numOfContents(numOfProcs);
    
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc == rank) {
            continue;
        }

        const std::size_t size = values[proc].size();
        assert(indexArrays[proc].size() == size * 2);
        numOfContents[proc] = size;

        const int sessionId = 0;
        sessionIds[proc] = sessionId;
        
        if (size > 0) {
            rComm.iSendData(sessionIds[proc], proc, TAG_MERGE_MATRIX_SESSION_ID);
            rComm.iSendData(numOfContents[proc], proc, TAG_MERGE_MATRIX_NUM_OF_CONTENTS + sessionId);
            rComm.iSendDataX(&(indexArrays[proc][0]), size * 2, proc, TAG_MERGE_MATRIX_INDECES + sessionId);
            rComm.iSendDataX(&(values[proc][0]), size, proc, TAG_MERGE_MATRIX_VALUES + sessionId);
        }
    }

    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc == rank) {
            continue;
        }

        rComm.wait(sessionIds[proc]);
        rComm.wait(numOfContents[proc]);
        rComm.wait(&(indexArrays[proc][0]));
        rComm.wait(&(values[proc][0]));
    }
}


void TlDistributeMatrix::mergeMatrixAsync_recv(bool isFinalize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    //const int numOfProcs = rComm.getNumOfProc();
    //const int rank = rComm.getRank();

    // セッション待ち受け
    int recvProc = 0;
#pragma omp critical(TlDistributeMatrix__mergeMatrixAsync_recv)
    {
    while (true) {
        if (this->mm_isWaitingSessionId_ != true) {
            rComm.iReceiveDataFromAnySource(this->mm_waitingSessionId_, TAG_MERGE_MATRIX_SESSION_ID);
            this->mm_isWaitingSessionId_ = true;
        }
        if (rComm.test(this->mm_waitingSessionId_, &recvProc) == true) {
            rComm.wait(this->mm_waitingSessionId_);
            this->mm_isWaitingSessionId_ = false;
            // task登録
            MergeMatrixRecvTask task(recvProc, this->mm_waitingSessionId_);
            this->mergeMatrixRecvTasks_.push_back(task);
        } else {
            break;
        }
    }
    
    MergeMatrixRecvTasks::iterator itEnd = this->mergeMatrixRecvTasks_.end();
    for (MergeMatrixRecvTasks::iterator it = this->mergeMatrixRecvTasks_.begin(); it != itEnd; ) {
        const int src = it->srcProc;
        if ((it->state & MM_RECV_NUM_OF_CONTENTS) == 0) {
            rComm.iReceiveData(it->numOfContents, src,
                               TAG_MERGE_MATRIX_NUM_OF_CONTENTS + it->sessionId);
            it->state |= MM_RECV_NUM_OF_CONTENTS;
        }
        if (((it->state & MM_RECV_NUM_OF_CONTENTS) != 0) &&
            ((it->state & MM_WAIT_NUM_OF_CONTENTS) == 0)) {
            if (rComm.wait(it->numOfContents) == true) {
                rComm.test(it->numOfContents);
                it->indeces.resize(it->numOfContents * 2);
                it->values.resize(it->numOfContents);
                it->state |= MM_WAIT_NUM_OF_CONTENTS;
            }
        }

        if (((it->state & MM_WAIT_NUM_OF_CONTENTS) != 0) &&
            ((it->state & MM_RECV_INDECES) == 0)) {
            rComm.iReceiveDataX(&(it->indeces[0]), it->numOfContents * 2, src,
                                TAG_MERGE_MATRIX_INDECES + it->sessionId);
            it->state |= MM_RECV_INDECES;
        }
        if (((it->state & MM_RECV_INDECES) != 0) &&
            ((it->state & MM_WAIT_INDECES) == 0)) {
            if (rComm.test(&(it->indeces[0])) == true) {
                rComm.wait(&(it->indeces[0]));
                it->state |= MM_WAIT_INDECES;
            }
        }

        if (((it->state & MM_WAIT_NUM_OF_CONTENTS) != 0) &&
            ((it->state & MM_RECV_VALUES) == 0)) {
            rComm.iReceiveDataX(&(it->values[0]), it->numOfContents, src,
                                TAG_MERGE_MATRIX_VALUES + it->sessionId);
            it->state |= MM_RECV_VALUES;
        }
        if (((it->state & MM_RECV_VALUES) != 0) &&
            ((it->state & MM_WAIT_VALUES) == 0)) {
            if (rComm.test(&(it->values[0])) == true) {
                rComm.wait(&(it->values[0]));
                it->state |= MM_WAIT_VALUES;
            }
        }

        if (((it->state & MM_WAIT_INDECES) != 0) &&
            ((it->state & MM_WAIT_VALUES)  != 0)) {
            const std::size_t numOfContents = it->numOfContents;
            for (std::size_t i = 0; i < numOfContents; ++i) {
                const index_type globalRow = it->indeces[i * 2    ];
                const index_type globalCol = it->indeces[i * 2 + 1];
                const double value = it->values[i];
                const int index = this->getIndex(globalRow, globalCol);
                assert(index != -1);
                assert((std::size_t)index < this->getNumOfMyElements());
                this->pData_[index] += value;
            }
            it->state |= MM_FINISHED;
        }

        if ((it->state & MM_FINISHED) != 0) {
            it = this->mergeMatrixRecvTasks_.erase(it);
            continue;
        }

        ++it;
    }
    }
    
    if (isFinalize == true) {
        if (this->mm_isWaitingSessionId_ == true) {
            rComm.cancel(this->mm_waitingSessionId_);
            this->mm_isWaitingSessionId_ = false;
        }
        rComm.checkNonBlockingCommunications();
    }
}


double TlDistributeMatrix::trace() const
{
    assert(this->m_nRows == this->m_nCols);

    double answer = 0.0;
    const int nGlobalCols = this->m_nCols;
    const int nRows = this->m_nMyRows;
    const int nCols = this->m_nMyCols;
    for (int c = 0; c < nCols; ++c) {
        const int nGlobalColIndex = this->m_ColIndexTable[c];
        if (nGlobalColIndex >= nGlobalCols) {
            continue;
        }

        for (int r = 0; r < nRows; ++r) {
            const int nGlobalRowIndex = this->m_RowIndexTable[r];
            if (nGlobalRowIndex == nGlobalColIndex) {
                const int index = r +  nRows * c; // row-major

                answer += this->pData_[index];
            }
        }
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);

    return answer;
}


TlDistributeMatrix& TlDistributeMatrix::operator=(const TlDistributeMatrix& rhs)
{
    if (&rhs != this) {
        assert(rhs.getNumOfRows() > 0);
        assert(rhs.getNumOfCols() > 0);
        this->m_nRows = rhs.m_nRows;
        this->m_nCols = rhs.m_nCols;
        this->m_nBlockSize = rhs.m_nBlockSize;

        this->initialize();
        std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(), this->pData_);
    }

    return (*this);
}


TlDistributeMatrix& TlDistributeMatrix::operator=(const TlDistributeSymmetricMatrix& rhs)
{
    *this = TlDistributeMatrix(rhs);
    return (*this);
}


void TlDistributeMatrix::initialize()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();

    TlScalapackContext::getData(this->m_nContext, this->m_nProc, this->m_nRank,
                                this->m_nProcGridRow, this->m_nProcGridCol);

    // my process position on the process matrix
    Cblacs_gridinfo(this->m_nContext, &(this->m_nProcGridRow), &(this->m_nProcGridCol),
                    &(this->m_nMyProcRow), &(this->m_nMyProcCol));

    // determine sizes of local matrix
    const int nStartRowProc = 0;
    const int nStartColProc = 0;
    this->m_nMyRows = std::max(1, numroc_(&(this->m_nRows), &(this->m_nBlockSize),
                                          &(this->m_nMyProcRow), &nStartRowProc, &(this->m_nProcGridRow)));
    this->m_nMyCols = std::max(1, numroc_(&(this->m_nCols), &(this->m_nBlockSize),
                                          &(this->m_nMyProcCol), &nStartColProc, &(this->m_nProcGridCol)));

    // make parameter, desca
    int nInfo = 0;
    descinit_(this->m_pDESC, &(this->m_nRows), &(this->m_nCols), &(this->m_nBlockSize), &(this->m_nBlockSize),
              &nStartRowProc, &nStartColProc, &(this->m_nContext), &(this->m_nMyRows), &nInfo);
    assert(nInfo == 0);

    // 行列データ用バッファメモリの確保
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
    this->pData_ = new double[this->getNumOfMyElements()];
    std::fill(this->pData_, this->pData_ + this->getNumOfMyElements(), 0.0);
    
    // 行方向のglobal_index v.s. local_indexのリストを作成
    {
        const int nMyRows = this->m_nMyRows;
        const int nBlockSize = this->m_nBlockSize;
        const int nBlockIndex = this->m_nMyProcRow * nBlockSize; // 各ローカル行列の最初のインデックス
        const int nIncrementBlockIndex = this->m_nProcGridRow * nBlockSize; // ブロック最初のインデックスの増分
        this->m_RowIndexTable.clear();
        this->m_RowIndexTable.reserve(nMyRows);
        for (int r = 0; r < nMyRows; ++r) {
            const div_t d = std::div(r, nBlockSize);
            const int i = nBlockIndex +(nIncrementBlockIndex*d.quot) +d.rem;
            if (i < this->m_nRows) {
                this->m_RowIndexTable.push_back(i);
            } else {
                break;
            }
        }
        std::vector<int>(this->m_RowIndexTable).swap(this->m_RowIndexTable);
    }

    // 列方向のglobal_index v.s. local_indexのリストを作成
    {
        const int nMyCols = this->m_nMyCols;
        const int nBlockSize = this->m_nBlockSize;
        const int nBlockIndex = this->m_nMyProcCol * this->m_nBlockSize; // 各ローカル行列の最初のインデックス
        const int nIncrementBlockIndex = this->m_nProcGridCol * this->m_nBlockSize; // ブロック最初のインデックスの増分
        this->m_ColIndexTable.clear();
        this->m_ColIndexTable.reserve(nMyCols);
        for (int c = 0; c < nMyCols; ++c) {
            const div_t d = std::div(c, nBlockSize);
            const int i = nBlockIndex +(nIncrementBlockIndex*d.quot) +d.rem;
            if (i < this->m_nCols) {
                this->m_ColIndexTable.push_back(i);
            } else {
                break;
            }
        }
        std::vector<int>(this->m_ColIndexTable).swap(this->m_ColIndexTable);
    }

    // getPartialMatrix通信用 初期化処理
    this->isDebugOut_GPM_ = false;

#pragma omp critical (TlDistributeMatrix_gpmServerTasks_)
    {
        this->gpmServerTasks_.clear();
    }
    this->isWaitingRequestHandShake_ = false;
    this->sessionId_ = 0;

    this->trafficControl_.clear();
    this->trafficControl_.resize(numOfProcs);
    for (int i = 0; i < numOfProcs; ++i) {
        this->trafficControl_[i] = NULL;
    }

    this->sessionTable_.resize(numOfProcs);
    for (int proc = 0; proc < numOfProcs; ++proc) {
        this->sessionTable_[proc].reset();
    }
    
    this->gpmClientTasks_.clear();

    // mergeMatrix用
    this->mm_isWaitingSessionId_ = false;
}


std::size_t TlDistributeMatrix::getNumOfMyElements() const {
    return (this->m_nMyRows * this->m_nMyCols);
}



TlMatrixObject::index_type TlDistributeMatrix::getNumOfRows() const
{
    return this->m_nRows;
}


TlMatrixObject::index_type TlDistributeMatrix::getNumOfCols() const
{
    return this->m_nCols;
}


int TlDistributeMatrix::getBlockSize() const {
    return this->m_nBlockSize;
}


int TlDistributeMatrix::getIndex(const index_type nGlobalRow, const index_type nGlobalCol) const
{
    int nAnswer = -1;

    std::vector<int>::const_iterator pRow = std::lower_bound(this->m_RowIndexTable.begin(), this->m_RowIndexTable.end(), nGlobalRow);
    if ((pRow != this->m_RowIndexTable.end()) && (*pRow == nGlobalRow)) {
        std::vector<int>::const_iterator pCol = std::lower_bound(this->m_ColIndexTable.begin(), this->m_ColIndexTable.end(), nGlobalCol);

        if ((pCol != this->m_ColIndexTable.end()) && (*pCol == nGlobalCol)) {
            const int localRow = pRow - this->m_RowIndexTable.begin();
            const int localCol = pCol - this->m_ColIndexTable.begin();

            nAnswer = localRow + localCol * (this->m_nMyRows);
        }
    }

    return nAnswer;
}


int TlDistributeMatrix::getProcIdForIndex(const index_type globalRow, const index_type globalCol) const
{
     const int rowBlockIndex = globalRow / this->m_nBlockSize;
     const int rowProcID = rowBlockIndex % this->m_nProcGridRow;
     const int colBlockIndex = globalCol / this->m_nBlockSize;
     const int colProcID = colBlockIndex % this->m_nProcGridCol;

     // "Row-major" for Cblacs_gridinit
     const int procMatrixIndex = rowProcID * this->m_nProcGridCol + colProcID;
     //const int procID = this->processMatrix_[procMatrixIndex];
     int procID = procMatrixIndex;

     return procID;
}


double TlDistributeMatrix::get(const index_type row, const index_type col) const
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //DataType& data_tmp = const_cast<DataType&>(this->data_);

    double dAnswer = 0.0;
    const int nGlobalRow = row +1;
    const int nGlobalCol = col +1;
    pdelget_("A", " ", &dAnswer, this->pData_, &nGlobalRow, &nGlobalCol, this->m_pDESC);

    return dAnswer;
}

void TlDistributeMatrix::set(const index_type row, const index_type col, const double dValue)
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));

    const int index = this->getIndex(row, col);
    if (index != -1) {
        this->pData_[index] = dValue;
    }
}


void TlDistributeMatrix::add(const index_type row, const index_type col, const double value)
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));

    const int index = this->getIndex(row, col);
    if (index != -1) {
        this->pData_[index] += value;
    }
}


void TlDistributeMatrix::addByList(const index_type* pIndexPairs,
                                   const double* pValues,
                                   const std::size_t size)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    rComm.checkNonBlockingCommunications();

    // 送信すべきインデックスリストの作成
    std::vector<std::vector<index_type> > sendIndexArrays(numOfProcs);
    std::vector<std::vector<double> > sendValues(numOfProcs);

    for (std::size_t i = 0; i < size; ++i) {
        const index_type globalRow = pIndexPairs[i*2   ];
        const index_type globalCol = pIndexPairs[i*2 +1];
        const double value = pValues[i];

        const int index = this->getIndex(globalRow, globalCol);
        if (index != -1) {
            this->pData_[index] += value;
        } else {
            const int targetProc = this->getProcIdForIndex(globalRow, globalCol);
            sendIndexArrays[targetProc].push_back(globalRow);
            sendIndexArrays[targetProc].push_back(globalCol);
            sendValues[targetProc].push_back(value);
        }
    }

    this->mergeMatrix_common(sendIndexArrays, sendValues);
    rComm.checkNonBlockingCommunications();
}


double TlDistributeMatrix::getLocal(const index_type row, const index_type col) const
{
    double answer = 0.0;
    const int index = this->getIndex(row, col);
    if (index != -1) {
        answer = this->pData_[index];
    }

    return answer;
}


std::vector<TlDistributeMatrix::index_type>
TlDistributeMatrix::getRowIndexTable() const
{
    return this->m_RowIndexTable;
}


std::vector<TlDistributeMatrix::index_type>
TlDistributeMatrix::getColIndexTable() const
{
    return this->m_ColIndexTable;
}


TlMatrix TlDistributeMatrix::getLocalMatrix() const
{
    const std::vector<index_type> rowIndexes = this->getRowIndexTable();
    const std::vector<index_type> colIndexes = this->getColIndexTable();
    const int numOfRowIndexes = rowIndexes.size();
    const int numOfColIndexes = colIndexes.size();
    TlMatrix answer(numOfRowIndexes, numOfColIndexes);
    
    for (int rowIndex = 0; rowIndex < numOfRowIndexes; ++rowIndex) {
        const index_type row = rowIndexes[rowIndex];
        for (int colIndex = 0; colIndex < numOfColIndexes; ++colIndex) {
            const index_type col = colIndexes[colIndex];
            answer.set(rowIndex, colIndex, this->getLocal(row, col));
        }
    }

    return answer;
}


TlDistributeMatrix& TlDistributeMatrix::operator+=(const TlDistributeMatrix& rhs)
{
    assert(this->m_nRows == rhs.m_nRows);
    assert(this->m_nCols == rhs.m_nCols);
    assert(this->m_nBlockSize == rhs.m_nBlockSize);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] += rhs.pData_[i];
    }

    return (*this);
}


TlDistributeMatrix& TlDistributeMatrix::operator-=(const TlDistributeMatrix& rhs)
{
    assert(this->m_nRows == rhs.m_nRows);
    assert(this->m_nCols == rhs.m_nCols);
    assert(this->m_nBlockSize == rhs.m_nBlockSize);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] -= rhs.pData_[i];
    }

    return (*this);
}


TlDistributeMatrix& TlDistributeMatrix::operator*=(const TlDistributeMatrix& rhs)
{
    TlDistributeMatrix tmp = *this;
    *this = tmp * rhs;

    return (*this);
}


TlDistributeMatrix& TlDistributeMatrix::operator*=(const TlDistributeSymmetricMatrix& rhs)
{
    TlDistributeMatrix tmp = *this;
    *this = tmp * rhs;

    return (*this);
}


TlDistributeMatrix& TlDistributeMatrix::operator*=(const double dCoef)
{
    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] *= dCoef;
    }

    return (*this);
}

TlDistributeMatrix& TlDistributeMatrix::operator/=(const double dCoef)
{
    return this->operator*=(1.0 / dCoef);
}


TlDistributeMatrix operator*(const TlDistributeMatrix& X, const TlDistributeMatrix& Y)
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //TlDistributeMatrix& Xtmp = const_cast<TlDistributeMatrix&>(X);
    //TlDistributeMatrix& Ytmp = const_cast<TlDistributeMatrix&>(Y);

    const int nXRow = X.getNumOfRows();
    const int nXCol = X.getNumOfCols();
    const int nYCol = Y.getNumOfCols();
    assert(X.getNumOfCols() == Y.getNumOfRows());

    const int nZRow = nXRow;
    const int nZCol = nYCol;
    TlDistributeMatrix Z(nZRow, nZCol);

    const double dAlpha = 1.0;
    const double dBeta = 1.0;
    const int nIX = 1;
    const int nJX = 1;
    const int nIY = 1;
    const int nJY = 1;
    const int nIZ = 1;

    pdgemm_("N","N",
            &nZRow, &nZCol, &nXCol, &dAlpha,
            X.pData_, &nIX, &nJX, X.m_pDESC,
            Y.pData_, &nIY, &nJY, Y.m_pDESC,
            &dBeta,
            Z.pData_, &nIZ, &nIZ, Z.m_pDESC);

    return Z;
}


TlDistributeVector operator*(const TlDistributeMatrix& A, const TlDistributeVector& X)
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //TlDistributeMatrix& Atmp = const_cast<TlDistributeMatrix&>(A);
    TlDistributeVector& Xtmp = const_cast<TlDistributeVector&>(X);

    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();

    assert(N == X.getSize());

    TlDistributeVector Y(M);

    const char TRANS = 'N';
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

    pdgemv_(&TRANS, &M, &N,
            &alpha, A.pData_, &IA, &JA, A.m_pDESC,
            &(Xtmp.data_[0]), &IX, &JX, X.m_pDESC, &INCX,
            &beta, &(Y.data_[0]), &IY, &JY, Y.m_pDESC, &INCY);

    return Y;
}


// TlDistributeVector operator*(const TlDistributeVector& X, const TlDistributeMatrix& A)
// {
//     // this const_cast is requiered for PGI compiler
//     // "error: expression must be an lvalue or a function designator"
//     TlDistributeMatrix& Atmp = const_cast<TlDistributeMatrix&>(A);
//     TlDistributeVector& Xtmp = const_cast<TlDistributeVector&>(X);

//     const int M = A.getNumOfRows();
//     const int N = A.getNumOfCols();

//     assert(M == X.getSize());

//     TlDistributeVector Y(N);

//     const char TRANS = 'T';
//     const double alpha = 1.0;
//     const double beta = 0.0;
//     const int IA = 1;
//     const int JA = 1;
//     const int IX = 1;
//     const int JX = 1;
//     const int INCX = 1;
//     const int IY = 1;
//     const int JY = 1;
//     const int INCY = 1;

//     pdgemv_(&TRANS, &M, &N,
//             &alpha, &(Atmp.data_[0]), &IA, &JA, A.m_pDESC,
//             &(Xtmp.data_[0]), &IX, &JX, X.m_pDESC, &INCX,
//             &beta, &(Y.data_[0]), &IY, &JY, Y.m_pDESC, &INCY);

//     return Y;
// }


TlDistributeMatrix operator*(const double X, const TlDistributeMatrix& Y)
{
    TlDistributeMatrix ans = Y;
    ans *= X;

    return ans;
}

TlDistributeMatrix operator*(const TlDistributeMatrix& X, const double Y)
{
    TlDistributeMatrix ans = X;
    ans *= Y;

    return ans;
}

TlDistributeMatrix operator+(const TlDistributeMatrix& X, const TlDistributeMatrix& Y)
{
    TlDistributeMatrix answer = X;
    answer += Y;

    return answer;
}

TlDistributeMatrix operator-(const TlDistributeMatrix& X, const TlDistributeMatrix& Y)
{
    TlDistributeMatrix answer = X;
    answer -= Y;

    return answer;
}

double TlDistributeMatrix::operator()(const int row, const int col) const
{
    return this->get(row, col);
}


double& TlDistributeMatrix::operator()(const index_type row, const index_type col)
{
    this->m_dTempVar = this->get(row, col);
    const int index = this->getIndex(row, col);
    if (index != -1) {
        return this->pData_[index];
    } else {
        return this->m_dTempVar; // 使い捨て
    }
}


bool TlDistributeMatrix::inverse()
{
#ifdef HAVE_SCALAPACK
    // using SCALAPACK
    return inverseByScaLapack(*this);
#else
    // without SCALAPACK
#error NOT found algebra package: need SCALAPACK library
#endif // HAVE_SCALAPACK  
}


TlVector TlDistributeMatrix::getDiagonalElements() const
{
    const index_type dim = std::min(this->getNumOfRows(), this->getNumOfCols());
    TlVector answer(dim);

    for (index_type i = 0; i < dim; ++i) {
        const double value = this->getLocal(i, i);
        answer.set(i, value);
    }
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);
    
    return answer;
}


// pddot?
const TlDistributeMatrix& TlDistributeMatrix::dot(const TlDistributeMatrix& X)
{
    assert(this->getNumOfRows() == X.getNumOfRows());
    assert(this->getNumOfCols() == X.getNumOfCols());

    const size_type size = this->getNumOfMyElements();

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
                this->pData_[index] *= X.pData_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->pData_[index] *= X.pData_[index];
        }
    }
#else
    // not use OpenMP
    for (size_type index = 0; index < size; ++index) {
        this->pData_[index] *= X.pData_[index];
    }
#endif // _OPENMP
    
    return (*this);
}


double TlDistributeMatrix::sum() const
{
    double answer = 0.0;
    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        answer += this->pData_[i];
    }
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);
    return answer;
}


bool TlDistributeMatrix::isLoadable(std::ifstream& ifs)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool answer = false;

    if (rComm.isMaster() == true) {
        answer = TlMatrix::isLoadable(ifs);
    }
    rComm.broadcast(answer);
    return answer;
}

bool TlDistributeMatrix::isLoadable(const std::string& rFilePath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool answer = false;

    if (rComm.isMaster() == true) {
        answer = TlMatrix::isLoadable(rFilePath);
    }
    rComm.broadcast(answer);
    return answer;
}

bool TlDistributeMatrix::load(const std::string& sFilePath)
{
    if (TlDistributeMatrix::isUsingPartialIO == true) {
        return this->loadLocal(sFilePath);
    }
    
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
        if (rComm.isMaster() == true) {
            std::cerr << "[error] could not open file. " << sFilePath << std::endl;
        }
        abort();
    }

    bAnswer = this->load(ifs);

    if (bAnswer != true) {
        if (rComm.isMaster() == true) {
            std::cerr << "TlDistributeMatrix::load() is not supported: " << sFilePath << std::endl;
        }
        std::abort();
    }

    if (rComm.isMaster() == true) {
        ifs.close();
    }

    return bAnswer;
}


bool TlDistributeMatrix::load(std::ifstream& ifs)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());

    const int numOfProcs = rComm.getNumOfProc();

    bool bAnswer = true;
    enum {RSFD, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};

    // read header
    int nType = 0;
    if (rComm.isMaster() == true) {
        ifs.read((char*)&(nType),         sizeof(int));
        ifs.read((char*)&(this->m_nRows), sizeof(int));
        ifs.read((char*)&(this->m_nCols), sizeof(int));
    }
    rComm.broadcast(bAnswer);
    rComm.broadcast(nType);

    switch (nType) {
    case RSFD:
        //     std::cerr << "load RSFD" << std::endl;
        break;
    default:
        std::cerr << "this matrix type is not supported. stop." << std::endl;
        bAnswer = false;
        break;
    }

    rComm.broadcast(this->m_nRows);
    rComm.broadcast(this->m_nCols);
    this->initialize();

    if (rComm.isMaster() == true) {
        static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
        std::vector<double> buf(bufferCount, 0.0);

        index_type currentRow = 0;
        index_type currentCol = 0;
        const index_type maxRow = this->m_nRows;
        const index_type maxCol = this->m_nCols;
        bool isFinished = false;

        std::vector<int> sizeLists(numOfProcs, 0);
        std::vector<std::vector<index_type> > rowColLists(numOfProcs);
        std::vector<std::vector<double> > valueLists(numOfProcs);
        std::vector<bool> isSendData(numOfProcs, false);
        
        while (isFinished == false) {
            // buffer分を一度に読み込み
            ifs.read((char*)(&buf[0]), sizeof(double) * bufferCount);

            // 各プロセスのバッファに振り分ける
            std::map<int, std::vector<index_type> > tmpRowColLists;
            std::map<int, std::vector<double> > tmpValueLists;
            for (std::size_t i = 0; i < bufferCount; ++i) {
                const int proc = this->getProcIdForIndex(currentRow, currentCol);
                if (proc == 0) {
                    // masterが持っている
                    const size_type index = this->getIndex(currentRow, currentCol);
                    assert(index != -1);
                    this->pData_[index] = buf[i];
                } else {
                    tmpRowColLists[proc].push_back(currentRow);
                    tmpRowColLists[proc].push_back(currentCol);
                    tmpValueLists[proc].push_back(buf[i]);
                }
                
                // count up
                ++currentCol;
                if (currentCol >= maxCol) {
                    currentCol = 0;
                    ++currentRow;
                    if (currentRow >= maxRow) {
                        isFinished = true;
                        break;
                    }
                }
            }

            // データを送信
            std::map<int, std::vector<index_type> >::const_iterator itEnd = tmpRowColLists.end();
            for (std::map<int, std::vector<index_type> >::const_iterator it = tmpRowColLists.begin(); it != itEnd; ++it) {
                const int proc = it->first;
                const int numOfContents = it->second.size() / 2;
                assert(std::size_t(numOfContents) == tmpValueLists[proc].size());

                if (isSendData[proc] == true) {
                    rComm.wait(sizeLists[proc]);
                    if (sizeLists[proc] > 0) {
                        rComm.wait(&(rowColLists[proc][0]));
                        rComm.wait(&(valueLists[proc][0]));
                    }
                    isSendData[proc] = false;
                }
                
                sizeLists[proc] = numOfContents;
                rowColLists[proc] = tmpRowColLists[proc];
                valueLists[proc] = tmpValueLists[proc];
                
                rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
                if (sizeLists[proc] > 0) {
                    rComm.iSendDataX(&(rowColLists[proc][0]), (sizeLists[proc] * 2), proc, TAG_LOAD_ROWCOLS);
                    rComm.iSendDataX(&(valueLists[proc][0]), sizeLists[proc], proc, TAG_LOAD_VALUES);
                }
                isSendData[proc] = true;
            }
        } // end while

        for (int proc = 1; proc < numOfProcs; ++proc) { // proc == 0 は送信しない
            if (isSendData[proc] == true) {
                rComm.wait(sizeLists[proc]);
                rComm.wait(&(rowColLists[proc][0]));
                rComm.wait(&(valueLists[proc][0]));
                isSendData[proc] = false;
            }
        }
        
        // 終了メッセージを全ノードに送る
        std::vector<int> endMsg(numOfProcs, 0);
        for (int proc = 1; proc < numOfProcs; ++proc) { // proc == 0 は送信しない
            rComm.iSendData(endMsg[proc], proc, TAG_LOAD_END);
        }      
        for (int proc = 1; proc < numOfProcs; ++proc) { // proc == 0 は送信しない
            rComm.wait(endMsg[proc]);
        }      
    } else {
        // slave
        const int root = 0;
        int sizeList = 0;
        std::vector<index_type> rowColList;
        std::vector<double> valueList;
        int endMsg = 0;
        
        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
        rComm.iReceiveData(endMsg, root, TAG_LOAD_END);
        bool isLoopBreak = false;
        while (isLoopBreak == false) {
            if (rComm.test(sizeList) == true) {
                rComm.wait(sizeList);
//                 std::cerr << TlUtils::format("RECV [%d] size=%d",
//                                              rComm.getRank(),
//                                              sizeList)
//                           << std::endl;
                if (sizeList > 0) {
                    rowColList.resize(sizeList * 2);
                    valueList.resize(sizeList);
                    rComm.receiveDataX(&(rowColList[0]), sizeList * 2, root, TAG_LOAD_ROWCOLS);
                    rComm.receiveDataX(&(valueList[0]), sizeList, root, TAG_LOAD_VALUES);

                    for (int i = 0; i < sizeList; ++i) {
                        const index_type row = rowColList[i * 2    ];
                        const index_type col = rowColList[i * 2 + 1];
                        const size_type index = this->getIndex(row, col);
                        
//                         std::cerr << TlUtils::format("RECV [%d] (%4d, %4d)",
//                                                      rComm.getRank(),
//                                                      row, col) << std::endl;
                        
                        assert(index != -1);
                        if (index == -1) {
                            abort();
                        }
                        this->pData_[index] = valueList[i];
                    }
                }

                rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
            }

            if (rComm.test(endMsg) == true) {
                rComm.wait(endMsg);
                rComm.cancel(sizeList);
//                 std::cerr << TlUtils::format("RECV [%d] END",
//                                                  rComm.getRank())
//                               << std::endl;
                isLoopBreak = true;
            }
        }
    }

    assert(rComm.checkNonBlockingCommunications());
    return bAnswer;
}


bool TlDistributeMatrix::save(const std::string& sFilePath) const
{
    //std::cerr << "TlDistributeMatrix::save() file=" << sFilePath << std::endl;
    if (TlDistributeMatrix::isUsingPartialIO == true) {
        return this->saveLocal(sFilePath);
    }
    
    bool bAnswer = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());

    if (rComm.isMaster() == true) {
        // master
        const int nGlobalRows = this->m_nRows;
        const int nGlobalCols = this->m_nCols;
        TlFileMatrix fm(sFilePath, nGlobalRows, nGlobalCols);

        // store local matrix
        {
            const int nRows = this->m_nMyRows;
            const int nCols = this->m_nMyCols;
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

                    const int index = r +  nRows * c; // row-major
                    fm.set(nGlobalRowIndex, nGlobalColIndex, this->pData_[index]);
                }
            }
        }

        // recive submatrix & write
        const int numOfSlaves = rComm.getNumOfProc() -1;
        for (int i = 0; i < numOfSlaves; ++i) {
//       std::cerr << "TlDistributeMatrix::save() loop i = " << i << std::endl;
            int msg = 0;
            int src = 0;
            rComm.receiveDataFromAnySource(msg, &src, TAG_SAVE_HANDSHAKE); // どのプロセスから依頼があったかを調査
            assert(msg == 0);
//       std::cerr << "request from " << src << std::endl;
            msg = 0;
            rComm.sendData(msg, src, TAG_SAVE_HANDSHAKE_OK); // そのプロセスにデータ送信を依頼

            index_type nRows;
            index_type nCols;
            std::vector<index_type> rowIndexTable;
            std::vector<index_type> colIndexTable;

            rComm.receiveData(nRows, src, TAG_SAVE_DATA_ROWS);
            rComm.receiveData(nCols, src, TAG_SAVE_DATA_COLS);
            rComm.receiveData(rowIndexTable, src, TAG_SAVE_DATA_ROWINDEX);
            rComm.receiveData(colIndexTable, src, TAG_SAVE_DATA_COLINDEX);
            const std::size_t bufSize = nRows * nCols;
            double* pBuf = new double[bufSize];
            rComm.iReceiveDataX(pBuf, bufSize, src, TAG_SAVE_DATA);
            rComm.wait(pBuf);

            const int nMaxRows = rowIndexTable.size();
            const int nMaxCols = colIndexTable.size();
            for (int r = 0; r < nMaxRows; ++r) {
                const int nGlobalRowIndex = rowIndexTable[r];
                if (nGlobalRowIndex >= nGlobalRows) {
                    continue;
                }

                for (int c = 0; c < nMaxCols; ++c) {
                    const int nGlobalColIndex = colIndexTable[c];
                    if (nGlobalColIndex >= nGlobalCols) {
                        continue;
                    }

                    const int index = r + nRows * c; // row-major
                    fm.set(nGlobalRowIndex, nGlobalColIndex, pBuf[index]);
                }
            }

            delete[] pBuf;
            pBuf = NULL;
        }
    } else {
        // slave: send submatrix
        int msg = 0;
        const int root = 0;
        rComm.sendData(msg, root, TAG_SAVE_HANDSHAKE);
        rComm.receiveData(msg, root, TAG_SAVE_HANDSHAKE_OK);
        assert(msg == 0);

        rComm.sendData(this->m_nMyRows, root, TAG_SAVE_DATA_ROWS);
        rComm.sendData(this->m_nMyCols, root, TAG_SAVE_DATA_COLS);
        rComm.sendData(this->m_RowIndexTable, root, TAG_SAVE_DATA_ROWINDEX);
        rComm.sendData(this->m_ColIndexTable, root, TAG_SAVE_DATA_COLINDEX);
        rComm.iSendDataX(this->pData_, this->getNumOfMyElements(), root, TAG_SAVE_DATA);
        rComm.wait(this->pData_);
    }

    rComm.broadcast(bAnswer);
    assert(rComm.checkNonBlockingCommunications());
    return bAnswer;
}


bool TlDistributeMatrix::saveLocal(const std::string& filePath) const
{
    //std::cerr << "TlDistributeMatrix::saveLocal() file=" << filePath << std::endl;
    
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //DataType& data_tmp = const_cast<DataType&>(this->data_);

    bool answer = true;
    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

    const int type = 16; // means Distributed(16) + RSFD(0)
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
    ofs.write(reinterpret_cast<const char*>(this->pData_), sizeof(double) * this->getNumOfMyElements());
    
    ofs.close();

    //std::cerr << "save: rank=" << this->m_nRank << std::endl;
    //std::cerr << "save: proc=" << this->m_nProc << std::endl;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.barrier();
    
    return answer;
}


bool TlDistributeMatrix::loadLocal(const std::string& filePath)
{
    //std::cerr << "TlDistributeMatrix::loadLocal() file=" << filePath << std::endl;
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

        if (type != 16) {
            std::cerr << TlUtils::format("file type mismatch(%d). ", type)
                      << __FILE__ << "," << __LINE__
                      << std::endl;
        }
        this->m_nRows = globalRow;
        this->m_nCols = globalCol;
        //std::cerr << "rank = " << rank << std::endl;
        //std::cerr << "proc = " << proc << std::endl;
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


bool TlDistributeMatrix::saveText(const std::string& sFilePath) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    //this->update();
    bool bAnswer = true;

    std::ofstream ofs;
    if (rComm.isMaster() == true) {
        ofs.open(sFilePath.c_str(), std::ofstream::out);
    }

    //if (ofs.good()){
    bAnswer = this->saveText(ofs);
    //} else {
    //bAnswer =false;
    //}

    if (rComm.isMaster() == true) {
        ofs.close();
    }

    return bAnswer;
}


bool TlDistributeMatrix::saveText(std::ofstream& ofs) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    //this->update();
    bool bAnswer = true;
    const int nRows = this->getNumOfRows();
    const int nCols = this->getNumOfCols();

    if (rComm.isMaster() == true) {
        ofs << "TEXT\n";
        ofs << nRows << "\n";
        ofs << nCols << "\n";
    }

    // print out LCAO coefficent
    for (int i=0; i < nRows; ++i) {
        for (int j=0; j < nCols; ++j) {
            const double tmp = this->get(i, j);
            if (rComm.isMaster() == true) {
                ofs << TlUtils::format(" %10.6lf", tmp);
            }
        }
        if (rComm.isMaster() == true) {
            ofs << "\n";
        }
    }
    if (rComm.isMaster() == true) {
        ofs << "\n";
    }

    return bAnswer;
}


////////////////////////////////////////////////////////////////////////
bool inverseByScaLapack(TlDistributeMatrix& X)
{
    //TlCommunicate& rComm = TlCommunicate::getInstance();
    //const double dStartTime = rComm.getTime();

    //X.print(std::cout);

    bool bAnswer = true;

    const int M = X.getNumOfRows();
    const int N = X.getNumOfCols();
    const int IA = 1;
    const int JA = 1;
    const int zero = 0;
    int info = 0;
    const int sizeOf_IPIV = numroc_(&X.m_nRows, &X.m_nBlockSize,
                                    &X.m_nMyProcRow, &zero, &X.m_nProcGridRow)+ X.m_nBlockSize;
    //std::cout << "IPIV=" << sizeOf_IPIV << std::endl;

    int* IPIV = new int[sizeOf_IPIV];
    //   std::fill_n(IPIV, sizeOf_IPIV, 0);

    pdgetrf_(&M, &N, X.pData_, &IA, &JA, X.m_pDESC, IPIV, &info);

    if (info == 0) {
        int LWORK = -1;
        int LIWORK = -1;
        double* WORK_SIZE = new double[1];
        int* IWORK_SIZE = new int[1];
        pdgetri_(&M, X.pData_, &IA, &JA, X.m_pDESC, IPIV,
                 WORK_SIZE, &LWORK, IWORK_SIZE, &LIWORK, &info);

        LWORK = static_cast<int>(WORK_SIZE[0]);
        //std::cout << "LWORK=" << LWORK << std::endl;
        double* WORK = new double[LWORK];
        LIWORK = IWORK_SIZE[0];
        //std::cout << "LIWORK=" << LIWORK << std::endl;
        int* IWORK = new int[LIWORK];
        delete[] WORK_SIZE;
        WORK_SIZE = NULL;
        delete[] IWORK_SIZE;
        IWORK_SIZE = NULL;

        pdgetri_(&M, X.pData_, &IA, &JA, X.m_pDESC, IPIV,
                 WORK, &LWORK, IWORK, &LIWORK, &info);

        if (info != 0) {
            std::cout << "pdgetri_ returns " << info << std::endl;
            bAnswer = false;
        }

        delete[] IWORK;
        IWORK = NULL;

        delete[] WORK;
        WORK = NULL;
    } else {
        std::cout << "pdgetrf_ returns " << info << std::endl;
        bAnswer = false;
    }

    delete[] IPIV;
    IPIV = NULL;

//   const double dEndTime = rComm.getTime();
//   if (rComm.isMaster() == true){
//     std::cout << "TlDistributeMatrix::inverseMatrix() time:" << (dEndTime - dStartTime) << std::endl;
//   }

    return bAnswer;
}

const TlDistributeMatrix& TlDistributeMatrix::transpose()
{
    const int M = this->m_nCols;
    const int N = this->m_nRows;
    const double alpha = 1.0;
    const double beta = 0.0;
    const int IA = 1;
    const int JA = 1;
    TlDistributeMatrix C(M, N);
    const int IC = 1;
    const int JC = 1;
    pdtran_(&M, &N,
            &alpha, &(this->pData_[0]), &IA, &JA, this->m_pDESC,
            &beta,  &(C.pData_[0]), &IC, &JC, C.m_pDESC);
    
    (*this) = C;
    
    return (*this);
}


