#include <cassert>
#include <iostream>
#include <fstream>
#include <limits>
#include <numeric>

#include "TlMath.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "TlFileSymmetricMatrix.h"
#include "scalapack.h"
#include "TlUtils.h"
#include "TlLogging.h"

TlDistributeSymmetricMatrix::TlDistributeSymmetricMatrix(const index_type dim)
    : TlDistributeMatrix(dim, dim)
{
}


TlDistributeSymmetricMatrix::TlDistributeSymmetricMatrix(const TlDistributeSymmetricMatrix& rhs)
    : TlDistributeMatrix((TlDistributeMatrix&)rhs)
{
}


TlDistributeSymmetricMatrix::TlDistributeSymmetricMatrix(const TlDistributeMatrix& rhs)
    : TlDistributeMatrix(rhs)
{
    assert(rhs.getNumOfRows() == rhs.getNumOfCols());
    // コピーされたバッファの下半分しか使わない
}


TlDistributeSymmetricMatrix::TlDistributeSymmetricMatrix(const TlDistributeVector& rhs,
                                                         const index_type dim)
    : TlDistributeMatrix(rhs, dim, dim)
{
}


TlDistributeSymmetricMatrix::~TlDistributeSymmetricMatrix()
{
}


void TlDistributeSymmetricMatrix::resize(const index_type size)
{
    assert(size > 0);
    TlDistributeMatrix::resize(size, size);
}


TlDistributeSymmetricMatrix& TlDistributeSymmetricMatrix::operator=(const TlDistributeSymmetricMatrix& rhs)
{
    if (&rhs != this) {
        TlDistributeMatrix::operator=(rhs);
    }

    return (*this);
}


TlDistributeSymmetricMatrix& TlDistributeSymmetricMatrix::operator=(const TlDistributeMatrix& rhs)
{
    if (&rhs != this) {
        TlDistributeSymmetricMatrix tmp(rhs);
        *this = tmp;
    }
    
    return (*this);
}


double TlDistributeSymmetricMatrix::get(index_type row, index_type col) const
{
    if (row < col) {
        std::swap(row, col);
    }
    return TlDistributeMatrix::get(row, col);
}


double TlDistributeSymmetricMatrix::operator()(index_type row, index_type col) const
{
    if (row < col) {
        std::swap(row, col);
    }
    return this->get(row, col);
}


void TlDistributeSymmetricMatrix::set(index_type row, index_type col, const double value)
{
    if (row < col) {
        std::swap(row, col);
    }
    TlDistributeMatrix::set(row, col, value);
}


void TlDistributeSymmetricMatrix::add(index_type row, index_type col, const double value)
{
    if (row < col) {
        std::swap(row, col);
    }
    TlDistributeMatrix::add(row, col, value);
}


double& TlDistributeSymmetricMatrix::operator()(index_type row, index_type col)
{
    if (row < col) {
        std::swap(row, col);
    }
    return TlDistributeMatrix::operator()(row, col);
}


double TlDistributeSymmetricMatrix::getLocal(index_type row, index_type col) const
{
    if (row < col) {
        std::swap(row, col);
    }
    return TlDistributeMatrix::getLocal(row, col);
}


TlSparseSymmetricMatrix TlDistributeSymmetricMatrix::getPartialMatrix(const double threshold) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    const index_type globalSize = this->m_nRows;
    assert(this->m_nRows == this->m_nCols);

    // 自分のデータのうち、有意なデータを抽出する
    TlSparseSymmetricMatrix m(globalSize);
    {
        const index_type numOfLocalRows = this->m_nMyRows;
        const index_type numOfLocalCols = this->m_nMyCols;
        for (index_type c = 0; c < numOfLocalCols; ++c) {
            const index_type globalColIndex = this->m_ColIndexTable[c];
            const size_type index_tmp = numOfLocalRows * c;

            for (index_type r = 0; r < numOfLocalRows; ++r) {
                const index_type globalRowIndex = this->m_RowIndexTable[r];
                if (globalRowIndex < globalColIndex) {
                    // for symmetry
                    continue;
                }

                const size_type index = r + index_tmp;
                const double value = this->pData_[index];
                if (std::fabs(value) > threshold) {
                    m.set(globalRowIndex, globalColIndex, value);
                }
            }
        }
    }

    // 集計
    std::size_t numOfElements = m.getSize();
    rComm.allReduce_SUM(numOfElements);
    const std::size_t averageNumOfElements = (numOfElements + proc -1) / proc;


//   // 平均よりも多い要素を持っているSparseMatrixは他に渡す
//   std::vector<int> sendProcList(proc, 0);
//   std::vector<int> recvProcList(proc, 0);
//   TlSparseSymmetricMatrix diffM(globalSize);
//   const int maxCycle = m.getSize() - averageNumOfElements;
//   if (maxCycle > 0) {
//     sendProcList[rank] = maxCycle;
//   } else if (maxCycle < 0) {
//     recvProcList[rank] = - maxCycle;
//   }
//   for (int i = 0; i < maxCycle; ++i) {
//     int row = 0;
//     int col = 0;
//     const double value = m.pop(&row, &col);
//     diffM.set(row, col, value);
//   }

//   // diffMの集計
//   int numOfDiffMElements = diffM.getSize();
//   rComm.allReduce_SUM(numOfDiffMElements);
//   const int numOfHolds = (numOfDiffMElements + proc -1) / proc; // この数だけ担当する
//   rComm.allReduce_SUM(sendProcList);
//   rComm.allReduce_SUM(recvProcList);
//   int numOfTransportCycle = std::accumulate(sendProcList.begin(), sendProcList.end(), 0);
//   assert(numOfTransportCycle == std::accumulate(recvProcList.begin(), recvProcList.end(), 0));

    // diffMの平均化
    const std::size_t myStartCount = averageNumOfElements * rank;
    const std::size_t myEndCount = std::min(averageNumOfElements * (rank +1), numOfElements);
//   std::cout << TlUtils::format("[%d] (%d ~ %d)", rank, myStartCount, myEndCount) << std::endl;

    std::size_t count = 0;
    //int c2 = 0;
    TlSparseSymmetricMatrix ans(globalSize);
    for (int i = 0; i < proc; ++i) {
        TlSparseSymmetricMatrix tmp(globalSize);
        if (i == rank) {
            tmp = m;

            // for debug
//       {
//  for (TlSparseSymmetricMatrix::const_iterator p = tmp.begin(); p != tmp.end(); ++p){
//    int row, col;
//    tmp.index(p->first, &row, &col);
//    std::cout << TlUtils::format("<<<<[%d] count=%d (%d, %d)", rank, c2, row, col) << std::endl;
//    c2++;
//  }
//       }

        }
        rComm.broadcast(tmp, i);

        TlSparseSymmetricMatrix::const_iterator pEnd = tmp.end();
        for (TlSparseSymmetricMatrix::const_iterator p = tmp.begin(); p != pEnd; ++p) {
            if ((myStartCount <= count) && (count < myEndCount)) {
                ans.set(*p);

//  int row, col;
//  ans.index(p->first, &row, &col);
//  std::cout << TlUtils::format(">>>>[%d] count=%d (%d, %d)", rank, count, row, col) << std::endl;
            }
            ++count;
        }

        rComm.barrier();
    }
    assert(count == numOfElements);

    // 最終成果物
    //m.merge(addM);

    // debug
//   {
//     rComm.barrier();
//     for (int i = 0; i < rComm.getNumOfProc(); ++i) {
//       if (i == rComm.getRank()) {
//  std::cout <<  i << " th matrix >>>>" << std::endl;
//  ans.print(std::cout);
//  std::cout << "<<<<\n" << std::endl;
//       }
//       rComm.barrier();
//     }
//   }

    return ans;
}


TlVector TlDistributeSymmetricMatrix::getPartialMatrix(int* pStartRow, int* pEndRow, int* pStartCol, int* pEndCol) const
{
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

    const int interval = (range + (proc -1)) / proc;
    const int localStart = interval * rank;
    const int localEnd = localStart + interval -1;

    TlVector P(interval);
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

                const int globalIndex = globalRowIndex + (2 * globalSize - globalColIndex)*(globalColIndex -1) / 2;
                if (localStart == globalIndex) {
                    *pStartRow = globalRowIndex;
                    *pStartCol = globalColIndex;
                }
                if (localEnd == globalIndex) {
                    *pEndRow = globalRowIndex;
                    *pEndCol = globalColIndex;
                }

                if ((localStart <= globalIndex) && (globalIndex <= localEnd)) {
                    const std::size_t index = r + numOfLocalRows * c; // row-major
                    const double value = pBuf[index];
                    P[(globalIndex - localStart)] = value;
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


void TlDistributeSymmetricMatrix::getPartialMatrix(TlSparseSymmetricMatrix& ioMatrix) const
{
    assert(this->getNumOfRows() == ioMatrix.getNumOfRows());
    assert(this->getNumOfCols() == ioMatrix.getNumOfCols());

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    for (int i = 0; i < nProc; ++i) {
        TlSparseSymmetricMatrix tmp;
        if (i == nRank) {
            tmp = ioMatrix;
        }
        rComm.broadcast(tmp, i);

        TlSparseSymmetricMatrix::iterator pEnd = tmp.end();
        for (TlSparseSymmetricMatrix::iterator p = tmp.begin(); p != pEnd; ++p) {
            //int nGlobalRow = p->first.row;
            //int nGlobalCol = p->first.col;
            unsigned long globalIndex = p->first;
            int nGlobalRow = 0;
            int nGlobalCol = 0;
            tmp.index(globalIndex, &nGlobalRow, &nGlobalCol);

            if (nGlobalRow < nGlobalCol) {
                std::swap(nGlobalRow, nGlobalCol);
            }
            const double dValue = this->get(nGlobalRow, nGlobalCol);
            p->second = dValue;
        }

        if (i == nRank) {
            ioMatrix = tmp;
        }
    }
}


// 複数スレッドから同時に呼び出さない
bool TlDistributeSymmetricMatrix::getSparseMatrixX(TlSparseSymmetricMatrix* pMatrix,
                                                    bool isFinalize) const
{
    return TlDistributeMatrix::getSparseMatrixX(pMatrix, isFinalize);
}


bool TlDistributeSymmetricMatrix::getPartialMatrixX(TlPartialSymmetricMatrix* pMatrix,
                                                    bool isFinalize) const
{
    return TlDistributeMatrix::getPartialMatrixX(pMatrix, isFinalize);
}


void TlDistributeSymmetricMatrix::getPartialMatrixX_registerTask(TlPartialMatrix* pMatrix) const
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

            if (startRow != startCol) {
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
            } else {
                for (index_type row = startRow; row < endRow; ++row) {
                    for (index_type col = startCol; col <= row; ++col) {
                        
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
            }
            
            
            this->gpmClientTasks_[pMatrix] = task;
        }
    }
}



void TlDistributeSymmetricMatrix::mergeSparseMatrix(const TlSparseSymmetricMatrix& M)
{
    assert(M.getNumOfRows() == this->getNumOfRows());
    assert(M.getNumOfCols() == this->getNumOfCols());
    
    // 送信すべきインデックスリストの作成
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    TlSparseMatrix::const_iterator itEnd = M.end();
    for (TlSparseMatrix::const_iterator it = M.begin(); it != itEnd; ++it) {
        TlSparseMatrix::KeyType key = it->first;
        index_type globalRow, globalCol;
        M.index(key, &globalRow, &globalCol);
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


void TlDistributeSymmetricMatrix::mergePartialMatrix(const TlPartialSymmetricMatrix& M)
{
    assert(M.getNumOfRows() == this->getNumOfRows());
    assert(M.getNumOfCols() == this->getNumOfCols());

    // 送信すべきインデックスリストの作成
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    const index_type startRow = M.getStartRow();
    const index_type startCol = M.getStartCol();
    const index_type endRow = startRow + M.getRowRange();
    const index_type endCol = startCol + M.getColRange();

    for (index_type globalRow = startRow; globalRow < endRow; ++globalRow) {
        for (index_type globalCol = startCol; globalCol < endCol; ++globalCol) {
            if (globalRow < globalCol) {
                continue;
            }

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
}


void TlDistributeSymmetricMatrix::mergePartialMatrix(const std::list<TlPartialSymmetricMatrix>& matrixArray)
{
    // 送信すべきインデックスリストの作成
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    std::vector<std::vector<index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    std::list<TlPartialSymmetricMatrix>::const_iterator itEnd = matrixArray.end();
    for (std::list<TlPartialSymmetricMatrix>::const_iterator it = matrixArray.begin();
         it != itEnd; ++it) {
        const TlPartialSymmetricMatrix& M = *it;
        assert(M.getNumOfRows() == this->getNumOfRows());
        assert(M.getNumOfCols() == this->getNumOfCols());

        const index_type startRow = M.getStartRow();
        const index_type startCol = M.getStartCol();
        const index_type endRow = startRow + M.getRowRange();
        const index_type endCol = startCol + M.getColRange();
        
        for (index_type globalRow = startRow; globalRow < endRow; ++globalRow) {
            for (index_type globalCol = startCol; globalCol < endCol; ++globalCol) {
                if (globalRow < globalCol) {
                    continue;
                }
                
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
    }
    

    this->mergeMatrix_common(indexArrays, values);
}


void TlDistributeSymmetricMatrix::mergeSparseMatrixAsync(const TlSparseMatrix* pMatrix, bool isFinalize)
{
    if (pMatrix != NULL) {
        assert(pMatrix->getNumOfRows() == this->getNumOfRows());
        assert(pMatrix->getNumOfCols() == this->getNumOfCols());
        
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int numOfProcs = rComm.getNumOfProc();
        
        // 送信すべきインデックスリストの作成
        std::vector<std::vector<index_type> > indexArrays(numOfProcs);
        std::vector<std::vector<double> > values(numOfProcs);
        
        TlSparseSymmetricMatrix::const_iterator itEnd = pMatrix->end();
        for (TlSparseSymmetricMatrix::const_iterator it = pMatrix->begin(); it != itEnd; ++it) {
            TlSparseSymmetricMatrix::KeyType key = it->first;
            index_type globalRow, globalCol;
            pMatrix->index(key, &globalRow, &globalCol);
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


double TlDistributeSymmetricMatrix::getMaxAbsoluteElement(int* pOutRow, int* pOutCol) const
{
    return TlDistributeMatrix::getMaxAbsoluteElement(pOutRow, pOutCol);
}


// const TlDistributeSymmetricMatrix& TlDistributeSymmetricMatrix::dot(const TlDistributeSymmetricMatrix& X)
// {
//     // 親メンバ関数と同じ
//     assert(this->getNumOfRows() == X.getNumOfRows());
//     assert(this->getNumOfCols() == X.getNumOfCols());

//     const size_type size = this->getNumOfMyElements();

// #ifdef _OPENMP
//     // use OpenMP
//     const size_type quot = size / MAX_LOOP;
//     const int rem = size - quot * MAX_LOOP;
// #pragma omp parallel
//     {
//         for (size_type block = 0; block < quot; ++block) {
//             const size_type index_base = block * MAX_LOOP;
// #pragma omp for
//             for (int i = 0; i < MAX_LOOP; ++i) {
//                 const size_type index = index_base + i;
//                 this->pData_[index] *= X.pData_[index];
//             }
//         }

//         const size_type index_base = quot * MAX_LOOP;
// #pragma omp for
//         for (int i = 0; i < rem; ++i) {
//             const size_type index = index_base + i;
//             this->pData_[index] *= X.pData_[index];
//         }
//     }
// #else
//     // not use OpenMP
//     for (size_type index = 0; index < size; ++index) {
//         this->data_[index] *= X.data_[index];
//     }
// #endif // _OPENMP
    
//     return (*this);
// }


double TlDistributeSymmetricMatrix::sum() const
{
    double answer = 0.0;

    // 対角項と下半分非対角項の和
    const int numOfLocalRows = this->m_nMyRows;
    const int numOfLocalCols = this->m_nMyCols;
    for (int localRowIndex = 0; localRowIndex < numOfLocalRows; ++localRowIndex) {
        const index_type row = this->m_RowIndexTable[localRowIndex];
        for (int localColIndex = 0; localColIndex < numOfLocalCols; ++localColIndex) {
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


bool TlDistributeSymmetricMatrix::isLoadable(std::ifstream& ifs)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool answer = false;

    if (rComm.isMaster() == true) {
        answer = TlSymmetricMatrix::isLoadable(ifs);
    }
    rComm.broadcast(answer);
    return answer;
}

bool TlDistributeSymmetricMatrix::isLoadable(const std::string& rFilePath)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool answer = false;

    if (rComm.isMaster() == true) {
        answer = TlSymmetricMatrix::isLoadable(rFilePath);
    }
    rComm.broadcast(answer);
    return answer;
}

bool TlDistributeSymmetricMatrix::load(const std::string& sFilePath)
{
    if (TlDistributeMatrix::isUsingPartialIO == true) {
        return this->loadLocal(sFilePath);
    }
    
    TlCommunicate& rComm = TlCommunicate::getInstance();

    std::ifstream ifs;
    if (rComm.isMaster() == true) {
        ifs.open(sFilePath.c_str(), std::ifstream::binary | std::ifstream::in);
    }

    bool bIsFail = false;
    if (rComm.isMaster() == true) {
        bIsFail = ifs.fail();
    }
    rComm.broadcast(bIsFail);
    if (bIsFail) {
#ifdef DEBUG
        std::cerr << "[error] TlDistributeSymmetricMatrix::load(): could not open file. " << sFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    bool bAnswer = this->load(ifs);

    if (bAnswer != true) {
        if (rComm.isMaster() == true) {
            std::cerr << "TlDistributeSymmetricMatrix::load() is not supported: " << sFilePath << std::endl;
            std::abort();
        }
    }

    if (rComm.isMaster() == true) {
        ifs.close();
    }

    rComm.broadcast(bAnswer);
    return bAnswer;
}


bool TlDistributeSymmetricMatrix::load(std::ifstream& ifs)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();

    bool bAnswer = true;
    enum {RSFD, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};

    // read header
    int nType = 0;
    if (rComm.isMaster() == true) {
        ifs.read((char*)&nType,           sizeof(int));
        ifs.read((char*)&(this->m_nRows), sizeof(int));
        ifs.read((char*)&(this->m_nCols), sizeof(int));
    }
    rComm.broadcast(bAnswer);
    rComm.broadcast(nType);

    switch (nType) {
    case RLHD:
        //std::cerr << "load RLHD" << std::endl;
        break;
    default:
        if (rComm.isMaster() == true) {
            std::cerr << "this matrix type is not supported. stop." << std::endl;
            std::cerr << "type=" << nType << ", row=" << this->m_nRows << ", col=" << this->m_nCols << std::endl;
        }
        bAnswer = false;
        break;
    }

    rComm.broadcast(this->m_nRows);
    rComm.broadcast(this->m_nCols);
    this->initialize();
    //std::cerr << "TlDistributeSymmetricMatrix::load() init." << std::endl;

    if (rComm.isMaster() == true) {
        static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
        std::vector<double> buf(bufferCount, 0.0);

        index_type currentRow = 0;
        index_type currentCol = 0;
        const index_type maxRow = this->m_nRows;
        //const index_type maxCol = this->m_nCols;
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
                if (currentCol > currentRow) {
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
                assert(numOfContents == tmpValueLists[proc].size());

                if (numOfContents == 0) {
                    // 送るリストがない場合は送らない
                    continue;
                }
                
                if (isSendData[proc] == true) {
                    rComm.wait(sizeLists[proc]);
                    rComm.wait(&(rowColLists[proc][0]));
                    rComm.wait(&(valueLists[proc][0]));
                    isSendData[proc] = false;
                }
                
                sizeLists[proc] = numOfContents;
                rowColLists[proc] = tmpRowColLists[proc];
                valueLists[proc] = tmpValueLists[proc];
                
                rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
                rComm.iSendDataX(&(rowColLists[proc][0]), (sizeLists[proc] * 2), proc, TAG_LOAD_ROWCOLS);
                rComm.iSendDataX(&(valueLists[proc][0]), sizeLists[proc], proc, TAG_LOAD_VALUES);
                isSendData[proc] = true;
            }
        } // end while

        for (int proc = 1; proc < numOfProcs; ++proc) { // proc == 0 は送信しない
            if (isSendData[proc] == true) {
                rComm.wait(sizeLists[proc]);
                rComm.wait(&(rowColLists[proc][0]));
                rComm.wait(&(valueLists[proc][0]));
                sizeLists[proc] = false;
            }
        }
        
        // 終了メッセージを全ノードに送る
        for (int proc = 1; proc < numOfProcs; ++proc) { // proc == 0 は送信しない
            if (isSendData[proc] == true) {
                rComm.wait(sizeLists[proc]);
                rComm.wait(&(rowColLists[proc][0]));
                rComm.wait(&(valueLists[proc][0]));
            }
        }
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
        bool isBreakLoop = false;
        while (isBreakLoop == false) {
            if (rComm.test(sizeList) == true) {
                rComm.wait(sizeList);
                // std::cerr << TlUtils::format("RECV [%d] sizeList=",
                //                              rComm.getRank(), sizeList)
                //           << std::endl;
                assert(sizeList > 0);

                rowColList.resize(sizeList * 2);
                valueList.resize(sizeList);
                rComm.iReceiveDataX(&(rowColList[0]), sizeList * 2, root, TAG_LOAD_ROWCOLS);
                rComm.iReceiveDataX(&(valueList[0]), sizeList, root, TAG_LOAD_VALUES);
                rComm.wait(&(rowColList[0]));
                rComm.wait(&(valueList[0]));
                rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
                
                for (int i = 0; i < sizeList; ++i) {
                    const index_type row = rowColList[i * 2    ];
                    const index_type col = rowColList[i * 2 + 1];
                    const size_type index = this->getIndex(row, col);
                    assert(index != -1);
                    this->pData_[index] = valueList[i];
                }
            }

            if (rComm.test(endMsg) == true) {
                rComm.wait(endMsg);
                rComm.cancel(sizeList);
                // std::cerr << TlUtils::format("RECV [%d] END",
                //                              rComm.getRank())
                //           << std::endl;
                isBreakLoop = true;
            }
        }
    }

    // std::cerr << TlUtils::format("[%d] END", rComm.getRank())
    //           << std::endl;
    return bAnswer;
}


bool TlDistributeSymmetricMatrix::save(const std::string& sFilePath) const
{
    if (TlDistributeMatrix::isUsingPartialIO == true) {
        return this->saveLocal(sFilePath);
    }
    
    bool bAnswer = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        // master
        assert(this->m_nRows == this->m_nCols);
        const int nGlobalSize = this->m_nRows;
        TlFileSymmetricMatrix fm(sFilePath, nGlobalSize);

        // store local matrix
        {
            const int nRows = this->m_nMyRows;
            const int nCols = this->m_nMyCols;
            for (int r = 0; r < nRows; ++r) {
                const int nGlobalRowIndex = this->m_RowIndexTable[r];
                if (nGlobalRowIndex >= nGlobalSize) {
                    continue;
                }

                for (int c = 0; c < nCols; ++c) {
                    const int nGlobalColIndex = this->m_ColIndexTable[c];
                    if (nGlobalColIndex >= nGlobalSize) {
                        continue;
                    }
                    if (nGlobalRowIndex < nGlobalColIndex) {
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
            int msg = 0;
            int src = 0;
            rComm.receiveDataFromAnySource(msg, &src, TAG_SAVE_HANDSHAKE);
            assert(msg == 0);
            msg = 0;
            rComm.sendData(msg, src, TAG_SAVE_HANDSHAKE_OK);

            index_type nRows;
            index_type nCols;
            std::vector<int> rowIndexTable;
            std::vector<int> colIndexTable;
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
                if (nGlobalRowIndex >= nGlobalSize) {
                    continue;
                }

                for (int c = 0; c < nMaxCols; ++c) {
                    const int nGlobalColIndex = colIndexTable[c];
                    if (nGlobalColIndex >= nGlobalSize) {
                        continue;
                    }

                    // for symmetry
                    if (nGlobalColIndex > nGlobalRowIndex) {
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
    return bAnswer;
}


bool TlDistributeSymmetricMatrix::saveLocal(const std::string& filePath) const
{
    //std::cerr << "TlDistributeSymmetricMatrix::saveLocal() file=" << filePath << std::endl;
    
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //DataType& data_tmp = const_cast<DataType&>(this->data_);

    bool answer = true;
    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

    const int type = 18; // means Distributed(16) + RLHD(2)
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


bool TlDistributeSymmetricMatrix::loadLocal(const std::string& filePath)
{
    //std::cerr << "TlDistributeSymmetricMatrix::loadLocal() file=" << filePath << std::endl;
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


bool TlDistributeSymmetricMatrix::diagonal(TlVector* pEigVal, TlDistributeMatrix* pEigVec,
                                           TlDistributeSymmetricMatrix::DIAGONAL_METHOD method)
{
#ifdef HAVE_SCALAPACK
    if (method == DIVIDE_AND_CONQUER) {
        return diagonalByScaLapack_DC(*this, pEigVal, pEigVec);
    } else {
        return diagonalByScaLapack_QR(*this, pEigVal, pEigVec);
    }
#else
    std::cerr << "sorry. this code is not implemented." << std::endl;
    abort();
#endif // HAVE_SCALAPACK
}

bool TlDistributeSymmetricMatrix::inverse()
{
#ifdef HAVE_SCALAPACK
    // using SCALAPACK
    return inverseByScaLapack(*this);
#else
    // without SCALAPACK
#error NOT found algebra package: need SCALAPACK library
#endif // HAVE_SCALAPACK  
}

TlDistributeMatrix multiplicationByScaLapack(const TlDistributeSymmetricMatrix& X,
                                             const TlDistributeMatrix& Y)
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //TlDistributeSymmetricMatrix& Xtmp = const_cast<TlDistributeSymmetricMatrix&>(X);
    //TlDistributeMatrix& Ytmp = const_cast<TlDistributeMatrix&>(Y);

    assert(X.getNumOfCols() == Y.getNumOfRows());

    //X.symmetrize();
    TlDistributeMatrix Z(X.m_nRows, Y.m_nCols);

    // use SCALAPACK
    const char SIDE = 'L';                  // L means "C := alpha*A*B + beta*C",
    // R means "C := alpha*B*A + beta*C"
    const char UPLO = 'L';                  // L means the lower triangular part of the symmetric matrix
    // U means the upper triangular part of the symmetric matrix
    const int M = Z.getNumOfRows();        // the number of rows of the matrix  C
    const int N = Z.getNumOfCols();        // the number of columns of the matrix C
    const double ALPHA = 1.0;               // ALPHA specifies the scalar alpha

    const double* A = X.pData_;

    const int IA = 1;
    const int JA = 1;

    //const double* B = &(const_cast<TlDistributeMatrix&>(Y).m_aMatrix[0]);    // DIMENSION (LDB, n)
    const double* B = Y.pData_;
    const int IB = 1;
    const int JB = 1;

    const double BETA = 0.0;                // BETA  specifies the scalar  beta
    double* C = Z.pData_;          // DIMENSION (LDC, n)
    const int IC = 1;
    const int JC = 1;

    pdsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &IA, &JA, X.m_pDESC,
            B, &IB, &JB, Y.m_pDESC, &BETA,
            C, &IC, &JC, Z.m_pDESC);

    return Z;
}

TlDistributeMatrix multiplicationByScaLapack(const TlDistributeMatrix& X, const TlDistributeSymmetricMatrix& Y)
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    //TlDistributeMatrix& Xtmp = const_cast<TlDistributeMatrix&>(X);
    //TlDistributeSymmetricMatrix& Ytmp = const_cast<TlDistributeSymmetricMatrix&>(Y);

    assert(X.getNumOfCols() == Y.getNumOfRows());

    //Y.symmetrize();
    TlDistributeMatrix Z(X.m_nRows, Y.m_nCols);

    // use SCALAPACK
    const char SIDE = 'R';                  // L means "C := alpha*A*B + beta*C",
    // R means "C := alpha*B*A + beta*C"
    const char UPLO = 'L';                  // L means the lower triangular part of the symmetric matrix
    // U means the upper triangular part of the symmetric matrix
    const int M = Z.getNumOfRows();        // the number of rows of the matrix  C
    const int N = Z.getNumOfCols();        // the number of columns of the matrix C
    const double ALPHA = 1.0;               // ALPHA specifies the scalar alpha

    //const double* A = &(const_cast<TlDistributeMatrix&>(X).m_aMatrix[0]);    // DIMENSION (LDA, ka)
    const double* A = Y.pData_;
    const int IA = 1;
    const int JA = 1;

    //const double* B = &(const_cast<TlDistributeSymmetricMatrix&>(Y).m_aMatrix[0]);    // DIMENSION (LDB, n)
    const double* B = X.pData_;    // DIMENSION (LDB, n)
    const int IB = 1;
    const int JB = 1;

    const double BETA = 0.0;                // BETA  specifies the scalar  beta
    double* C = Z.pData_;          // DIMENSION (LDC, n)
    const int IC = 1;
    const int JC = 1;

    pdsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &IA, &JA, Y.m_pDESC,
            B, &IB, &JB, X.m_pDESC, &BETA,
            C, &IC, &JC, Z.m_pDESC);

    return Z;
}


bool diagonalByScaLapack_QR(TlDistributeSymmetricMatrix& inMatrix,
                            TlVector* outEigVal, TlDistributeMatrix* outEigVec)
{
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
    double* W = outEigVal->data_;

    outEigVec->resize(N, N);
    double* Z = outEigVec->pData_;

    const int IZ = 1;
    const int JZ = 1;
    const int* DESCZ = outEigVec->m_pDESC;

    int LWORK = -1;
    double* pWorkSizeCheck = new double[1];
    int info = 0;

    pdsyev_(JOBZ, UPLO, &N, pBufA,
            &IA, &JA, DESCA, W,
            Z,  &IZ,  &JZ, DESCZ, pWorkSizeCheck, &LWORK, &info);

    if (info != 0) {
        std::cout << "pdsyev_ error @1st call: " << info << std::endl;
        return false;
    }

    LWORK = (int)pWorkSizeCheck[0];
    double* pWork = new double[LWORK];

    delete[] pWorkSizeCheck;
    pWorkSizeCheck = NULL;

    pdsyev_(JOBZ, UPLO, &N, pBufA,
            &IA, &JA, DESCA, W,
            Z,  &IZ,  &JZ, DESCZ, pWork, &LWORK, &info);

    delete[] pWork;
    pWork = NULL;

    if (info != 0) {
        std::cout << "pdsyev_ error @2nd call: " << info << std::endl;
    }
    return ((info == 0) ? true : false);
}


// Divide-and-Conquer Algorithm
bool diagonalByScaLapack_DC(TlDistributeSymmetricMatrix& inMatrix,
                            TlVector* outEigVal, TlDistributeMatrix* outEigVec)
{
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
    double* W = outEigVal->data_;

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

    pdsyevd_(JOBZ, UPLO, &N, pBufA,
             &IA, &JA, DESCA, W,
             Z,  &IZ,  &JZ, DESCZ, pWorkSizeCheck, &LWORK,
             IWORK, &LIWORK, &info);

    if (info != 0) {
        std::cout << "pdsyev_ error @1st call: " << info << std::endl;
        return false;
    }

    LWORK = (int)pWorkSizeCheck[0];
    double* pWork = new double[LWORK];

    delete[] pWorkSizeCheck;
    pWorkSizeCheck = NULL;

    pdsyevd_(JOBZ, UPLO, &N, pBufA,
             &IA, &JA, DESCA, W,
             Z,  &IZ,  &JZ, DESCZ, pWork, &LWORK,
             IWORK, &LIWORK, &info);

    delete[] IWORK;
    IWORK = NULL;

    delete[] pWork;
    pWork = NULL;

    if (info != 0) {
        std::cout << "pdsyev_ error @2nd call: " << info << std::endl;
    }
    return ((info == 0) ? true : false);
}


bool inverseByScaLapack(TlDistributeSymmetricMatrix& X)
{
    TlDistributeMatrix Y(X);
    bool bAnswer = Y.inverse();
    if (bAnswer == true) {
        X = Y;
    } else {
        std::cout << "inverseByScaLapack return false." << std::endl;
        abort();
    }

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

//   return true;

    return bAnswer;
}

// void TlDistributeSymmetricMatrix::symmetrize() const {
//   const int nDim = this->m_nRows;
//   for (int row = 0; row < nDim; ++row){
//     const int nGlobalRow = row +1;
//     for (int col = 0; col < row; ++col){
//       assert(row > col);
//       const double dValue = this->get(row, col);

//       // copy (row, col) to (col, row)
//       const int nGlobalCol = col +1;
//       pdelset_(&(this->data_[0]), &nGlobalCol, &nGlobalRow, this->m_pDESC, &dValue);
//     }
//   }
// }

TlDistributeMatrix operator*(const TlDistributeSymmetricMatrix& X, const TlDistributeSymmetricMatrix& Y)
{
    TlDistributeMatrix Z = Y;
    return multiplicationByScaLapack(X, Z);
}

TlDistributeMatrix operator*(const TlDistributeSymmetricMatrix& X, const TlDistributeMatrix& Y)
{
    return multiplicationByScaLapack(X, Y);
}


TlDistributeMatrix operator*(const TlDistributeMatrix& X, const TlDistributeSymmetricMatrix& Y)
{
    return multiplicationByScaLapack(X, Y);
}


TlDistributeSymmetricMatrix operator*(const double X, const TlDistributeSymmetricMatrix& Y)
{
    TlDistributeSymmetricMatrix ans = Y;
    ans *= X;

    return ans;
}


TlDistributeVector operator*(const TlDistributeSymmetricMatrix& A, const TlDistributeVector& X)
{
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    TlDistributeSymmetricMatrix& Atmp = const_cast<TlDistributeSymmetricMatrix&>(A);
    TlDistributeVector& Xtmp = const_cast<TlDistributeVector&>(X);

    const int N = A.getNumOfCols();
    assert(N == X.getSize());

    TlDistributeVector Y(N);

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

    pdsymv_(&UPLO, &N,
            &alpha, Atmp.pData_, &IA, &JA, A.m_pDESC,
            &(Xtmp.data_[0]), &IX, &JX, X.m_pDESC, &INCX,
            &beta, &(Y.data_[0]), &IY, &JY, Y.m_pDESC, &INCY);

    return Y;
}


TlDistributeSymmetricMatrix operator*(const TlDistributeSymmetricMatrix& X, const double Y)
{
    TlDistributeSymmetricMatrix ans = X;
    ans *= Y;

    return ans;
}


// Harbrecht, Peter, Schneider, 2011
TlDistributeMatrix 
TlDistributeSymmetricMatrix::choleskyFactorization(const double threshold) const
{
    TlLogging& log = TlLogging::getInstance();

    const index_type N = this->getNumOfRows();
    TlVector d = this->getDiagonalElements();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    TlDistributeMatrix L(N, N);
    index_type m = 0;
    double sum_ll = 0.0;

    while (error > threshold) {
        std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                       pivot.end());
        const index_type i = it - pivot.begin();
        std::swap(pivot[m], pivot[i]);
        
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(m, pivot[m], l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;
        
        for (index_type i = m +1; i < N; ++i) {
            double sum_ll = 0.0;
            for (index_type j = 0; j < m; ++j) {
                sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
            }
            const double l_m_pi = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
            L.set(m, pivot[i], l_m_pi);
            
            d[pivot[i]] -= l_m_pi * l_m_pi;
        }
            
        error = 0.0;
        for (index_type i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }

        log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
        ++m;
    }

    L.transpose();
    L.resize(N, m);

    return L;
}


// TlDistributeMatrix 
// TlDistributeSymmetricMatrix::choleskyFactorization_mod0(const double threshold) const
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     TlLogging& log = TlLogging::getInstance();

//     const index_type N = this->getNumOfRows();
//     TlVector d = this->getDiagonalElements();
//     double error = d.sum();
//     std::vector<TlVector::size_type> pivot(N);
//     for (index_type i = 0; i < N; ++i) {
//         pivot[i] = i;
//     }

//     TlDistributeMatrix L(N, N);
//     index_type m = 0;
//     double sum_ll = 0.0;

//     while (error > threshold) {
//         std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
//                                                                        pivot.end());
//         const index_type i = it - pivot.begin();
//         std::swap(pivot[m], pivot[i]);
        
//         const double l_m_pm = std::sqrt(d[pivot[m]]);
//         L.set(m, pivot[m], l_m_pm);
        
//         const double inv_l_m_pm = 1.0 / l_m_pm;
        
//         // L(0:m, pivot[m])を一時保管
//         const TlVector L_x_pm = L.getColVector(pivot[m]); // TODO: バッファの摂り過ぎ

//         for (index_type i = m +1; i < N; ++i) {
//             double sum_ll = 0.0;
//             for (index_type j = 0; j < m; ++j) {
//                 //sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
//                 sum_ll += L_x_pm[j] * L.getLocal(j, pivot[i]);
//             }
//             rComm.allReduce_SUM(sum_ll);
//             const double l_m_pi = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
//             L.set(m, pivot[i], l_m_pi);
            
//             d[pivot[i]] -= l_m_pi * l_m_pi;
//         }
            
//         error = 0.0;
//         for (index_type i = m +1; i < N; ++i) {
//             error += d[pivot[i]];
//         }

//         log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
//         ++m;
//     }

//     L.transpose();
//     L.resize(N, m);

//     return L;
// }


TlDistributeMatrix 
TlDistributeSymmetricMatrix::choleskyFactorization_mod(const double threshold) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlLogging& log = TlLogging::getInstance();

    const index_type N = this->getNumOfRows();
    TlVector d = this->getDiagonalElements();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    TlDistributeMatrix L(N, N);
    index_type m = 0;
    double sum_ll = 0.0;

    while (error > threshold) {
        std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                       pivot.end());
        const index_type i = it - pivot.begin();
        std::swap(pivot[m], pivot[i]);
        
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(m, pivot[m], l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;
        
        // L(0:m, pivot[m])を一時保管
        const TlVector L_x_pm = L.getColVector(pivot[m]); // TODO: バッファの取り過ぎ

        for (index_type i = m +1; i < N; ++i) {
            double sum_ll = 0.0;
            for (index_type j = 0; j < m; ++j) {
                //sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
                sum_ll += L_x_pm[j] * L.getLocal(j, pivot[i]);
            }
            rComm.allReduce_SUM(sum_ll);

            double l_m_pi = 0.0;
            {
                // double G_pm_pi_exact = this->get(pivot[m], pivot[i]);
                // double l_m_pi_exact = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
                // double G_pm_pi_local = this->getLocal(pivot[m], pivot[i]);

                // for (int p = 0; p < rComm.getNumOfProcs(); ++p) {
                //     if (p == rComm.getRank()) {
                //         std::cerr << TlUtils::format("[%d] G(%d, %d)=% f(% f), holder=%d, index=%d",
                //                                      rComm.getRank(),
                //                                      pivot[m], pivot[i],
                //                                      G_pm_pi_exact, G_pm_pi_local,
                //                                      this->getProcIdForIndex(pivot[m], pivot[i]),
                //                                      this->getIndex(pivot[m], pivot[i]))
                //                   << std::endl;
                //     }
                //     rComm.barrier();
                // }

                int holder = this->getProcIdForIndex(pivot[m], pivot[i]);
                if (holder == rComm.getRank()) {
                    assert(this->getIndex(pivot[m], pivot[i]) != -1);
                    //assert(std::fabs(this->getLocal(pivot[m], pivot[i]) - G_pm_pi_exact) < 1.0E-5);
                    l_m_pi = (this->getLocal(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
                }
                rComm.broadcast(l_m_pi, holder);
            }
            //assert(std::fabs(l_m_pi - l_m_pi_exact) < 1.0E-5);
            L.set(m, pivot[i], l_m_pi);
            
            d[pivot[i]] -= l_m_pi * l_m_pi;
        }
            
        error = 0.0;
        for (index_type i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }

        log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
        ++m;
    }

    L.transpose();
    L.resize(N, m);

    return L;
}


// TlDistributeMatrix 
// TlDistributeSymmetricMatrix::choleskyFactorization_mod(const double threshold) const
// {
//     TlLogging& log = TlLogging::getInstance();
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     const index_type N = this->getNumOfRows();
//     TlVector d = this->getDiagonalElements();
//     double error = d.sum();
//     std::vector<TlVector::size_type> pivot(N);
//     for (index_type i = 0; i < N; ++i) {
//         pivot[i] = i;
//     }

//     TlDistributeMatrix L(N, N);
//     index_type m = 0;
//     double sum_ll = 0.0;

//     while (error > threshold) {
//         std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
//                                                                        pivot.end());
//         const index_type i = it - pivot.begin();
//         std::swap(pivot[m], pivot[i]);
        
//         const double l_m_pm = std::sqrt(d[pivot[m]]);
//         L.set(m, pivot[m], l_m_pm);
        
//         const double inv_l_m_pm = 1.0 / l_m_pm;

//         // L(0:m, pivot[m])を一時保管
//         const TlVector L_x_pm = L.getColVector(pivot[m]); // TODO: バッファの摂り過ぎ

//         for (index_type i = m +1; i < N; ++i) {
//             double sum_ll = 0.0;
//             for (index_type j = 0; j < m; ++j) {
//                 //sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
//                 sum_ll += L_x_pm[j] * L.getLocal(j, pivot[i]);
//             }
//             rComm.allReduce_SUM(sum_ll);

//             double l_m_pi = 0.0;
//             l_m_pi = (this->getLocal(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
//             rComm.allReduce_SUM(l_m_pi);
//             L.set(m, pivot[i], l_m_pi);
            
//             d[pivot[i]] -= l_m_pi * l_m_pi;
//         }
            
//         error = 0.0;
//         for (index_type i = m +1; i < N; ++i) {
//             error += d[pivot[i]];
//         }

//         log.info(TlUtils::format("cholesky: m=%d, err=%e", m, error));
//         ++m;
//     }

//     L.transpose();
//     L.resize(N, m);

//     return L;
// }
