#include "tl_dense_symmetric_matrix_impl_scalapack.h"
#include "TlCommunicate.h"
#include "scalapack.h"
#include "tl_dense_symmetric_matrix_io.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_matrix_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ImplScalapack::TlDenseSymmetricMatrix_ImplScalapack(
    const TlMatrixObject::index_type dim)
    : TlDenseGeneralMatrix_ImplScalapack(dim, dim) {}

TlDenseSymmetricMatrix_ImplScalapack::TlDenseSymmetricMatrix_ImplScalapack(
    const TlDenseSymmetricMatrix_ImplScalapack& rhs)
    : TlDenseGeneralMatrix_ImplScalapack(rhs) {}

// dim = rows
TlDenseSymmetricMatrix_ImplScalapack::TlDenseSymmetricMatrix_ImplScalapack(
    const TlDenseGeneralMatrix_ImplScalapack& rhs)
    : TlDenseGeneralMatrix_ImplScalapack(rhs) {
  if (rhs.getNumOfRows() != rhs.getNumOfCols()) {
    this->log_.critical(TlUtils::format("dims are not consistent: %d != %d",
                                        rhs.getNumOfRows(),
                                        rhs.getNumOfCols()));
    throw;
  }
  assert(rhs.getNumOfRows() == rhs.getNumOfCols());
  // コピーされたバッファの下半分しか使わない
}

// TlDenseSymmetricMatrix_ImplScalapack::TlDenseSymmetricMatrix_ImplScalapack(
//     const TlDenseVector_ImplScalapack& v, const TlMatrixObject::index_type
//     dim) : TlDenseGeneralMatrix_ImplScalapack(dim, dim) {
//   this->restore(v);
// }

TlDenseSymmetricMatrix_ImplScalapack::~TlDenseSymmetricMatrix_ImplScalapack() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
double TlDenseSymmetricMatrix_ImplScalapack::get(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
  // if (row < col) {
  //   std::swap(row, col);
  // }

  return TlDenseGeneralMatrix_ImplScalapack::get(row, col);
}

void TlDenseSymmetricMatrix_ImplScalapack::set(TlMatrixObject::index_type row,
                                               TlMatrixObject::index_type col,
                                               const double value) {
  TlDenseGeneralMatrix_ImplScalapack::set(row, col, value);
  TlDenseGeneralMatrix_ImplScalapack::set(col, row, value);
}

void TlDenseSymmetricMatrix_ImplScalapack::add(TlMatrixObject::index_type row,
                                               TlMatrixObject::index_type col,
                                               const double value) {
  TlDenseGeneralMatrix_ImplScalapack::add(row, col, value);
  if (row != col) {
    TlDenseGeneralMatrix_ImplScalapack::add(col, row, value);
  }
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ImplScalapack::eig(
    TlDenseVector_ImplLapack* pEigVal,
    TlDenseGeneralMatrix_ImplScalapack* pEigVec,
    TlDenseSymmetricMatrix_ImplScalapack::DIAGONAL_METHOD method) const {
  if (method == DIVIDE_AND_CONQUER) {
    return diagonalByScaLapack_DC(*this, pEigVal, pEigVec);
  } else {
    return diagonalByScaLapack_QR(*this, pEigVal, pEigVec);
  }
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrix_ImplScalapack::load(const std::string& sFilePath) {
  // if (TlDenseGeneralMatrix_ImplScalapack::isUsingPartialIO == true) {
  //   return this->loadLocal(sFilePath);
  // }

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
    this->log_.critical(TlUtils::format("Not supported file type: %s @%s.%d",
                                        sFilePath.c_str(), __FILE__, __LINE__));
    std::abort();
  }

  if (rComm.isMaster() == true) {
    fs.close();
  }

  rComm.broadcast(bAnswer);
  return bAnswer;
}

bool TlDenseSymmetricMatrix_ImplScalapack::load(std::fstream& fs) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());

  const int numOfProcs = rComm.getNumOfProc();

  // read header
  bool answer = true;
  TlMatrixObject::MatrixType matType;
  TlMatrixObject::index_type dim;
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

  this->rows_ = dim;
  this->cols_ = dim;
  assert(this->getNumOfRows() == this->getNumOfCols());
  this->initialize();

  if (rComm.isMaster() == true) {
    static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
    std::vector<double> buf(bufferCount, 0.0);

    TlMatrixObject::index_type row = 0;
    TlMatrixObject::index_type col = 0;
    TlMatrixObject::size_type count = 0;
    const TlMatrixObject::size_type maxCount = dim * (dim + 1) / 2;
    bool isFinished = false;

    std::vector<TlMatrixObject::index_type> sizeLists(numOfProcs, 0);
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
        int rank = 0;

        this->getGlobalRowCol2LocalRowCol(row, col, &rank);
        elements[rank].push_back(
            TlMatrixObject::MatrixElement(row, col, buf[i]));

        this->getGlobalRowCol2LocalRowCol(col, row, &rank);
        elements[rank].push_back(
            TlMatrixObject::MatrixElement(col, row, buf[i]));

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
        const TlMatrixObject::index_type sizeList = elements[0].size();
        for (TlMatrixObject::index_type i = 0; i < sizeList; ++i) {
          const TlMatrixObject::index_type globalRow = elements[0][i].row;
          const TlMatrixObject::index_type globalCol = elements[0][i].col;

          int rank = 0;
          TlMatrixObject::index_type myRow, myCol;
          this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow,
                                            &myCol);
          assert(rank == this->rank_);
          const TlMatrixObject::size_type index =
              this->getLocalIndex(myRow, myCol);
          assert((0 <= index) && (index < this->getNumOfMyElements()));
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
    // sizeList=-1 で終了
    std::vector<int> endMsg(numOfProcs, -1);
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.iSendData(endMsg[proc], proc, TAG_LOAD_SIZE);
    }
    for (int proc = 1; proc < numOfProcs; ++proc) {  // proc == 0 は送信しない
      rComm.wait(endMsg[proc]);
    }
  } else {
    // slave
    const int root = 0;
    TlMatrixObject::index_type sizeList = 0;
    std::vector<TlMatrixObject::MatrixElement> elements;
    int endMsg = 0;

    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
    bool isBreakLoop = false;
    while (isBreakLoop == false) {
      if (rComm.test(sizeList) == true) {
        rComm.wait(sizeList);
        if (sizeList >= 0) {
          if (sizeList > 0) {
            elements.resize(sizeList);
            rComm.receiveDataX(&(elements[0]), sizeList, root, TAG_LOAD_VALUES);

            for (TlMatrixObject::index_type i = 0; i < sizeList; ++i) {
              const TlMatrixObject::index_type globalRow = elements[i].row;
              const TlMatrixObject::index_type globalCol = elements[i].col;

              int rank = 0;
              TlMatrixObject::index_type myRow, myCol;
              this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank,
                                                &myRow, &myCol);
              assert(rank == this->rank_);
              const TlMatrixObject::size_type index =
                  this->getLocalIndex(myRow, myCol);
              this->pData_[index] = elements[i].value;
            }
          }
          rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
        } else {
          assert(sizeList < 0);
          isBreakLoop = true;
        }
      }
    }
  }

  assert(rComm.checkNonBlockingCommunications());
  return answer;
}

bool TlDenseSymmetricMatrix_ImplScalapack::save(
    const std::string& sFilePath) const {
  // if (TlDenseGeneralMatrix_ImplScalapack::isUsingPartialIO == true) {
  //   return this->saveLocal(sFilePath);
  // }

  bool answer = true;

  TlCommunicate& rComm = TlCommunicate::getInstance();
  assert(rComm.checkNonBlockingCommunications());

  if (rComm.isMaster() == true) {
    // master
    const int globalSize = this->getNumOfRows();
    TlFileSymmetricMatrix fm(sFilePath, globalSize);

    // store local matrix
    {
      const std::vector<TlMatrixObject::MatrixElement> buf =
          this->getMatrixElementsInLocal2();
      TlDenseGeneralMatrix_ImplScalapack::saveElements(&fm, buf);
    }

    // recive submatrix & write
    const int numOfSlaves = rComm.getNumOfProc() - 1;
    for (int i = 0; i < numOfSlaves; ++i) {
      int src = 0;
      int tag = TAG_SAVE_HANDSHAKE;
      TlMatrixObject::size_type bufSize = 0;
      rComm.iReceiveDataFromAnySource(bufSize, tag);
      rComm.wait(&bufSize, &src);

      std::vector<TlMatrixObject::MatrixElement> buf(bufSize);
      rComm.receiveDataX(&(buf[0]), bufSize, src, TAG_SAVE_DATA);

      TlDenseGeneralMatrix_ImplScalapack::saveElements(&fm, buf);
    }
  } else {
    // slave: send submatrix
    const int root = 0;
    const std::vector<TlMatrixObject::MatrixElement> buf =
        this->getMatrixElementsInLocal2();
    const TlMatrixObject::size_type bufSize = buf.size();
    rComm.sendData(bufSize, root, TAG_SAVE_HANDSHAKE);
    rComm.iSendDataX(&(buf[0]), buf.size(), root, TAG_SAVE_DATA);
    rComm.wait(&(buf[0]));
  }

  rComm.broadcast(answer);
  assert(rComm.checkNonBlockingCommunications());
  return answer;
}

std::vector<TlMatrixObject::MatrixElement>
TlDenseSymmetricMatrix_ImplScalapack::getMatrixElementsInLocal2() const {
  const TlMatrixObject::size_type numOfMyElements = this->getNumOfMyElements();
  std::vector<TlMatrixObject::MatrixElement> answer;
  answer.reserve(numOfMyElements);

  const TlMatrixObject::index_type localRows = this->myRows_;
  const TlMatrixObject::index_type localCols = this->myCols_;
  // TlMatrixObject::size_type count = 0;
  for (TlMatrixObject::index_type localRow = 0; localRow < localRows;
       ++localRow) {
    const TlMatrixObject::index_type globalRow =
        this->getLocal2Global_row(localRow);
    for (TlMatrixObject::index_type localCol = 0; localCol < localCols;
         ++localCol) {
      const TlMatrixObject::index_type globalCol =
          this->getLocal2Global_col(localCol);

      if (globalRow >= globalCol) {
        answer.push_back(TlMatrixObject::MatrixElement(
            globalRow, globalCol, this->getLocal(globalRow, globalCol)));
      }
    }
  }

  return answer;
}

// bool TlDenseSymmetricMatrix_ImplScalapack::saveLocal(
//     const std::string& filePath) const {
//   // std::cerr << "TlDenseSymmetricMatrix_blacs::saveLocal() file=" <<
//   filePath
//   // << std::endl;
//
//   // this const_cast is requiered for PGI compiler
//   // "error: expression must be an lvalue or a function designator"
//   // DataType& data_tmp = const_cast<DataType&>(this->data_);
//
//   bool answer = true;
//   std::ofstream ofs;
//   ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);
//
//   const int type = 18;  // means Distributed(16) + RLHD(2)
//   const index_type globalRow = this->m_nRows;
//   const index_type globalCol = this->m_nCols;
//   const index_type myRows = this->m_nMyRows;
//   const index_type myCols = this->m_nCols;
//   ofs.write(reinterpret_cast<const char*>(&type), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&globalRow), sizeof(index_type));
//   ofs.write(reinterpret_cast<const char*>(&globalCol), sizeof(index_type));
//   ofs.write(reinterpret_cast<const char*>(&myRows), sizeof(index_type));
//   ofs.write(reinterpret_cast<const char*>(&myCols), sizeof(index_type));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nRank), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nProc), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nProcGridRow),
//   sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nProcGridCol),
//   sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nMyProcRow), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nMyProcCol), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(&this->m_nBlockSize), sizeof(int));
//   ofs.write(reinterpret_cast<const char*>(this->pData_),
//             sizeof(double) * this->getNumOfMyElements());
//
//   ofs.close();
//
//   // std::cerr << "save: rank=" << this->m_nRank << std::endl;
//   // std::cerr << "save: proc=" << this->m_nProc << std::endl;
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   rComm.barrier();
//
//   return answer;
// }

// bool TlDenseSymmetricMatrix_ImplScalapack::loadLocal(
//     const std::string& filePath) {
//   // std::cerr << "TlDenseSymmetricMatrix_blacs::loadLocal() file=" <<
//   filePath
//   // << std::endl;
//   bool answer = false;
//
//   std::ifstream ifs;
//   ifs.open(filePath.c_str());
//   if (ifs.fail() != true) {
//     int type = 0;
//     index_type globalRow = 0;
//     index_type globalCol = 0;
//     index_type myRows = 0;
//     index_type myCols = 0;
//     int rank = 0;
//     int proc = 0;
//     int procGridRow = 0;
//     int procGridCol = 0;
//     int myProcRow = 0;
//     int myProcCol = 0;
//     int blockSize = 0;
//     ifs.read((char*)&type, sizeof(int));
//     ifs.read((char*)&globalRow, sizeof(index_type));
//     ifs.read((char*)&globalCol, sizeof(index_type));
//     ifs.read((char*)&myRows, sizeof(index_type));
//     ifs.read((char*)&myCols, sizeof(index_type));
//     ifs.read((char*)&rank, sizeof(int));
//     ifs.read((char*)&proc, sizeof(int));
//     ifs.read((char*)&procGridRow, sizeof(int));
//     ifs.read((char*)&procGridCol, sizeof(int));
//     ifs.read((char*)&myProcRow, sizeof(int));
//     ifs.read((char*)&myProcCol, sizeof(int));
//     ifs.read((char*)&blockSize, sizeof(int));
//
//     if (type != 18) {
//       std::cerr << TlUtils::format("file type mismatch(%d). ", type) <<
//       __FILE__
//                 << "," << __LINE__ << std::endl;
//     }
//     this->m_nRows = globalRow;
//     this->m_nCols = globalCol;
//     // std::cerr << "rank = " << rank << std::endl;
//     // std::cerr << "proc = " << proc << std::endl;
//     assert(rank == this->m_nRank);
//     assert(proc == this->m_nProc);
//     assert(procGridRow == this->m_nProcGridRow);
//     assert(procGridCol == this->m_nProcGridCol);
//     assert(myProcRow == this->m_nMyProcRow);
//     assert(myProcCol == this->m_nMyProcCol);
//     assert(blockSize == this->m_nBlockSize);
//     this->initialize();
//
//     ifs.read((char*)this->pData_, sizeof(double) *
//     this->getNumOfMyElements());
//
//     answer = true;
//   }
//
//   TlCommunicate& rComm = TlCommunicate::getInstance();
//   rComm.barrier();
//
//   return answer;
// }

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
TlMatrixObject::size_type TlDenseSymmetricMatrix_ImplScalapack::index(
    TlMatrixObject::index_type row, TlMatrixObject::index_type col) const {
  // This class treats 'U' type matrix.
  // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; in Fortran
  // index = row + col * (col +1) /2; // in C/C++ (in row >= col)
  assert((0 <= row) && (row < this->getNumOfRows()));
  assert((0 <= col) && (col < this->getNumOfCols()));

  if (row < col) {
    std::swap(row, col);
  }
  const TlMatrixObject::size_type index = row * (row + 1) / 2 + col;
  assert(index < this->getNumOfElements());

  return index;
}

// ---------------------------------------------------------------------------
bool diagonalByScaLapack_QR(
    const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
    TlDenseVector_ImplLapack* outEigVal,
    TlDenseGeneralMatrix_ImplScalapack* outEigVec) {
  TlLogging& logger = TlLogging::getInstance();
  assert(outEigVal != NULL);
  assert(outEigVec != NULL);

  const char* JOBZ = "V";
  const char* UPLO = "L";

  assert(inMatrix.getNumOfRows() == inMatrix.getNumOfCols());
  const int N = inMatrix.getNumOfRows();

  double* pBufA = inMatrix.pData_;

  const int IA = 1;
  const int JA = 1;
  const int* DESCA = inMatrix.pDESC_;

  outEigVal->resize(N);
  double* W = outEigVal->data();

  outEigVec->resize(N, N);
  double* Z = outEigVec->pData_;

  const int IZ = 1;
  const int JZ = 1;
  const int* DESCZ = outEigVec->pDESC_;

  int LWORK = -1;
  double* pWorkSizeCheck = new double[1];
  int info = 0;

  pdsyev_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ,
          pWorkSizeCheck, &LWORK, &info);

  if (info != 0) {
    logger.critical(TlUtils::format("pdsyev_ error @1st call: %d", info));
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
    logger.critical(TlUtils::format("pdsyev_ error @2nd call: %d", info));
  }
  return ((info == 0) ? true : false);
}

// Divide-and-Conquer Algorithm
bool diagonalByScaLapack_DC(
    const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
    TlDenseVector_ImplLapack* outEigVal,
    TlDenseGeneralMatrix_ImplScalapack* outEigVec) {
  TlLogging& logger = TlLogging::getInstance();
  assert(outEigVal != NULL);
  assert(outEigVec != NULL);

  const char* JOBZ = "V";
  const char* UPLO = "L";

  assert(inMatrix.getNumOfRows() == inMatrix.getNumOfCols());
  const int N = inMatrix.getNumOfRows();

  double* pBufA = inMatrix.pData_;

  const int IA = 1;
  const int JA = 1;
  const int* DESCA = inMatrix.pDESC_;

  outEigVal->resize(N);
  double* W = outEigVal->data();

  outEigVec->resize(N, N);
  double* Z = outEigVec->pData_;

  const int IZ = 1;
  const int JZ = 1;
  const int* DESCZ = outEigVec->pDESC_;

  int LWORK = -1;
  double* pWorkSizeCheck = new double[1];

  const int NPCOL = inMatrix.procGridCol_;
  const int LIWORK = 7 * N + 8 * NPCOL + 2;
  int* IWORK = new int[LIWORK];
  int info = 0;

  pdsyevd_(JOBZ, UPLO, &N, pBufA, &IA, &JA, DESCA, W, Z, &IZ, &JZ, DESCZ,
           pWorkSizeCheck, &LWORK, IWORK, &LIWORK, &info);

  if (info != 0) {
    logger.critical(TlUtils::format("pdsyev_ error @1st call: %d", info));
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
    logger.critical(TlUtils::format("pdsyev_ error @2nd call: %d", info));
  }
  return ((info == 0) ? true : false);
}
