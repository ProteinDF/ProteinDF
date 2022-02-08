#include "tl_dense_general_matrix_impl_scalapack.h"

#include "TlCommunicate.h"
#include "scalapack.h"
#include "tl_dense_general_matrix_io.h"
#include "tl_dense_general_matrix_object.h"
#include "tl_dense_vector_impl_scalapack.h"
#include "tl_matrix_utils.h"
#include "tl_scalapack_context.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

const std::size_t TlDenseGeneralMatrix_ImplScalapack::FILE_BUFFER_SIZE = 100 * 1024 * 1024;  // 100 MB
bool TlDenseGeneralMatrix_ImplScalapack::isUsingPartialIO = false;

/// @file
/// OpenMPに関する注意事項
///
/// #pragma omp critical (TlDenseGeneralMatrix_ImplScalapack_gpmServerTasks_) は
/// this->gpmServerTasks_ メンバ変数の書き換え処理を行う箇所に設置している。
/// 同様に #pragma omp critical
/// (TlDenseGeneralMatrix_ImplScalapack_gpmClientTasks_) は
/// this->gpmClientTasks_ メンバ変数の書き換え処理を行う箇所に設置している。

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplScalapack::TlDenseGeneralMatrix_ImplScalapack(const TlMatrixObject::index_type row,
                                                                       const TlMatrixObject::index_type col)
    : context_(0), rows_(row), cols_(col), blockSize_(TlScalapackContext::getBlockSize()), pData_(NULL) {
    this->initialize();
}

TlDenseGeneralMatrix_ImplScalapack::TlDenseGeneralMatrix_ImplScalapack(const TlDenseGeneralMatrix_ImplScalapack& rhs)
    : context_(0), rows_(rhs.rows_), cols_(rhs.cols_), blockSize_(rhs.blockSize_), pData_(NULL) {
    this->initialize();
    if (rhs.getNumOfMyElements() > 0) {
        std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(), this->pData_);
    }
}

TlDenseGeneralMatrix_ImplScalapack::~TlDenseGeneralMatrix_ImplScalapack() {
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getNumOfRows() const {
    return this->rows_;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getNumOfCols() const {
    return this->cols_;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplScalapack::getNumOfElements() const {
    const TlMatrixObject::size_type rows = this->getNumOfRows();
    const TlMatrixObject::size_type cols = this->getNumOfCols();
    return rows * cols;
}

void TlDenseGeneralMatrix_ImplScalapack::resize(const TlMatrixObject::index_type row,
                                                const TlMatrixObject::index_type col) {
    assert(row > 0);
    assert(col > 0);
    if ((row == this->getNumOfRows()) && (col == this->getNumOfCols())) {
        // do not need operation.
        return;
    }

    // backup old condition
    TlDenseGeneralMatrix_ImplScalapack tmp(*this);

    // new size and zero clear
    this->rows_ = row;
    this->cols_ = col;
    this->initialize();

    // copy
    const int copyMaxRows = std::min(row, tmp.getNumOfRows());
    const int copyMaxCols = std::min(col, tmp.getNumOfCols());
    assert(this->getBlockSize() == tmp.getBlockSize());
    {
        // ブロックサイズが同じでプロセスグリッドが同一なら
        // ローカルコピーで十分。
        std::vector<TlMatrixObject::index_type>::const_iterator rBegin = tmp.m_RowIndexTable.begin();
        std::vector<TlMatrixObject::index_type>::const_iterator cBegin = tmp.m_ColIndexTable.begin();
        std::vector<TlMatrixObject::index_type>::const_iterator rEnd =
            std::lower_bound(tmp.m_RowIndexTable.begin(), tmp.m_RowIndexTable.end(), copyMaxRows);
        std::vector<TlMatrixObject::index_type>::const_iterator cEnd =
            std::lower_bound(tmp.m_ColIndexTable.begin(), tmp.m_ColIndexTable.end(), copyMaxCols);
        for (std::vector<TlMatrixObject::index_type>::const_iterator r = rBegin; r != rEnd; ++r) {
            const TlMatrixObject::index_type row = *r;
            for (std::vector<TlMatrixObject::index_type>::const_iterator c = cBegin; c != cEnd; ++c) {
                const TlMatrixObject::index_type col = *c;
                const TlMatrixObject::index_type oldIndex =
                    std::distance(rBegin, r) + std::distance(cBegin, c) * tmp.myRows_;
                assert(0 <= oldIndex);
                assert(oldIndex < tmp.getNumOfMyElements());

                const double value = tmp.pData_[oldIndex];
                this->set(row, col, value);
            }
        }
    }
}

double TlDenseGeneralMatrix_ImplScalapack::get(const TlMatrixObject::index_type globalRow,
                                               const TlMatrixObject::index_type globalCol) const {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    // DataType& data_tmp = const_cast<DataType&>(this->data_);

    double answer = 0.0;
    const int fortranGlobalRow = globalRow + 1;
    const int fortranGlobalCol = globalCol + 1;
    pdelget_("A", " ", &answer, this->pData_, &fortranGlobalRow, &fortranGlobalCol, this->pDESC_);

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::set(const TlMatrixObject::index_type globalRow,
                                             const TlMatrixObject::index_type globalCol, const double value) {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        assert(this->getIndex0(globalRow, globalCol) != -1);
        assert(this->getIndex0(globalRow, globalCol) == index);
        this->pData_[index] = value;
    }
}

void TlDenseGeneralMatrix_ImplScalapack::add(const TlMatrixObject::index_type globalRow,
                                             const TlMatrixObject::index_type globalCol, const double value) {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        this->pData_[index] += value;
    }
}

double TlDenseGeneralMatrix_ImplScalapack::getLocal(const TlMatrixObject::index_type globalRow,
                                                    const TlMatrixObject::index_type globalCol) const {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    double answer = 0.0;

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        answer = this->pData_[index];
    }

    return answer;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator=(
    const TlDenseGeneralMatrix_ImplScalapack& rhs) {
    if (&rhs != this) {
        assert(rhs.getNumOfRows() > 0);
        assert(rhs.getNumOfCols() > 0);
        this->rows_ = rhs.rows_;
        this->cols_ = rhs.cols_;
        this->blockSize_ = rhs.blockSize_;

        this->initialize();
        std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(), this->pData_);
    }

    return (*this);
}

TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator+=(
    const TlDenseGeneralMatrix_ImplScalapack& rhs) {
    assert(this->rows_ == rhs.rows_);
    assert(this->cols_ == rhs.cols_);
    assert(this->blockSize_ == rhs.blockSize_);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] += rhs.pData_[i];
    }

    return (*this);
}

TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator-=(
    const TlDenseGeneralMatrix_ImplScalapack& rhs) {
    assert(this->rows_ == rhs.rows_);
    assert(this->cols_ == rhs.cols_);
    assert(this->blockSize_ == rhs.blockSize_);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] -= rhs.pData_[i];
    }

    return (*this);
}

TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator*=(const double coef) {
    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] *= coef;
    }

    return (*this);
}

TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator/=(const double coef) {
    return this->operator*=(1.0 / coef);
}

TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::operator*=(
    const TlDenseGeneralMatrix_ImplScalapack& rhs) {
    TlDenseGeneralMatrix_ImplScalapack tmp = *this;
    *this = tmp * rhs;

    return (*this);
}

TlDenseGeneralMatrix_ImplScalapack TlDenseGeneralMatrix_ImplScalapack::inverse() const {
    TlDenseGeneralMatrix_ImplScalapack A = *this;

    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();
    const int IA = 1;
    const int JA = 1;
    const int zero = 0;
    int info = 0;
    const int sizeOf_IPIV = numroc_(&A.rows_, &A.blockSize_, &A.m_nMyProcRow, &zero, &A.procGridRow_) + A.blockSize_;

    int* IPIV = new int[sizeOf_IPIV];

    pdgetrf_(&M, &N, A.pData_, &IA, &JA, A.pDESC_, IPIV, &info);

    if (info == 0) {
        int LWORK = -1;
        int LIWORK = -1;
        double* WORK_SIZE = new double[1];
        int* IWORK_SIZE = new int[1];
        pdgetri_(&M, A.pData_, &IA, &JA, A.pDESC_, IPIV, WORK_SIZE, &LWORK, IWORK_SIZE, &LIWORK, &info);

        LWORK = static_cast<int>(WORK_SIZE[0]);
        double* WORK = new double[LWORK];
        LIWORK = IWORK_SIZE[0];
        int* IWORK = new int[LIWORK];
        delete[] WORK_SIZE;
        WORK_SIZE = NULL;
        delete[] IWORK_SIZE;
        IWORK_SIZE = NULL;

        pdgetri_(&M, A.pData_, &IA, &JA, A.pDESC_, IPIV, WORK, &LWORK, IWORK, &LIWORK, &info);

        if (info != 0) {
            std::cout << "pdgetri_ returns " << info << std::endl;
            throw;
        }

        delete[] IWORK;
        IWORK = NULL;

        delete[] WORK;
        WORK = NULL;
    } else {
        std::cout << "pdgetrf_ returns " << info << std::endl;
        throw;
    }

    delete[] IPIV;
    IPIV = NULL;

    return A;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
std::vector<double> TlDenseGeneralMatrix_ImplScalapack::diagonals() const {
    const TlMatrixObject::index_type dim = std::min(this->getNumOfRows(), this->getNumOfCols());
    std::vector<double> answer(dim, 0.0);

    for (TlMatrixObject::index_type i = 0; i < dim; ++i) {
        const double value = this->getLocal(i, i);
        answer[i] = value;
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(&(answer[0]), dim);

    return answer;
}

double TlDenseGeneralMatrix_ImplScalapack::sum() const {
    double answer = 0.0;
    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        answer += this->pData_[i];
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);
    return answer;
}

double TlDenseGeneralMatrix_ImplScalapack::trace() const {
    double answer = 0.0;
    const int nGlobalCols = this->cols_;
    const int nRows = this->myRows_;
    const int nCols = this->myCols_;
    for (int c = 0; c < nCols; ++c) {
        const int nGlobalColIndex = this->m_ColIndexTable[c];
        if (nGlobalColIndex >= nGlobalCols) {
            continue;
        }

        for (int r = 0; r < nRows; ++r) {
            const int nGlobalRowIndex = this->m_RowIndexTable[r];
            if (nGlobalRowIndex == nGlobalColIndex) {
                const int index = r + nRows * c;  // row-major

                answer += this->pData_[index];
            }
        }
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(answer);

    return answer;
}

double TlDenseGeneralMatrix_ImplScalapack::getRMS() const {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double sum2 = 0.0;
    const TlMatrixObject::size_type maxIndex = this->getNumOfMyElements();
    for (TlMatrixObject::size_type i = 0; i < maxIndex; ++i) {
        double tmp = this->pData_[i];
        sum2 += tmp * tmp;
    }
    rComm.allReduce_SUM(sum2);

    const double elements = this->getNumOfRows() * this->getNumOfCols();

    const double rms = std::sqrt(sum2 / elements);
    return rms;
}

double TlDenseGeneralMatrix_ImplScalapack::getMaxAbsoluteElement(TlMatrixObject::index_type* pOutRow,
                                                                 TlMatrixObject::index_type* pOutCol) const {
    TlMatrixObject::index_type row = 0;
    TlMatrixObject::index_type col = 0;
    double maxValue = this->getLocalMaxAbsoluteElement(&row, &col);

    // communication
    std::vector<TlMatrixObject::index_type> rowcol(2);
    rowcol[0] = row;
    rowcol[1] = col;
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int proc = rComm.getNumOfProc();
        const int rank = rComm.getRank();
        for (int i = 1; i < proc; ++i) {
            if (i == rank) {
                rComm.sendData(maxValue);
                rComm.sendData(rowcol);
            } else if (rComm.isMaster() == true) {
                double rec_maxValue = 0.0;
                std::vector<TlMatrixObject::index_type> rec_rowcol;
                rComm.receiveData(rec_maxValue, i);
                rComm.receiveData(rec_rowcol, i);
                if (maxValue < rec_maxValue) {
                    maxValue = rec_maxValue;
                    row = rowcol[0];
                    col = rowcol[1];
                }
            }
        }

        rComm.broadcast(maxValue);
        rComm.broadcast(rowcol);
    }

    if (pOutRow != NULL) {
        *pOutRow = rowcol[0];
    }
    if (pOutCol != NULL) {
        *pOutCol = rowcol[1];
    }

    return maxValue;
}

TlDenseGeneralMatrix_ImplScalapack TlDenseGeneralMatrix_ImplScalapack::transpose() const {
    TlDenseGeneralMatrix_ImplScalapack answer = *this;
    answer.transposeInPlace();

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::transposeInPlace() {
    const int M = this->cols_;
    const int N = this->rows_;
    const double alpha = 1.0;
    const double beta = 0.0;
    const int IA = 1;
    const int JA = 1;
    TlDenseGeneralMatrix_ImplScalapack C(M, N);
    const int IC = 1;
    const int JC = 1;
    pdtran_(&M, &N, &alpha, &(this->pData_[0]), &IA, &JA, this->pDESC_, &beta, &(C.pData_[0]), &IC, &JC, C.pDESC_);

    (*this) = C;
}

TlDenseGeneralMatrix_ImplScalapack TlDenseGeneralMatrix_ImplScalapack::dot(
    const TlDenseGeneralMatrix_ImplScalapack& X) const {
    TlDenseGeneralMatrix_ImplScalapack answer = *this;
    answer.dotInPlace(X);

    return answer;
}

// pddot?
const TlDenseGeneralMatrix_ImplScalapack& TlDenseGeneralMatrix_ImplScalapack::dotInPlace(
    const TlDenseGeneralMatrix_ImplScalapack& X) {
    assert(this->getNumOfRows() == X.getNumOfRows());
    assert(this->getNumOfCols() == X.getNumOfCols());

    const TlMatrixObject::size_type size = this->getNumOfMyElements();
    for (TlMatrixObject::size_type index = 0; index < size; ++index) {
        this->pData_[index] *= X.pData_[index];
    }

    return (*this);
}

// 複数スレッドから同時に呼び出さない
bool TlDenseGeneralMatrix_ImplScalapack::getSparseMatrix(TlSparseMatrix* pMatrix, bool isFinalize) const {
    bool answer = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (pMatrix != NULL) {
// Client task通信処理
#pragma omp critical(TlDenseGeneralMatrix_Scalapack__getSparseMatrixX)
        {
            this->getSparseMatrix_registerTask(pMatrix);
            answer = this->getPartialMatrix_ClientTasks(pMatrix);
        }
    } else {
        // Server処理
        this->getPartialMatrix_ServerTasks(isFinalize);
    }

    if (isFinalize == true) {
        rComm.checkNonBlockingCommunications();
        assert(this->gpmClientTasks_.size() == 0);

#pragma omp critical(TlDenseGeneralMatrix_blacs_gpmClientTasks_)
        { this->gpmClientTasks_.clear(); }
    }

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::mergeSparseMatrix(const TlSparseMatrix& M) {
    assert(M.getNumOfRows() == this->getNumOfRows());
    assert(M.getNumOfCols() == this->getNumOfCols());

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    rComm.checkNonBlockingCommunications();

    // 送信すべきインデックスリストの作成
    std::vector<std::vector<TlMatrixObject::index_type> > indexArrays(numOfProcs);
    std::vector<std::vector<double> > values(numOfProcs);

    TlSparseMatrix::const_iterator itEnd = M.end();
    for (TlSparseMatrix::const_iterator it = M.begin(); it != itEnd; ++it) {
        TlMatrixObject::index_type globalRow = it->first.row;
        TlMatrixObject::index_type globalCol = it->first.col;
        const double value = it->second;

        int rank = 0;
        TlMatrixObject::index_type myRow = 0;
        TlMatrixObject::index_type myCol = 0;
        this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
        if (rank == this->rank_) {
            const TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
            this->pData_[index] += value;
        } else {
            indexArrays[rank].push_back(globalRow);
            indexArrays[rank].push_back(globalCol);
            values[rank].push_back(value);
        }
    }

    this->mergeMatrix_common(indexArrays, values);
    rComm.checkNonBlockingCommunications();
}

std::vector<TlMatrixObject::index_type> TlDenseGeneralMatrix_ImplScalapack::getRowIndexTable() const {
    return this->m_RowIndexTable;
}

std::vector<TlMatrixObject::index_type> TlDenseGeneralMatrix_ImplScalapack::getColIndexTable() const {
    return this->m_ColIndexTable;
}

void TlDenseGeneralMatrix_ImplScalapack::getLocalMatrix(TlDenseGeneralMatrixObject* pOutMatrix) const {
    const std::vector<TlMatrixObject::index_type> rowIndexes = this->getRowIndexTable();
    const std::vector<TlMatrixObject::index_type> colIndexes = this->getColIndexTable();
    const int numOfRowIndexes = rowIndexes.size();
    const int numOfColIndexes = colIndexes.size();
    pOutMatrix->resize(numOfRowIndexes, numOfColIndexes);

    for (int rowIndex = 0; rowIndex < numOfRowIndexes; ++rowIndex) {
        const TlMatrixObject::index_type row = rowIndexes[rowIndex];
        for (int colIndex = 0; colIndex < numOfColIndexes; ++colIndex) {
            const TlMatrixObject::index_type col = colIndexes[colIndex];
            pOutMatrix->set(rowIndex, colIndex, this->getLocal(row, col));
        }
    }
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseGeneralMatrix_ImplScalapack::load(const std::string& sFilePath) {
    // if (TlDenseGeneralMatrix_blacs::isUsingPartialIO == true) {
    //   return this->loadLocal(sFilePath);
    // }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool bAnswer = false;

    std::ifstream ifs;
    if (rComm.isMaster() == true) {
        ifs.open(sFilePath.c_str(), std::ios::in | std::ios::binary);
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
            std::cerr << "TlDenseGeneralMatrix_ImplScalapack::load() is not "
                         "supported: "
                      << sFilePath << std::endl;
        }
        std::abort();
    }

    if (rComm.isMaster() == true) {
        ifs.close();
    }

    return bAnswer;
}

bool TlDenseGeneralMatrix_ImplScalapack::load(std::ifstream& ifs) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());

    const int numOfProcs = rComm.getNumOfProc();

    // read header
    bool answer = true;
    TlMatrixObject::HeaderInfo headerInfo;
    TlMatrixObject::MatrixType matrixType = TlMatrixObject::UNDEFINED;
    std::vector<TlMatrixObject::index_type> rowcol(2);
    if (rComm.isMaster() == true) {
        const bool isLoadable = TlMatrixUtils::getHeaderInfo(ifs, &headerInfo);
        assert((headerInfo.matrixType == TlMatrixObject::CSFD) || (headerInfo.matrixType == TlMatrixObject::RSFD));
        matrixType = headerInfo.matrixType;
        rowcol[0] = headerInfo.numOfRows;
        rowcol[1] = headerInfo.numOfCols;
        // std::cerr << TlUtils::format("%d (%d, %d)", (int)matrixType, rowcol[0], rowcol[1]) << std::endl;
    }
    rComm.broadcast(&(rowcol[0]), 2, 0);

    this->rows_ = rowcol[0];
    this->cols_ = rowcol[1];
    this->initialize();

    if (rComm.isMaster() == true) {
        static const std::size_t bufferCount = FILE_BUFFER_SIZE / sizeof(double);
        std::vector<double> buf(bufferCount, 0.0);

        TlMatrixObject::index_type row = 0;
        TlMatrixObject::index_type col = 0;
        TlMatrixObject::size_type count = 0;
        const TlMatrixObject::size_type maxCount = this->getNumOfRows() * this->getNumOfCols();
        bool isFinished = false;

        std::vector<TlMatrixObject::index_type> sizeLists(numOfProcs, 0);
        std::vector<std::vector<TlMatrixObject::MatrixElement> > elementsBuf(numOfProcs);
        std::vector<bool> isSendData(numOfProcs, false);

        while (isFinished == false) {
            // buffer分を一度に読み込み
            ifs.read((char*)(&buf[0]), sizeof(double) * bufferCount);

            // 各プロセスのバッファに振り分ける
            std::vector<std::vector<TlMatrixObject::MatrixElement> > elements(numOfProcs);
            for (std::size_t i = 0; i < bufferCount; ++i) {
                if (matrixType == TlMatrixObject::RSFD) {
                    const std::ldiv_t d = std::ldiv(count, this->getNumOfCols());
                    row = d.quot;
                    col = d.rem;
                } else {
                    assert(matrixType == TlMatrixObject::CSFD);
                    const std::ldiv_t d = std::ldiv(count, this->getNumOfRows());
                    col = d.quot;
                    row = d.rem;
                }
                // const int proc = this->getProcIdForIndex(row, col);
                int rank = 0;
                this->getGlobalRowCol2LocalRowCol(row, col, &rank);
                elements[rank].push_back(TlMatrixObject::MatrixElement(row, col, buf[i]));

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
                    rComm.iSendDataX(&(elementsBuf[proc][0]), sizeLists[proc], proc, TAG_LOAD_VALUES);
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
                    TlMatrixObject::index_type myRow;
                    TlMatrixObject::index_type myCol;
                    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
                    assert(rank == this->rank_);
                    TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
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
        // int endMsg = 0;

        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
        // rComm.iReceiveData(endMsg, root, TAG_LOAD_END);
        bool isLoopBreak = false;
        while (isLoopBreak == false) {
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
                            TlMatrixObject::index_type myRow;
                            TlMatrixObject::index_type myCol;
                            this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
                            assert(rank == this->rank_);
                            TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
                            this->pData_[index] = elements[i].value;
                        }
                    }
                    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
                } else {
                    assert(sizeList < 0);
                    isLoopBreak = true;
                }
            }

            // if (rComm.test(endMsg) == true) {
            //   rComm.wait(endMsg);
            //   rComm.cancel(sizeList);
            //   isLoopBreak = true;
            // }
        }
    }

    assert(rComm.checkNonBlockingCommunications());
    return answer;
}

bool TlDenseGeneralMatrix_ImplScalapack::save(const std::string& filePath) const {
    // if (TlDenseGeneralMatrix_blacs::isUsingPartialIO == true) {
    //   return this->saveLocal(sFilePath);
    // }

    bool answer = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());

    if (rComm.isMaster() == true) {
        // master
        const TlMatrixObject::index_type globalRows = this->getNumOfRows();
        const TlMatrixObject::index_type globalCols = this->getNumOfCols();
        TlFileGenericMatrix fm(filePath, globalRows, globalCols);

        // store local matrix
        {
            const std::vector<TlMatrixObject::MatrixElement> buf = this->getMatrixElementsInLocal2();
            this->saveElements(&fm, buf);
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

            this->saveElements(&fm, buf);
        }
    } else {
        // slave: send submatrix
        const int root = 0;
        const std::vector<TlMatrixObject::MatrixElement> buf = this->getMatrixElementsInLocal2();
        const TlMatrixObject::size_type bufSize = buf.size();
        rComm.sendData(bufSize, root, TAG_SAVE_HANDSHAKE);
        rComm.iSendDataX(&(buf[0]), buf.size(), root, TAG_SAVE_DATA);
        rComm.wait(&(buf[0]));
    }

    rComm.broadcast(answer);
    assert(rComm.checkNonBlockingCommunications());
    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::setUsingPartialIO(bool isUsePIO) {
    TlDenseGeneralMatrix_ImplScalapack::isUsingPartialIO = isUsePIO;
}

// void TlDenseGeneralMatrix_ImplScalapack::dump(
//     TlDenseVector_ImplScalapack* v) const {
//   v->resize(this->getNumOfRows() * this->getNumOfCols());
//   std::copy(this->pData_, this->pData_ + this->getNumOfMyElements(),
//             v->vector_);
// }

// void TlDenseGeneralMatrix_ImplScalapack::restore(
//     const TlDenseVector_ImplScalapack& v) {
//   const std::size_t bufSize = this->getNumOfMyElements();
//   std::copy(v.vector_, v.vector_ + bufSize, this->pData_);
// }

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseGeneralMatrix_ImplScalapack::initialize() {
    TlScalapackContext::getData(this->context_, this->proc_, this->rank_, this->procGridRow_, this->procGridCol_);

    // my process position on the process matrix
    Cblacs_gridinfo(this->context_, &(this->procGridRow_), &(this->procGridCol_), &(this->m_nMyProcRow),
                    &(this->m_nMyProcCol));

    // determine sizes of local matrix
    const int nStartRowProc = 0;
    const int nStartColProc = 0;
    this->myRows_ = std::max(
        1, numroc_(&(this->rows_), &(this->blockSize_), &(this->m_nMyProcRow), &nStartRowProc, &(this->procGridRow_)));
    this->myCols_ = std::max(
        1, numroc_(&(this->cols_), &(this->blockSize_), &(this->m_nMyProcCol), &nStartColProc, &(this->procGridCol_)));

    // make parameter, desca
    int nInfo = 0;
    descinit_(this->pDESC_, &(this->rows_), &(this->cols_), &(this->blockSize_), &(this->blockSize_), &nStartRowProc,
              &nStartColProc, &(this->context_), &(this->myRows_), &nInfo);
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
        const int nMyRows = this->myRows_;
        const int nBlockSize = this->blockSize_;
        const int nBlockIndex = this->m_nMyProcRow * nBlockSize;           // 各ローカル行列の最初のインデックス
        const int nIncrementBlockIndex = this->procGridRow_ * nBlockSize;  // ブロック最初のインデックスの増分
        this->m_RowIndexTable.clear();
        this->m_RowIndexTable.reserve(nMyRows);
        for (int r = 0; r < nMyRows; ++r) {
            const div_t d = std::div(r, nBlockSize);
            const int i = nBlockIndex + (nIncrementBlockIndex * d.quot) + d.rem;
            assert(i == this->getLocal2Global_row(r));
            if (i < this->rows_) {
                this->m_RowIndexTable.push_back(i);
            } else {
                break;
            }
        }
        std::vector<int>(this->m_RowIndexTable).swap(this->m_RowIndexTable);
    }

    // 列方向のglobal_index v.s. local_indexのリストを作成
    {
        const int nMyCols = this->myCols_;
        const int nBlockSize = this->blockSize_;
        const int nBlockIndex = this->m_nMyProcCol * this->blockSize_;           // 各ローカル行列の最初のインデックス
        const int nIncrementBlockIndex = this->procGridCol_ * this->blockSize_;  // ブロック最初のインデックスの増分
        this->m_ColIndexTable.clear();
        this->m_ColIndexTable.reserve(nMyCols);
        for (int c = 0; c < nMyCols; ++c) {
            const div_t d = std::div(c, nBlockSize);
            const int i = nBlockIndex + (nIncrementBlockIndex * d.quot) + d.rem;
            assert(i == this->getLocal2Global_col(c));
            if (i < this->cols_) {
                this->m_ColIndexTable.push_back(i);
            } else {
                break;
            }
        }
        std::vector<int>(this->m_ColIndexTable).swap(this->m_ColIndexTable);
    }

    //   // getPartialMatrix通信用 初期化処理
    //   this->isDebugOut_GPM_ = false;
    //
    // #pragma omp critical(TlDenseGeneralMatrix_ImplScalapack_gpmServerTasks_)
    //   { this->gpmServerTasks_.clear(); }
    //   this->isWaitingRequestHandShake_ = false;
    //   this->sessionId_ = 0;
    //
    //   this->trafficControl_.clear();
    //   this->trafficControl_.resize(numOfProcs);
    //   for (int i = 0; i < numOfProcs; ++i) {
    //     this->trafficControl_[i] = NULL;
    //   }
    //
    //   this->sessionTable_.resize(numOfProcs);
    //   for (int proc = 0; proc < numOfProcs; ++proc) {
    //     this->sessionTable_[proc].reset();
    //   }
    //
    //   this->gpmClientTasks_.clear();
    //
    //   // mergeMatrix用
    //   this->mm_isWaitingSessionId_ = false;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplScalapack::getNumOfMyElements() const {
    const TlMatrixObject::size_type rows = this->myRows_;
    const TlMatrixObject::size_type cols = this->myCols_;

    return rows * cols;
}

int TlDenseGeneralMatrix_ImplScalapack::getBlockSize() const {
    return this->blockSize_;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_ImplScalapack::getIndex0(
    const TlMatrixObject::index_type nGlobalRow, const TlMatrixObject::index_type nGlobalCol) const {
    TlMatrixObject::size_type answer = -1;

    std::vector<TlMatrixObject::index_type>::const_iterator pRow =
        std::lower_bound(this->m_RowIndexTable.begin(), this->m_RowIndexTable.end(), nGlobalRow);
    if ((pRow != this->m_RowIndexTable.end()) && (*pRow == nGlobalRow)) {
        std::vector<TlMatrixObject::index_type>::const_iterator pCol =
            std::lower_bound(this->m_ColIndexTable.begin(), this->m_ColIndexTable.end(), nGlobalCol);

        if ((pCol != this->m_ColIndexTable.end()) && (*pCol == nGlobalCol)) {
            const TlMatrixObject::size_type localRow = pRow - this->m_RowIndexTable.begin();
            const TlMatrixObject::size_type localCol = pCol - this->m_ColIndexTable.begin();

            answer = localRow + localCol * (this->myRows_);
        }
    }

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::getGlobalRowCol2LocalRowCol(const TlMatrixObject::index_type globalRow,
                                                                     const TlMatrixObject::index_type globalCol,
                                                                     int* pRank, TlMatrixObject::index_type* pLocalRow,
                                                                     TlMatrixObject::index_type* pLocalCol) const {
    const std::div_t rowBlockDiv = std::div(globalRow, this->blockSize_);
    const std::div_t colBlockDiv = std::div(globalCol, this->blockSize_);
    const int rowBlockIndex = rowBlockDiv.quot;
    const int colBlockIndex = colBlockDiv.quot;
    const int rowLocalIndex = rowBlockDiv.rem;
    const int colLocalIndex = colBlockDiv.rem;

    const std::div_t rowIndexResult = std::div(rowBlockIndex, this->procGridRow_);
    const std::div_t colIndexResult = std::div(colBlockIndex, this->procGridCol_);
    int rowProcCycle = rowIndexResult.quot;
    int colProcCycle = colIndexResult.quot;
    int rowProcID = rowIndexResult.rem;
    int colProcID = colIndexResult.rem;

    if (pRank != NULL) {
        // "Row-major"
        const int rank = rowProcID * this->procGridCol_ + colProcID;
        *pRank = rank;
    }
    if (pLocalRow != NULL) {
        *pLocalRow = rowProcCycle * this->blockSize_ + rowLocalIndex;
    }
    if (pLocalCol != NULL) {
        *pLocalCol = colProcCycle * this->blockSize_ + colLocalIndex;
    }
}

// void TlDenseGeneralMatrix_ImplScalapack::getLocalRowCol2GlobalRowCol(
//     const TlMatrixObject::index_type localRow,
//     const TlMatrixObject::index_type localCol,
//     TlMatrixObject::index_type* pGlobalRow,
//     TlMatrixObject::index_type* pGlobalCol) const {
//   const int blockSize = this->getBlockSize();
//
//   const int offsetRow = this->m_nMyProcRow * blockSize;
//   const int offsetCol = this->m_nMyProcCol * blockSize;
//   const int blocksRow = this->procGridRow_ * blockSize;
//   const int blocksCol = this->procGridCol_ * blockSize;
//
//   const std::div_t divRow = std::div(localRow, blockSize);
//   const std::div_t divCol = std::div(localCol, blockSize);
//   *pGlobalRow = offsetRow + blocksRow * divRow.quot + divRow.rem;
//   *pGlobalCol = offsetCol + blocksCol * divCol.quot + divCol.rem;
// }

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getLocal2Global_row(
    const TlMatrixObject::index_type localIndex) const {
    return this->getLocal2Global(localIndex, this->m_nMyProcRow, this->procGridRow_);
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getLocal2Global_col(
    const TlMatrixObject::index_type localIndex) const {
    return this->getLocal2Global(localIndex, this->m_nMyProcCol, this->procGridCol_);
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getLocal2Global(
    const TlMatrixObject::index_type localIndex, const int procID, const int numOfProcs) const {
    const int blockSize = this->getBlockSize();

    const int offset = procID * blockSize;
    const int distance = numOfProcs * blockSize;

    const std::div_t div_result = std::div(localIndex, blockSize);
    const int globalIndex = offset + distance * div_result.quot + div_result.rem;

    return globalIndex;
}

TlMatrixObject::index_type TlDenseGeneralMatrix_ImplScalapack::getLocalIndex(
    const TlMatrixObject::index_type localRow, const TlMatrixObject::index_type localCol) const {
    const TlMatrixObject::size_type index = localRow + this->myRows_ * localCol;
    return index;
}

// int TlDenseGeneralMatrix_ImplScalapack::getProcIdForIndex(
//     const TlMatrixObject::index_type globalRow,
//     const TlMatrixObject::index_type globalCol) const {
//   const int rowBlockIndex = globalRow / this->blockSize_;
//   const int rowProcID = rowBlockIndex % this->procGridRow_;
//   const int colBlockIndex = globalCol / this->blockSize_;
//   const int colProcID = colBlockIndex % this->procGridCol_;
//
//   // "Row-major" for Cblacs_gridinit
//   const int procMatrixIndex = rowProcID * this->procGridCol_ + colProcID;
//   // const int procID = this->processMatrix_[procMatrixIndex];
//   // int procID = procMatrixIndex;
//
//   return procID;
// }

double TlDenseGeneralMatrix_ImplScalapack::getLocalMaxAbsoluteElement(TlMatrixObject::index_type* pOutRow,
                                                                      TlMatrixObject::index_type* pOutCol) const {
    double dMaxValue = 0.0;
    std::size_t index = 0;
    const TlMatrixObject::size_type nMaxIndex = this->getNumOfMyElements();
    for (TlMatrixObject::size_type i = 0; i < nMaxIndex; ++i) {
        const double v = std::fabs(this->pData_[i]);
        if (dMaxValue < v) {
            dMaxValue = v;
            index = i;
        }
    }

    if ((pOutRow != NULL) || (pOutCol != NULL)) {
        ldiv_t t = std::ldiv(index, this->myRows_);
        TlMatrixObject::index_type myRow = t.rem;
        TlMatrixObject::index_type myCol = t.quot;
        TlMatrixObject::index_type nRow = 0;
        TlMatrixObject::index_type nCol = 0;
        if ((myRow < static_cast<TlMatrixObject::index_type>(this->m_RowIndexTable.size())) &&
            (myCol < static_cast<TlMatrixObject::index_type>(this->m_ColIndexTable.size()))) {
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

void TlDenseGeneralMatrix_ImplScalapack::getPartialMatrix_ServerTasks(const bool isFinalize) const {
#pragma omp critical(TlDenseGeneralMatrix_Scalapack__getPartialMatrix_ServerTasks)
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();

        int requestProc = 0;
        if (this->isWaitingRequestHandShake_ == false) {
            // セッション開始
            rComm.iReceiveDataFromAnySource(this->sessionId_, TAG_GPM_SESSION_ID);
            this->isWaitingRequestHandShake_ = true;

            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] SRV [**]: recv session id from any.", rComm.getRank()) << std::endl;
            }
        }

        if ((this->isWaitingRequestHandShake_ == true) && (rComm.test(this->sessionId_, &requestProc) == true)) {
            // セッション番号受信
            rComm.wait(this->sessionId_);
            this->isWaitingRequestHandShake_ = false;

            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] SRV [OK]: recv session id from [%d].", rComm.getRank(), requestProc)
                          << std::endl;
            }

            // taskの登録
            GPM_ServerTask task;
            task.state = 0;
            task.sessionId = this->sessionId_;
            task.requestProc = requestProc;
            task.numOfComponents = 0;

            this->gpmServerTasks_.push_back(task);
        }

        GpmServerTasksType::iterator itEnd = this->gpmServerTasks_.end();
        for (GpmServerTasksType::iterator it = this->gpmServerTasks_.begin(); it != itEnd; ++it) {
            GPM_ServerTask& task = *it;
            const int proc = task.requestProc;
            if ((task.state & GPM_SERVER_RECV_NUM_OF_COMPONENTS) == 0) {
                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format(
                                     "[%d] SRV [**]: recv num of components "
                                     "from [%d]",
                                     rComm.getRank(), proc)
                              << std::endl;
                }

                rComm.iReceiveData(task.numOfComponents, proc, TAG_GPM_NUM_OF_COMPONENTS + task.sessionId);
                task.state |= GPM_SERVER_RECV_NUM_OF_COMPONENTS;
                continue;
            }

            if (((task.state & GPM_SERVER_RECV_NUM_OF_COMPONENTS) != 0) &&
                ((task.state & GPM_SERVER_WAIT_NUM_OF_COMPONENTS) == 0)) {
                if (rComm.test(task.numOfComponents) == true) {
                    rComm.wait(task.numOfComponents);
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format(
                                         "[%d] SRV [OK]: recv num of components from "
                                         "[%d]. "
                                         "num = %d",
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
                    std::cerr << TlUtils::format("[%d] SRV [**]: recv components from [%d]", rComm.getRank(), proc)
                              << std::endl;
                }

                assert(task.numOfComponents > 0);
                task.components.resize(task.numOfComponents);
                rComm.iReceiveDataX((int*)&(task.components[0]), task.numOfComponents, proc,
                                    TAG_GPM_COMPONENTS + task.sessionId);
                task.state |= GPM_SERVER_RECV_COMPONENTS;
                continue;
            }

            if (((task.state & GPM_SERVER_RECV_COMPONENTS) != 0) && ((task.state & GPM_SERVER_WAIT_COMPONENTS) == 0)) {
                if (rComm.test((int*)&(task.components[0])) == true) {
                    rComm.wait((int*)&(task.components[0]));
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] SRV [OK]: recv components from [%d]", rComm.getRank(), proc)
                                  << std::endl;
                    }

                    task.state |= GPM_SERVER_WAIT_COMPONENTS;
                    continue;
                }
            }

            if (((task.state & GPM_SERVER_WAIT_COMPONENTS) != 0) &&
                ((task.state & GPM_SERVER_SEND_ELEMENT_VALUES) == 0)) {
                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] SRV [**]: send elements to [%d].", rComm.getRank(), proc)
                              << std::endl;
                }

                const int size = task.numOfComponents / 2;
                task.elementValues.resize(size);
                for (int i = 0; i < size; ++i) {
                    const TlMatrixObject::index_type globalRow = task.components[i * 2];
                    const TlMatrixObject::index_type globalCol = task.components[i * 2 + 1];

                    int rank = 0;
                    TlMatrixObject::index_type myRow = 0;
                    TlMatrixObject::index_type myCol = 0;
                    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
                    assert(rank == this->rank_);
                    TlMatrixObject::size_type localIndex = this->getLocalIndex(myRow, myCol);
                    const double value = this->pData_[localIndex];
                    task.elementValues[i] = value;
                }
                rComm.iSendDataX((double*)&(task.elementValues[0]), size, proc,
                                 TAG_GPM_ELEMENT_VALUES + task.sessionId);
                task.state |= GPM_SERVER_SEND_ELEMENT_VALUES;
                continue;
            }

            if (((task.state & GPM_SERVER_SEND_ELEMENT_VALUES) != 0) &&
                ((task.state & GPM_SERVER_WAIT_ELEMENT_VALUES) == 0)) {
                if (rComm.test((double*)&(task.elementValues[0])) == true) {
                    rComm.wait((double*)&(task.elementValues[0]));
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] SRV [OK]: send elements to [%d].", rComm.getRank(), proc)
                                  << std::endl;
                    }
                    task.state |= GPM_SERVER_WAIT_ELEMENT_VALUES;
                    continue;
                }
            }
        }

        // 終了したtaskをリストから除く
        // GpmServerTasksType::iterator itEnd = this->gpmServerTasks_.end();
        for (GpmServerTasksType::iterator it = this->gpmServerTasks_.begin(); it != itEnd;) {
            if ((it->state & GPM_SERVER_WAIT_ELEMENT_VALUES) != 0) {
                it = this->gpmServerTasks_.erase(it);
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

void TlDenseGeneralMatrix_ImplScalapack::getSparseMatrix_registerTask(TlSparseMatrix* pMatrix) const {
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
            task.components.resize(numOfProcs, std::vector<TlMatrixObject::index_type>());
            task.numOfComponents.resize(numOfProcs);
            task.elementValues.resize(numOfProcs);
            task.isFinished = false;

            // 要求された要素をどのプロセスが持っているかをチェック
            for (TlSparseMatrix::iterator it = pMatrix->begin(); it != pMatrix->end(); ++it) {
                TlMatrixObject::index_type globalRow = it->first.row;
                TlMatrixObject::index_type globalCol = it->first.col;

                int rank = 0;
                TlMatrixObject::index_type myRow = 0;
                TlMatrixObject::index_type myCol = 0;
                this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);
                // const TlMatrixObject::size_type localIndex =
                // this->getIndex(row, col);
                if (rank == this->rank_) {
                    // 自分が持っていれば格納すべき行列に値を代入する
                    pMatrix->set(globalRow, globalCol, this->getLocal(globalRow, globalCol));
                } else {
                    // 他プロセスが持っていれば、そのプロセスに要求するためのリストを作成する
                    // const int proc = this->getProcIdForIndex(row, col);
                    //                     std::cerr << TlUtils::format("[%d]
                    //                     check elements([%d] has (%d, %d)",
                    //                                                   rComm.getRank(),
                    //                                                   proc,
                    //                                                   row,
                    //                                                   col)
                    //                                << std::endl;
                    task.components[rank].push_back(globalRow);
                    task.components[rank].push_back(globalCol);
                }
            }
#pragma omp critical(TlDenseGeneralMatrix_Scalapack_gpmClientTasks_)
            { this->gpmClientTasks_[pMatrix] = task; }
        }
    }
}

int TlDenseGeneralMatrix_ImplScalapack::getPartialMatrix_getSessionId(const int proc) const {
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

bool TlDenseGeneralMatrix_ImplScalapack::getPartialMatrix_ClientTasks(TlMatrixObject* pMatrix) const {
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
                    task.state[proc] =
                        (GPM_CLIENT_COUNT_COMPONENTS | GPM_CLIENT_SEND_SESSION_ID | GPM_CLIENT_WAIT_SESSION_ID |
                         GPM_CLIENT_SEND_NUM_OF_COMPONENTS | GPM_CLIENT_WAIT_NUM_OF_COMPONENTS |
                         GPM_CLIENT_SEND_COMPONENTS | GPM_CLIENT_WAIT_COMPONENTS | GPM_CLIENT_RECV_ELEMENT_VALUES |
                         GPM_CLIENT_WAIT_ELEMENT_VALUES);
                }
                continue;
            }

            if ((task.state[proc] & GPM_CLIENT_SEND_SESSION_ID) == 0) {
                if (this->trafficControl_[proc] == NULL) {
                    // 相手にセッション番号を送信する
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format(
                                         "[%d] CLI [**]: send session id to "
                                         "[%d]. num = %d",
                                         rComm.getRank(), proc, task.numOfComponents[proc])
                                  << std::endl;
                    }

                    rComm.iSendData(task.sessionIds[proc], proc, TAG_GPM_SESSION_ID);
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
                            std::cerr << TlUtils::format(
                                             "[%d] CLI [OK]: send session id "
                                             "to [%d].",
                                             rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif  // NDEBUG

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
                        std::cerr << TlUtils::format(
                                         "[%d] CLI [**]: send num of "
                                         "components to [%d]. "
                                         "num = %d",
                                         rComm.getRank(), proc, task.numOfComponents[proc])
                                  << std::endl;
                    }
#endif  // NDEBUG

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
                            std::cerr << TlUtils::format(
                                             "[%d] CLI [OK]: send num of "
                                             "components to [%d].",
                                             rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif  // NDEBUG

                        task.state[proc] |= GPM_CLIENT_WAIT_NUM_OF_COMPONENTS;
                        continue;
                    }
                }
            }

            if ((task.state[proc] & GPM_CLIENT_SEND_COMPONENTS) == 0) {
                if (this->trafficControl_[proc] == pMatrix) {
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: send components to [%d].", rComm.getRank(), proc)
                                  << std::endl;
                    }

                    // 相手に要求要素組を送信する
                    rComm.iSendDataX((int*)&(task.components[proc][0]), task.numOfComponents[proc], proc,
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
                            std::cerr << TlUtils::format(
                                             "[%d] CLI [OK]: send components "
                                             "to [%d].",
                                             rComm.getRank(), proc)
                                      << std::endl;
                        }
#endif  // NDEBUG

                        task.state[proc] |= GPM_CLIENT_WAIT_COMPONENTS;
                        continue;
                    }
                }
            }

            if ((task.state[proc] & GPM_CLIENT_RECV_ELEMENT_VALUES) == 0) {
                if (this->trafficControl_[proc] == pMatrix) {
#ifndef NDEBUG
                    if (this->isDebugOut_GPM_ == true) {
                        std::cerr << TlUtils::format("[%d] CLI [**]: recv elements from [%d].", rComm.getRank(), proc)
                                  << std::endl;
                    }
#endif  // NDEBUG

                    // 行列データを待ち受ける
                    const std::size_t size = task.numOfComponents[proc] / 2;
                    task.elementValues[proc].resize(size);
                    rComm.iReceiveDataX((double*)&(task.elementValues[proc][0]), size, proc,
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
                            std::cerr << TlUtils::format(
                                             "[%d] CLI [OK]: recv elements "
                                             "from [%d].",
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
                    const TlMatrixObject::index_type row = task.components[proc][i * 2];
                    const TlMatrixObject::index_type col = task.components[proc][i * 2 + 1];
                    const double value = task.elementValues[proc][i];

                    pMatrix->set(row, col, value);
                    //                     if (this->isDebugOut_GPM_ == true) {
                    //                         std::cerr <<
                    //                         TlUtils::format("[%d] CLI [OK]:
                    //                         set (%d, %d) = %f.",
                    //                                                      rComm.getRank(),
                    //                                                      row,
                    //                                                      col,
                    //                                                      value)
                    //                                   << std::endl;
                    //                     }
                }
            }

#ifndef NDEBUG
            if (this->isDebugOut_GPM_ == true) {
                std::cerr << TlUtils::format("[%d] CLI [##]: task finished.", rComm.getRank()) << std::endl;
            }
#endif  // NDEBUG

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
#pragma omp critical(TlDenseGeneralMatrix_blacs_gpmClientTasks_)
                { this->gpmClientTasks_.erase(it); }

                if (this->isDebugOut_GPM_ == true) {
                    std::cerr << TlUtils::format("[%d] CLI [##]: task deleted.", rComm.getRank()) << std::endl;
                }

                answer = true;
            }
        }
    }

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::mergeMatrix_common(
    const std::vector<std::vector<TlDenseVectorObject::index_type> >& indexArrays,
    const std::vector<std::vector<double> >& values) {
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
            rComm.iSendDataX((TlMatrixObject::index_type*)&(indexArrays[proc][0]), size, proc,
                             TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES);
            rComm.iSendDataX((double*)&(values[proc][0]), size / 2, proc, TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES);
        }
    }

    // 行列データを送信する処理 ------------------------------------------------
    // 送られるリストの処理
    std::vector<std::size_t> numOfInputIndeces(numOfProcs, 0);
    std::vector<std::vector<TlMatrixObject::index_type> > inputIndexList(numOfProcs);
    std::vector<std::vector<double> > inputValues(numOfProcs);
    enum {
        RECV_NUM_OF_INDECES = 1,
        WAIT_NUM_OF_INDECES = 2,
        RECV_INDECES = 4,
        WAIT_INDECES = 8,
        RECV_VALUES = 16,
        WAIT_VALUES = 32,
        FINISHED = 64
    };
    std::vector<unsigned int> state(numOfProcs, 0);
    state[rank] = (RECV_NUM_OF_INDECES | WAIT_NUM_OF_INDECES | RECV_INDECES | WAIT_INDECES | RECV_VALUES | WAIT_VALUES |
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
                rComm.iReceiveData(numOfInputIndeces[proc], proc, TAG_TLDISTRIBUTE_MATRIX__MERGE__NUM_OF_INDECES);
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
                    rComm.iReceiveDataX((TlMatrixObject::index_type*)&(inputIndexList[proc][0]), size, proc,
                                        TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES);
                    state[proc] |= RECV_INDECES;
                } else {
                    assert(size == 0);
                    state[proc] |= (RECV_INDECES | WAIT_INDECES | RECV_VALUES | WAIT_VALUES | FINISHED);
                }
            }

            if (((state[proc] & RECV_INDECES) == RECV_INDECES) && ((state[proc] & WAIT_INDECES) != WAIT_INDECES)) {
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
                    rComm.iReceiveDataX((double*)&(inputValues[proc][0]), (size / 2), proc,
                                        TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES);
                    state[proc] |= RECV_VALUES;
                }
            }

            if (((state[proc] & RECV_VALUES) == RECV_VALUES) && ((state[proc] & WAIT_VALUES) != WAIT_VALUES)) {
                if (rComm.test(&(inputValues[proc][0])) == true) {
                    rComm.wait(&(inputValues[proc][0]));
                    state[proc] |= WAIT_VALUES;
                }
            }

            if (((state[proc] & WAIT_VALUES) == WAIT_VALUES) && ((state[proc] & WAIT_INDECES) == WAIT_INDECES)) {
                const std::size_t size = numOfInputIndeces[proc];
                const int max_i = size / 2;
                for (int i = 0; i < max_i; ++i) {
                    const TlMatrixObject::index_type globalRow = inputIndexList[proc][i * 2];
                    const TlMatrixObject::index_type globalCol = inputIndexList[proc][i * 2 + 1];
                    const double value = inputValues[proc][i];

                    int rank = 0;
                    TlMatrixObject::index_type myRow = 0;
                    TlMatrixObject::index_type myCol = 0;
                    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow, &myCol);

                    assert(rank == this->rank_);
                    TlMatrixObject::size_type index = this->getLocalIndex(myRow, myCol);
                    assert(index < this->getNumOfMyElements());
                    this->pData_[index] += value;
                }

                state[proc] |= FINISHED;
                ++terminateProc;
            }
        }

        if (terminateProc >= numOfProcs - 1) {
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

// -----------------------------------------------------------------------------
// std::vector<TlMatrixObject::MatrixElement>
// TlDenseGeneralMatrix_ImplScalapack::getMatrixElementsInLocal() const {
//   const TlMatrixObject::size_type numOfMyElements =
//   this->getNumOfMyElements(); std::vector<TlMatrixObject::MatrixElement>
//   answer(numOfMyElements);
//
//   TlMatrixObject::size_type count = 0;
//   const TlMatrixObject::index_type numOfRowIndeces =
//       this->m_RowIndexTable.size();
//   const TlMatrixObject::index_type numOfColIndeces =
//       this->m_ColIndexTable.size();
//   for (TlMatrixObject::index_type i = 0; i < numOfRowIndeces; ++i) {
//     const TlMatrixObject::index_type globalRow = this->m_RowIndexTable[i];
//
//     for (TlMatrixObject::index_type j = 0; j < numOfColIndeces; ++j) {
//       const TlMatrixObject::index_type globalCol = this->m_ColIndexTable[j];
//
//       answer[count] = TlMatrixObject::MatrixElement(
//           globalRow, globalCol, this->getLocal(globalRow, globalCol));
//       ++count;
//     }
//   }
//   assert(numOfMyElements == count);
//
//   return answer;
// }

std::vector<TlMatrixObject::MatrixElement> TlDenseGeneralMatrix_ImplScalapack::getMatrixElementsInLocal2() const {
    const TlMatrixObject::size_type numOfMyElements = this->getNumOfMyElements();
    std::vector<TlMatrixObject::MatrixElement> answer(numOfMyElements);

    const TlMatrixObject::index_type localRows = this->myRows_;
    const TlMatrixObject::index_type localCols = this->myCols_;
    TlMatrixObject::size_type count = 0;
    for (TlMatrixObject::index_type localRow = 0; localRow < localRows; ++localRow) {
        const TlMatrixObject::index_type globalRow = this->getLocal2Global_row(localRow);
        for (TlMatrixObject::index_type localCol = 0; localCol < localCols; ++localCol) {
            const TlMatrixObject::index_type globalCol = this->getLocal2Global_col(localCol);

            answer[count] = TlMatrixObject::MatrixElement(globalRow, globalCol, this->getLocal(globalRow, globalCol));
            ++count;
        }
    }

    assert(numOfMyElements == count);

    return answer;
}

void TlDenseGeneralMatrix_ImplScalapack::saveElements(
    TlDenseMatrix_IO_object* pFileMatrix, const std::vector<TlMatrixObject::MatrixElement>& elements) const {
    std::vector<TlMatrixObject::MatrixElement>::const_iterator itEnd = elements.end();
    for (std::vector<TlMatrixObject::MatrixElement>::const_iterator it = elements.begin(); it != itEnd; ++it) {
        pFileMatrix->set(it->row, it->col, it->value);
    }
}

// ---------------------------------------------------------------------------
// friends
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_ImplScalapack operator*(const TlDenseGeneralMatrix_ImplScalapack& X,
                                             const TlDenseGeneralMatrix_ImplScalapack& Y) {
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    // TlDenseGeneralMatrix_ImplScalapack& Xtmp =
    // const_cast<TlDenseGeneralMatrix_ImplScalapack&>(X);
    // TlDenseGeneralMatrix_ImplScalapack& Ytmp =
    // const_cast<TlDenseGeneralMatrix_ImplScalapack&>(Y);

    const int nXRow = X.getNumOfRows();
    const int nXCol = X.getNumOfCols();
    const int nYCol = Y.getNumOfCols();
    assert(X.getNumOfCols() == Y.getNumOfRows());

    const int nZRow = nXRow;
    const int nZCol = nYCol;
    TlDenseGeneralMatrix_ImplScalapack Z(nZRow, nZCol);

    const double dAlpha = 1.0;
    const double dBeta = 1.0;
    const int nIX = 1;
    const int nJX = 1;
    const int nIY = 1;
    const int nJY = 1;
    const int nIZ = 1;

    pdgemm_("N", "N", &nZRow, &nZCol, &nXCol, &dAlpha, X.pData_, &nIX, &nJX, X.pDESC_, Y.pData_, &nIY, &nJY, Y.pDESC_,
            &dBeta, Z.pData_, &nIZ, &nIZ, Z.pDESC_);

    return Z;
}

TlDenseVector_ImplScalapack operator*(const TlDenseGeneralMatrix_ImplScalapack& A,
                                      const TlDenseVector_ImplScalapack& X) {
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    // TlDenseGeneralMatrix_ImplScalapack& Atmp =
    // const_cast<TlDenseGeneralMatrix_ImplScalapack&>(A);
    TlDenseVector_ImplScalapack& Xtmp = const_cast<TlDenseVector_ImplScalapack&>(X);

    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();

    assert(N == X.getSize());

    TlDenseVector_ImplScalapack Y(M);

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

    pdgemv_(&TRANS, &M, &N, &alpha, A.pData_, &IA, &JA, A.pDESC_, Xtmp.pData_, &IX, &JX, Xtmp.pDESC_, &INCX, &beta,
            Y.pData_, &IY, &JY, Y.pDESC_, &INCY);

    return Y;
}

TlDenseVector_ImplScalapack operator*(const TlDenseVector_ImplScalapack& X,
                                      const TlDenseGeneralMatrix_ImplScalapack& A) {
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    // TlDenseGeneralMatrix_ImplScalapack& Atmp =
    // const_cast<TlDenseGeneralMatrix_ImplScalapack&>(A);
    TlDenseVector_ImplScalapack& Xtmp = const_cast<TlDenseVector_ImplScalapack&>(X);

    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();

    assert(N == X.getSize());

    TlDenseVector_ImplScalapack Y(M);

    const char TRANS = 'T';
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

    pdgemv_(&TRANS, &M, &N, &alpha, A.pData_, &IA, &JA, A.pDESC_, Xtmp.pData_, &IX, &JX, Xtmp.pDESC_, &INCX, &beta,
            Y.pData_, &IY, &JY, Y.pDESC_, &INCY);

    return Y;
}

// ---------------------------------------------------------------------------
