#include <cassert>
#include <cstdlib>
//#include <bitset>
//#include <list>
//#include <vector>

#include "scalapack.h"
#include "tl_dense_scalapack_object.h"
#include "tl_scalapack_context.h"

const std::size_t TlDenseScalapackObject::FILE_BUFFER_SIZE =
    100 * 1024 * 1024;  // 100 MB

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseScalapackObject::TlDenseScalapackObject(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col)
    : context_(0),
      rows_(row),
      cols_(col),
      blockSize_(TlScalapackContext::getBlockSize()),
      pData_(NULL) {
    this->initialize();
}

TlDenseScalapackObject::TlDenseScalapackObject(
    const TlDenseScalapackObject& rhs)
    : context_(0),
      rows_(rhs.rows_),
      cols_(rhs.cols_),
      blockSize_(rhs.blockSize_),
      pData_(NULL) {
    this->initialize();
    if (rhs.getNumOfMyElements() > 0) {
        std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(),
                  this->pData_);
    }
}

TlDenseScalapackObject::~TlDenseScalapackObject() {
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseScalapackObject::getNumOfRows() const {
    return this->rows_;
}

TlMatrixObject::index_type TlDenseScalapackObject::getNumOfCols() const {
    return this->cols_;
}

TlMatrixObject::size_type TlDenseScalapackObject::getNumOfElements() const {
    const TlMatrixObject::size_type rows = this->getNumOfRows();
    const TlMatrixObject::size_type cols = this->getNumOfCols();
    return rows * cols;
}

void TlDenseScalapackObject::resize(const TlMatrixObject::index_type row,
                                    const TlMatrixObject::index_type col) {
    assert(row > 0);
    assert(col > 0);
    if ((row == this->getNumOfRows()) && (col == this->getNumOfCols())) {
        // do not need operation.
        return;
    }

    // backup old condition
    TlDenseScalapackObject tmp(*this);

    // new size and zero clear
    this->rows_ = row;
    this->cols_ = col;
    this->initialize();

    // copy
    const TlMatrixObject::index_type maxLocalRow =
        std::min(this->getNumOfMyRows(), tmp.getNumOfMyRows());
    const TlMatrixObject::index_type maxLocalCol =
        std::min(this->getNumOfMyCols(), tmp.getNumOfMyCols());
    for (TlMatrixObject::index_type localRow = 0; localRow < maxLocalRow;
         ++localRow) {
        for (TlMatrixObject::index_type localCol = 0; localCol < maxLocalCol;
             ++localCol) {
            const TlMatrixObject::size_type localIndex =
                this->getLocalIndex(localRow, localCol);
            this->pData_[localIndex] = tmp.pData_[localIndex];
        }
    }
}

double TlDenseScalapackObject::get(
    const TlMatrixObject::index_type globalRow,
    const TlMatrixObject::index_type globalCol) const {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    // DataType& data_tmp = const_cast<DataType&>(this->data_);

    double answer = 0.0;
    const int fortranGlobalRow = globalRow + 1;
    const int fortranGlobalCol = globalCol + 1;
    pdelget_("A", " ", &answer, this->pData_, &fortranGlobalRow,
             &fortranGlobalCol, this->pDESC_);

    return answer;
}

void TlDenseScalapackObject::set(const TlMatrixObject::index_type globalRow,
                                 const TlMatrixObject::index_type globalCol,
                                 const double value) {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow,
                                      &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index =
            this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        this->pData_[index] = value;
    }
}

void TlDenseScalapackObject::add(const TlMatrixObject::index_type globalRow,
                                 const TlMatrixObject::index_type globalCol,
                                 const double value) {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow,
                                      &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index =
            this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        this->pData_[index] += value;
    }
}

void TlDenseScalapackObject::mul(const TlMatrixObject::index_type globalRow,
                                 const TlMatrixObject::index_type globalCol,
                                 const double value) {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow,
                                      &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index =
            this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        this->pData_[index] *= value;
    }
}

double TlDenseScalapackObject::getLocal(
    const TlMatrixObject::index_type globalRow,
    const TlMatrixObject::index_type globalCol) const {
    assert((0 <= globalRow) && (globalRow < this->getNumOfRows()));
    assert((0 <= globalCol) && (globalCol < this->getNumOfCols()));

    double answer = 0.0;

    int rank = 0;
    TlMatrixObject::index_type myRow = 0;
    TlMatrixObject::index_type myCol = 0;
    this->getGlobalRowCol2LocalRowCol(globalRow, globalCol, &rank, &myRow,
                                      &myCol);
    if (rank == this->rank_) {
        const TlMatrixObject::size_type index =
            this->getLocalIndex(myRow, myCol);
        assert(0 <= index && index < this->getNumOfMyElements());
        answer = this->pData_[index];
    }

    return answer;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseScalapackObject& TlDenseScalapackObject::operator=(
    const TlDenseScalapackObject& rhs) {
    if (&rhs != this) {
        assert(rhs.getNumOfRows() > 0);
        assert(rhs.getNumOfCols() > 0);
        this->rows_ = rhs.rows_;
        this->cols_ = rhs.cols_;
        this->blockSize_ = rhs.blockSize_;

        this->initialize();
        std::copy(rhs.pData_, rhs.pData_ + rhs.getNumOfMyElements(),
                  this->pData_);
    }

    return (*this);
}

TlDenseScalapackObject& TlDenseScalapackObject::operator+=(
    const TlDenseScalapackObject& rhs) {
    assert(this->rows_ == rhs.rows_);
    assert(this->cols_ == rhs.cols_);
    assert(this->blockSize_ == rhs.blockSize_);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] += rhs.pData_[i];
    }

    return (*this);
}

TlDenseScalapackObject& TlDenseScalapackObject::operator-=(
    const TlDenseScalapackObject& rhs) {
    assert(this->rows_ == rhs.rows_);
    assert(this->cols_ == rhs.cols_);
    assert(this->blockSize_ == rhs.blockSize_);

    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] -= rhs.pData_[i];
    }

    return (*this);
}

TlDenseScalapackObject& TlDenseScalapackObject::operator*=(const double coef) {
    const std::size_t bufSize = this->getNumOfMyElements();
    for (std::size_t i = 0; i < bufSize; ++i) {
        this->pData_[i] *= coef;
    }

    return (*this);
}

TlDenseScalapackObject& TlDenseScalapackObject::operator/=(const double coef) {
    return this->operator*=(1.0 / coef);
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
void TlDenseScalapackObject::initialize() {
    TlScalapackContext::getData(this->context_, this->proc_, this->rank_,
                                this->procGridRow_, this->procGridCol_);

    // my process position on the process matrix
    Cblacs_gridinfo(this->context_, &(this->procGridRow_),
                    &(this->procGridCol_), &(this->m_nMyProcRow),
                    &(this->m_nMyProcCol));

    // determine sizes of local matrix
    const int nStartRowProc = 0;
    const int nStartColProc = 0;
    this->myRows_ = std::max(
        1, numroc_(&(this->rows_), &(this->blockSize_), &(this->m_nMyProcRow),
                   &nStartRowProc, &(this->procGridRow_)));
    this->myCols_ = std::max(
        1, numroc_(&(this->cols_), &(this->blockSize_), &(this->m_nMyProcCol),
                   &nStartColProc, &(this->procGridCol_)));

    // make parameter, desca
    int nInfo = 0;
    descinit_(this->pDESC_, &(this->rows_), &(this->cols_), &(this->blockSize_),
              &(this->blockSize_), &nStartRowProc, &nStartColProc,
              &(this->context_), &(this->myRows_), &nInfo);
    assert(nInfo == 0);

    // 行列データ用バッファメモリの確保
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
    this->pData_ = new double[this->getNumOfMyElements()];
    std::fill(this->pData_, this->pData_ + this->getNumOfMyElements(), 0.0);
}

TlMatrixObject::index_type TlDenseScalapackObject::getNumOfMyRows() const {
    return this->myRows_;
}

TlMatrixObject::index_type TlDenseScalapackObject::getNumOfMyCols() const {
    return this->myCols_;
}

TlMatrixObject::size_type TlDenseScalapackObject::getNumOfMyElements() const {
    const TlMatrixObject::size_type rows = this->getNumOfMyRows();
    const TlMatrixObject::size_type cols = this->getNumOfMyCols();

    return rows * cols;
}

void TlDenseScalapackObject::getGlobalRowCol2LocalRowCol(
    const TlMatrixObject::index_type globalRow,
    const TlMatrixObject::index_type globalCol, int* pRank,
    TlMatrixObject::index_type* pLocalRow,
    TlMatrixObject::index_type* pLocalCol) const {
    const std::div_t rowBlockDiv = std::div(globalRow, this->blockSize_);
    const std::div_t colBlockDiv = std::div(globalCol, this->blockSize_);
    const int rowBlockIndex = rowBlockDiv.quot;
    const int colBlockIndex = colBlockDiv.quot;
    const int rowLocalIndex = rowBlockDiv.rem;
    const int colLocalIndex = colBlockDiv.rem;

    const std::div_t rowIndexResult =
        std::div(rowBlockIndex, this->procGridRow_);
    const std::div_t colIndexResult =
        std::div(colBlockIndex, this->procGridCol_);
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

TlMatrixObject::index_type TlDenseScalapackObject::getLocal2Global_row(
    const TlMatrixObject::index_type localIndex) const {
    return this->getLocal2Global(localIndex, this->m_nMyProcRow,
                                 this->procGridRow_);
}

TlMatrixObject::index_type TlDenseScalapackObject::getLocal2Global_col(
    const TlMatrixObject::index_type localIndex) const {
    return this->getLocal2Global(localIndex, this->m_nMyProcCol,
                                 this->procGridCol_);
}

TlMatrixObject::index_type TlDenseScalapackObject::getLocal2Global(
    const TlMatrixObject::index_type localIndex, const int procID,
    const int numOfProcs) const {
    const int blockSize = TlScalapackContext::getBlockSize();

    const int offset = procID * blockSize;
    const int distance = numOfProcs * blockSize;

    const std::div_t div_result = std::div(localIndex, blockSize);
    const int globalIndex =
        offset + distance * div_result.quot + div_result.rem;

    return globalIndex;
}

TlMatrixObject::index_type TlDenseScalapackObject::getLocalIndex(
    const TlMatrixObject::index_type localRow,
    const TlMatrixObject::index_type localCol) const {
    const TlMatrixObject::size_type index = localRow + this->myRows_ * localCol;
    return index;
}

std::vector<TlMatrixObject::MatrixElement>
TlDenseScalapackObject::getMatrixElementsInLocal() const {
    const TlMatrixObject::size_type numOfMyElements =
        this->getNumOfMyElements();
    std::vector<TlMatrixObject::MatrixElement> answer(numOfMyElements);

    const TlMatrixObject::index_type localRows = this->myRows_;
    const TlMatrixObject::index_type localCols = this->myCols_;
    TlMatrixObject::size_type count = 0;
    for (TlMatrixObject::index_type localRow = 0; localRow < localRows;
         ++localRow) {
        const TlMatrixObject::index_type globalRow =
            this->getLocal2Global_row(localRow);
        for (TlMatrixObject::index_type localCol = 0; localCol < localCols;
             ++localCol) {
            const TlMatrixObject::index_type globalCol =
                this->getLocal2Global_col(localCol);

            answer[count] = TlMatrixObject::MatrixElement(
                globalRow, globalCol, this->getLocal(globalRow, globalCol));
            ++count;
        }
    }

    assert(numOfMyElements == count);

    return answer;
}
