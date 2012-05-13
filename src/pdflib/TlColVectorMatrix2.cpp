#include <iostream>
#include "TlColVectorMatrix2.h"


TlColVectorMatrix2::TlColVectorMatrix2(const index_type row,
                                       const index_type col,
                                       int allProcs, int rank)
    : allProcs_(allProcs), rank_(rank) {
    this->resize(row, col);
}


TlColVectorMatrix2::TlColVectorMatrix2(const TlColVectorMatrix2& rhs) 
    : allProcs_(rhs.allProcs_), rank_(rhs.rank_) {
    this->resize(rhs.getNumOfRows(), rhs.getNumOfCols());
    this->data_ = rhs.data_;
}

TlColVectorMatrix2::~TlColVectorMatrix2()
{
}


TlColVectorMatrix2& TlColVectorMatrix2::operator=(const TlColVectorMatrix2& rhs)
{
    if (this != &rhs) {
        this->allProcs_ = rhs.allProcs_;
        this->rank_ = rhs.rank_;
        this->resize(rhs.getNumOfRows(), rhs.getNumOfCols());
        this->data_ = rhs.data_;
    }
}


void TlColVectorMatrix2::resize(const index_type newRows,
                                const index_type newCols)
{
    this->numOfRows_ = newRows;
    this->numOfCols_ = newCols;

    const div_t turns = std::div(newCols, this->allProcs_);
    const index_type localCols = turns.quot + 1;
    this->data_.resize(localCols);

    if (this->reserveRows_ < newRows) {
        this->reserveRows_ = newRows;
    }
    for (index_type i = 0; i < localCols; ++i) {
        this->data_[i].reserve(this->reserveRows_);
        this->data_[i].resize(newRows);
    }
}


void TlColVectorMatrix2::reserve_rows(const index_type newReserves) {
    this->reserveRows_ = std::max(this->getNumOfRows(), newReserves);
}


void TlColVectorMatrix2::set(index_type row, index_type col, double value)
{
    const div_t turns = std::div(col, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type col = turns.quot;
        this->data_[col][row] = value;
    }
}


TlVector TlColVectorMatrix2::getColVector(index_type col) const
{
    TlVector answer;
    const div_t turns = std::div(col, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type col = turns.quot;
        answer = TlVector(this->data_[col]);
    }

    return answer;
}


TlMatrixObject::index_type 
TlColVectorMatrix2::getColVector(index_type col,
                                 double *pBuf,
                                 index_type maxRowSize) const
{
    index_type copySize = 0;
    const div_t turns = std::div(col, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type col = turns.quot;
        copySize = std::min(this->getNumOfCols(), maxRowSize);
        std::copy(this->data_[col].begin(),
                  this->data_[col].begin() + copySize,
                  pBuf);
    }

    return copySize;
}


int TlColVectorMatrix2::getPEinChargeByCol(const index_type col) const
{
    assert((0 <= col) && (col < this->getNumOfCols()));
    const div_t turns = std::div(col, this->allProcs_);
    return turns.rem;
}


void TlColVectorMatrix2::save(const std::string& basename) const 
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();

    std::ofstream ofs;
    const std::string path = TlUtils::format("%s.part%d.mat",
                                             basename.c_str(),
                                             this->rank_);
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    ofs.write(reinterpret_cast<const char*>(&numOfRows), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&numOfCols), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&(this->allProcs_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->rank_)), sizeof(int));

    // data
    const index_type numOfLocalCols = this->data_.size();
    for (int i = 0; i < numOfLocalCols; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i][0])), sizeof(double) * numOfRows);
    }
    
    ofs.close();
}


void TlColVectorMatrix2::load(const std::string& basename)
{
    std::ifstream ifs;
    const std::string path = TlUtils::format("%s.paart%d.mat",
                                             basename.c_str(),
                                             this->rank_);
    ifs.open(path.c_str());

    // header
    index_type numOfRows = 0;
    index_type numOfCols = 0;
    int allProcs = 0;
    int rank = 0;
    ifs.read((char*)&numOfRows, sizeof(index_type));
    ifs.read((char*)&numOfCols, sizeof(index_type));
    ifs.read((char*)&allProcs, sizeof(int));
    ifs.read((char*)&rank, sizeof(int));
    this->resize(numOfRows, numOfCols);
    
    if ((allProcs != this->allProcs_) || 
        (rank != this->rank_)) {
        std::cerr << "something wrong!" << std::endl;
    }

    // data
    const div_t turns = std::div(this->getNumOfRows(), this->allProcs_);
    index_type numOfLocalCols = turns.quot;
    if (this->rank_ < turns.rem) {
        ++numOfLocalCols;
    }

    for (int i = 0; i < numOfLocalCols; ++i) {
        ifs.read((char*)&(this->data_[i][0]), sizeof(double) * this->getNumOfRows());
    }

    ifs.close();
}


