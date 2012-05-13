#include <iostream>
#include "TlRowVectorMatrix2.h"
#include "TlCommunicate.h"

TlRowVectorMatrix2::TlRowVectorMatrix2(const index_type row,
                                       const index_type col,
                                       int allProcs, int rank)
    : allProcs_(allProcs), rank_(rank) {
    this->resize(row, col);
}


TlRowVectorMatrix2::TlRowVectorMatrix2(const TlRowVectorMatrix2& rhs) 
    : allProcs_(rhs.allProcs_), rank_(rhs.rank_) {
    this->resize(rhs.getNumOfRows(), rhs.getNumOfCols());
    this->data_ = rhs.data_;
}


TlRowVectorMatrix2::~TlRowVectorMatrix2()
{
}


TlRowVectorMatrix2& TlRowVectorMatrix2::operator=(const TlRowVectorMatrix2& rhs)
{
    if (this != &rhs) {
        this->allProcs_ = rhs.allProcs_;
        this->rank_ = rhs.rank_;
        this->resize(rhs.getNumOfRows(), rhs.getNumOfCols());
        this->data_ = rhs.data_;
    }

    return *this;
}


void TlRowVectorMatrix2::resize(const index_type newRows,
                                const index_type newCols)
{
    this->numOfRows_ = newRows;
    this->numOfCols_ = newCols;

    const div_t turns = std::div(newRows, this->allProcs_);
    const index_type localRows = turns.quot + 1;
    this->data_.resize(localRows);

    if (this->reserveCols_ < newCols) {
        this->reserveCols_ = newCols;
    }
    for (index_type i = 0; i < localRows; ++i) {
        this->data_[i].reserve(this->reserveCols_);
        this->data_[i].resize(newCols);
    }
}


void TlRowVectorMatrix2::reserve_cols(const index_type newReserves) {
    this->reserveCols_ = std::max(this->getNumOfCols(), newReserves);
}


void TlRowVectorMatrix2::set(index_type row, index_type col, double value)
{
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        this->data_[row][col] = value;
    }
}


TlVector TlRowVectorMatrix2::getRowVector(index_type row) const
{
    TlVector answer;
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        answer = TlVector(this->data_[row]);
    }

    return answer;
}


TlMatrixObject::index_type 
TlRowVectorMatrix2::getRowVector(index_type row,
                                double *pBuf,
                                index_type maxColSize) const
{
    index_type copySize = 0;
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        copySize = std::min(this->getNumOfCols(), maxColSize);
        std::copy(this->data_[row].begin(),
                  this->data_[row].begin() + copySize,
                  pBuf);
    }

    return copySize;
}


int TlRowVectorMatrix2::getPEinChargeByRow(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    const div_t turns = std::div(row, this->allProcs_);
    return turns.rem;
}


TlMatrix TlRowVectorMatrix2::getTlMatrix() const 
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfRows, this->allProcs_);
    const index_type localRows = turns.quot + ((this->rank_ < turns.rem) ? 1 : 0);
    for (index_type r = 0; r < localRows; ++r) {
        const index_type row = r * this->allProcs_ + this->rank_;
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, this->data_[r][col]);
        }
    }
    
    return answer;
}


// TlColVectorMatrix2 TlRowVectorMatrix2::getColVectorMatrix() const 
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     const index_type numOfRows = this->getNumOfRows();
//     const index_type numOfCols = this->getNumOfCols();
//     const int numOfProcs = this->allProcs_;
//     const int myRank = this->rank_;
//     TlColVectorMatrix2 answer(numOfRows, numOfCols, numOfProcs, myRank);

//     const div_t turns = std::div(numOfRows, numOfProcs);
//     const index_type localRows = turns.quot + 1;
//     std::vector<double> buf(localRows * numOfCols);;
//     for (int i = 0; i < numOfProcs; ++i) {
//         if (i == myRank) {
//             for (index_type j = 0; j < localRows; ++j) {
//                 std::copy(this->data_[j].begin(),
//                           this->data_[j].begin() + numOfCols,
//                           buf.begin() + numOfCols * j);
//             }
//         }
//         rComm.broadcast(&(buf[0]), localRows * numOfCols, i);

//         // set
//         for (index_type j = 0; j < localRows; ++j) {
//             index_type row = numOfProcs * j + i;
//             if (row < numOfRows) {
//                 for (index_type col = 0; col < numOfCols; ++col) {
//                     answer.set(row, col, buf[numOfCols * j + col]);
//                 }
//             }
//         }
//     }
    
//     return answer;
// }


void TlRowVectorMatrix2::save(const std::string& basename) const 
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
    const index_type numOfLocalRows = this->data_.size();
    for (int i = 0; i < numOfLocalRows; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i][0])), sizeof(double) * numOfCols);
    }
    
    ofs.close();
}


void TlRowVectorMatrix2::load(const std::string& basename)
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
    
    assert((allProcs == this->allProcs_) && (rank == this->rank_));

    // data
    const div_t turns = std::div(this->getNumOfRows(), this->allProcs_);
    index_type numOfLocalRows = turns.quot;
    if (this->rank_ < turns.rem) {
        ++numOfLocalRows;
    }

    for (int i = 0; i < numOfLocalRows; ++i) {
        ifs.read((char*)&(this->data_[i][0]), sizeof(double) * this->getNumOfCols());
    }

    ifs.close();
}


