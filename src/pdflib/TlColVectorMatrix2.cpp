#include <iostream>
#include "TlColVectorMatrix2.h"
#include "TlMemManager.h"

TlColVectorMatrix2::TlColVectorMatrix2(const index_type row,
                                       const index_type col,
                                       int allProcs, int rank,
                                       bool isUsingMemManager)
    : numOfRows_(0), numOfCols_(0), reserveRows_(0),
      allProcs_(allProcs), rank_(rank),
      numOfLocalCols_(0),
      isUsingMemManager_(isUsingMemManager) {
    this->resize(row, col);
}


TlColVectorMatrix2::~TlColVectorMatrix2()
{
    const index_type numOfLocalCols = this->numOfLocalCols_;

    if (this->isUsingMemManager_ == true) {
        TlMemManager& rMemManager = TlMemManager::getInstance();
        const index_type reserveRows = this->reserveRows_;
        for (index_type i = 0; i < numOfLocalCols; ++i) {
            rMemManager.deallocate((char*)this->data_[i], sizeof(double)*reserveRows);
            this->data_[i] = NULL;
        }
    } else {
        for (index_type i = 0; i < numOfLocalCols; ++i) {
            delete[] this->data_[i];
            this->data_[i] = NULL;
        }
    }
    this->data_.clear();
}


void TlColVectorMatrix2::resize(const index_type newRows,
                                const index_type newCols)
{
    this->numOfCols_ = newCols;
    index_type prevNumOfLocalCols = this->numOfLocalCols_;
    const div_t turns = std::div(newCols, this->allProcs_);
    index_type newNumOfLocalCols = turns.quot;
    if (this->rank_ < turns.rem) {
        newNumOfLocalCols += 1;
    }
    this->numOfLocalCols_ = newNumOfLocalCols;
    if (newNumOfLocalCols > prevNumOfLocalCols) {
        this->data_.resize(newNumOfLocalCols, NULL);
    } else if (newNumOfLocalCols < prevNumOfLocalCols) {
        if (this->isUsingMemManager_ == true) {
            TlMemManager& rMemManager = TlMemManager::getInstance();
            const index_type reserveRows = this->reserveRows_;
            for (index_type i = newNumOfLocalCols; i < prevNumOfLocalCols; ++i) {
                rMemManager.deallocate((char*)this->data_[i], sizeof(double)*reserveRows);
                this->data_[i] = NULL;
            }
        } else {
            for (index_type i = newNumOfLocalCols; i < prevNumOfLocalCols; ++i) {
                delete this->data_[i];
                this->data_[i] = NULL;
            }
        }
        this->data_.resize(newNumOfLocalCols);
    }
    
    this->numOfRows_ = newRows;
    this->reserve_rows(newRows);
}


void TlColVectorMatrix2::reserve_rows(const index_type newReserves) {
    const index_type prevReserveRows = this->reserveRows_;
    const index_type newReserveRows = std::max(this->getNumOfRows(), newReserves);
    const index_type numOfRows = this->getNumOfRows();

    if (prevReserveRows < newReserveRows) {
        this->reserveRows_ = newReserveRows;
        const index_type numOfLocalCols = this->numOfLocalCols_;
        for (index_type i = 0; i < numOfLocalCols; ++i) {
            double* pNew = NULL;
            if (this->isUsingMemManager_ == true) {
                TlMemManager& rMemManager = TlMemManager::getInstance();
                pNew = (double*)rMemManager.allocate(sizeof(double)*newReserveRows);
            } else {
                pNew = new double[newReserveRows];
            }

            for (index_type j = 0; j < newReserveRows; ++j) {
                pNew[j] = 0.0;
            }
            
            if (this->data_[i] != NULL) {
                for (index_type j = 0; j < numOfRows; ++j) {
                    pNew[j] = this->data_[i][j];
                }
                if (this->isUsingMemManager_ == true) {
                    TlMemManager& rMemManager = TlMemManager::getInstance();
                    rMemManager.deallocate((char*)this->data_[i], sizeof(double)*prevReserveRows);
                } else {
                    delete[] this->data_[i];
                }
                this->data_[i] = NULL;
            }

            this->data_[i] = pNew;
        }
    }
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
        answer = TlVector(this->data_[col], this->getNumOfRows());
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
        copySize = std::min(this->getNumOfRows(), maxRowSize);
        std::copy(this->data_[col],
                  this->data_[col] + copySize,
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
    const std::string path = TlUtils::format("%s.part%d.mat",
                                             basename.c_str(),
                                             this->rank_);
    //std::cerr << "TlColVectorMatrix2::load() path=" << path << ", " << this->allProcs_ << std::endl;
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
    
    // std::cerr << TlUtils::format("ERROR: TlColVectorMatrix2::load() %s %d/%d (%dx%d; %d/%d)",
    //                              path.c_str(),
    //                              this->rank_, this->allProcs_,
    //                              numOfRows, numOfCols, rank, allProcs)
    //           << std::endl;
    assert((allProcs == this->allProcs_) && (rank == this->rank_));

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



TlMatrix TlColVectorMatrix2::getTlMatrix() const 
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    const index_type numOfLocalCols = this->numOfLocalCols_;
    for (index_type c = 0; c < numOfLocalCols; ++c) {
        const index_type col = c * this->allProcs_ + this->rank_;
        for (index_type row = 0; row < numOfRows; ++row) {
            answer.set(row, col, this->data_[c][row]);
        }
    }
    
    return answer;
}


