#include "TlCommunicate.h"
#include "TlRowVectorMatrix.h"


TlRowVectorMatrix::TlRowVectorMatrix(const index_type row,
                                     const index_type col)
{
    this->globalRows_ = 0;
    this->globalCols_ = 0;
    this->reserveCols_ = 10;
    this->resize(row, col);
}


TlRowVectorMatrix::~TlRowVectorMatrix()
{
}


void TlRowVectorMatrix::resize(const index_type new_globalRows,
                               const index_type new_globalCols)
{
    const index_type prev_globalRows = this->globalRows_;

    if (prev_globalRows != new_globalRows) {
        // remake PE table
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int numOfProcs = rComm.getNumOfProcs();
        const int myRank = rComm.getRank();

        this->row_PE_table_.resize(new_globalRows);
        int pe = 0;
        for (index_type r = 0; r < new_globalRows; ++r) {
            this->row_PE_table_[r] = pe;
            
            ++pe;
            if (pe >= numOfProcs) {
                pe = 0;
            }
        }

        std::vector<index_type> localRowTable;
        localRowTable.reserve(new_globalRows / numOfProcs +1);
        for (index_type r = 0; r < new_globalRows; ++r) {
            if (this->row_PE_table_[r] == myRank) {
                localRowTable.push_back(r);
            }
        }

        const std::size_t prev_localRowTableSize = this->data_.size();
        const std::size_t new_localRowTableSize = localRowTable.size();
        if (prev_localRowTableSize < new_localRowTableSize) {
            this->data_.resize(new_localRowTableSize);
            for (std::size_t i = prev_localRowTableSize; i < new_localRowTableSize; ++i) {
                this->data_[i].row = localRowTable[i];
            }
        } else if (prev_localRowTableSize < new_localRowTableSize) {
            this->data_.resize(new_localRowTableSize);
        }
    }

    std::vector<RowVector>::iterator itEnd = this->data_.end();
    for (std::vector<RowVector>::iterator it = this->data_.begin(); it != itEnd; ++it) {
        it->cols.reserve(this->reserveCols_);
        it->cols.resize(new_globalCols);
    }

    this->globalRows_ = new_globalRows;
    this->globalCols_ = new_globalCols;
}


void TlRowVectorMatrix::reserve_cols(const index_type cols) {
    this->reserveCols_ = std::max(globalCols_, cols);
}


void TlRowVectorMatrix::set(index_type row, index_type col, double value)
{
    std::vector<RowVector>::iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), RowVector(row));

    if ((it != this->data_.end()) && (it->row == row)) {
        it->cols[col] = value;
    }
}


TlVector TlRowVectorMatrix::getRowVector(index_type row) const
{
    TlVector answer;
    std::vector<RowVector>::const_iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), RowVector(row));
    if (it != this->data_.end()) {
        assert(it->row == row);
        answer = it->cols;
    }

    return answer;
}


TlMatrixObject::index_type 
TlRowVectorMatrix::getRowVector(index_type row,
                                double *pBuf,
                                index_type maxColSize) const
{
    index_type copySize = 0;
    std::vector<RowVector>::const_iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), RowVector(row));
    if (it != this->data_.end()) {
        assert(it->row == row);
        assert(it->cols.size() == this->getNumOfCols());
        copySize = std::min(this->getNumOfCols(), maxColSize);
        std::copy(it->cols.begin(), it->cols.begin() + copySize,
                  pBuf);
    }
    return copySize;
}


int TlRowVectorMatrix::getPEinChargeByRow(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    return this->row_PE_table_[row];
}


TlMatrix TlRowVectorMatrix::getTlMatrix() const 
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlMatrix answer(row, col);
    
    std::vector<RowVector>::const_iterator itEnd = this->data_.end();
    for (std::vector<RowVector>::const_iterator it = this->data_.begin(); it != itEnd; ++it) {
        const index_type r = it->row;
        for (index_type c = 0; c < col; ++c) {
            answer.set(r, c, it->cols[c]);
        }
    }

    rComm.allReduce_SUM(answer);
    return answer;
}


TlDistributeMatrix TlRowVectorMatrix::getTlDistributeMatrix() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlDistributeMatrix answer(row, col);

    int blocks = 0;
    std::vector<index_type> rows;
    std::vector<double> vtr;
    for (int pe = 0; pe < numOfProcs; ++pe) {
        if (pe == myRank) {
            blocks = this->data_.size();
            vtr.resize(col * blocks);
            for (int i = 0; i < blocks; ++i) {
                rows[i] = this->data_[i].row;
                for (index_type c = 0; c < col; ++c) {
                    vtr[col * i + c] = this->data_[i].cols[c];
                }
            }
        }

        rComm.broadcast(blocks, pe);
        rComm.broadcast(rows, pe);
        rComm.broadcast(vtr, pe);

        // set
        for (int i = 0; i < blocks; ++i) {
            const index_type r = rows[i];
            for (index_type c = 0; c < col; ++c) {
                answer.set(r, c, vtr[col * i + c]);
            }
        }
    }

    return answer;
}


void TlRowVectorMatrix::save(const std::string& basename,
                                          const int PE) const
{
    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();

    std::ofstream ofs;
    const std::string path = TlUtils::format("%s.%d.dm",
                                             basename.c_str(),
                                             PE);
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    int numOfMyRows = this->row_PE_table_.size();
    ofs.write(reinterpret_cast<const char*>(&PE), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&row), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&col), sizeof(index_type));

    // data
    ofs.write(reinterpret_cast<const char*>(&numOfMyRows), sizeof(int));
    for (int i = 0; i < numOfMyRows; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i].row)), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i].cols[0])), sizeof(double) * col);
    }
    
    ofs.close();
}


void TlRowVectorMatrix::load(const std::string& basename,
                                          const int PE)
{
    std::ifstream ifs;
    const std::string path = TlUtils::format("%s.%d.dm",
                                             basename.c_str(),
                                             PE);
    ifs.open(path.c_str());

    // header
    int inPE = 0;
    index_type row = 0;
    index_type col = 0;
    ifs.read((char*)&inPE, sizeof(int));
    ifs.read((char*)&row, sizeof(index_type));
    ifs.read((char*)&col, sizeof(index_type));
    this->resize(row, col);

    // data
    int numOfMyRows = 0;
    ifs.read((char*)&numOfMyRows, sizeof(numOfMyRows));
    assert(numOfMyRows == this->row_PE_table_.size());
    for (int i = 0; i < numOfMyRows; ++i) {
        index_type row_index = 0;
        ifs.read((char*)&row_index, sizeof(index_type));
        assert(this->data_[i].row == row_index);
        assert(this->data_[i].cols.size() == col);
        ifs.read((char*)&(this->data_[i].cols[0]), sizeof(double)*col);
    }

    ifs.close();
}


TlColVectorMatrix TlRowVectorMatrix::getTlColVectorMatrix() const 
{
    TlColVectorMatrix answer(this->getNumOfRows(), this->getNumOfCols());

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    const int rowCycles = numOfRows / numOfProcs;
    const int colCycles = numOfCols / numOfProcs;
    std::vector<double> buf((rowCycles +1) * numOfCols);
    for (int pe = 0; pe < numOfProcs; ++pe) {
        // calc comm number of data parts
        int numOfParts = rowCycles;
        if ((rowCycles * numOfProcs + pe) > numOfRows) {
            ++numOfParts;
        }

        // prepare comm data
        if (pe == myRank) {
            assert(numOfParts == this->data_.size());
            for (int i = 0; i < rowCycles; ++i) {
                std::copy(this->data_[i].cols.begin(),
                          this->data_[i].cols.begin() + numOfCols,
                          buf.begin() + numOfCols * i);
            }
        }
        
        rComm.broadcast(&(buf[0]), numOfCols * numOfParts, pe);

        // set to matrix
        for (int i = 0; i < numOfParts; ++i) {
            const index_type row = i * numOfProcs + pe;

            for (int j = 0; j < (colCycles -1); ++j) {
                const index_type col = j * numOfProcs + myRank;
                const double value = buf[col];
                answer.set(row, col, value);
            }
            {
                const index_type last_col = colCycles * numOfProcs + myRank;
                if (last_col < numOfCols) {
                    const double value = buf[last_col];
                    answer.set(row, last_col, value);
                }
            }
        }
    }

    return answer;
}
