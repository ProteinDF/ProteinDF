#include "TlCommunicate.h"
#include "TlColVectorMatrix.h"


TlColVectorMatrix::TlColVectorMatrix(const index_type row,
                                     const index_type col)
{
    this->globalRows_ = 0;
    this->globalCols_ = 0;
    this->reserveRows_ = 10;
    this->resize(row, col);
}


TlColVectorMatrix::~TlColVectorMatrix()
{
}


void TlColVectorMatrix::resize(const index_type new_globalRows,
                               const index_type new_globalCols)
{
    const index_type prev_globalCols = this->globalCols_;

    if (prev_globalCols != new_globalCols) {
        // remake PE table
        TlCommunicate& rComm = TlCommunicate::getInstance();
        const int numOfProcs = rComm.getNumOfProcs();
        const int myRank = rComm.getRank();

        this->col_PE_table_.resize(new_globalCols);
        int pe = 0;
        for (index_type c = 0; c < new_globalCols; ++c) {
            this->col_PE_table_[c] = pe;
            
            ++pe;
            if (pe >= numOfProcs) {
                pe = 0;
            }
        }

        std::vector<index_type> localColTable;
        localColTable.reserve(new_globalCols / numOfProcs +1);
        for (index_type c = 0; c < new_globalCols; ++c) {
            if (this->col_PE_table_[c] == myRank) {
                localColTable.push_back(c);
            }
        }

        const std::size_t prev_localColTableSize = this->data_.size();
        const std::size_t new_localColTableSize = localColTable.size();
        if (prev_localColTableSize < new_localColTableSize) {
            this->data_.resize(new_localColTableSize);
            for (std::size_t i = prev_localColTableSize; i < new_localColTableSize; ++i) {
                this->data_[i].col = localColTable[i];
            }
        } else if (prev_localColTableSize < new_localColTableSize) {
            this->data_.resize(new_localColTableSize);
        }
    }

    std::vector<ColVector>::iterator itEnd = this->data_.end();
    for (std::vector<ColVector>::iterator it = this->data_.begin(); it != itEnd; ++it) {
        it->rows.reserve(this->reserveRows_);
        it->rows.resize(new_globalRows);
    }

    this->globalRows_ = new_globalRows;
    this->globalCols_ = new_globalCols;
}


void TlColVectorMatrix::reserve_rows(const index_type rows) {
    this->reserveRows_ = std::max(globalRows_, rows);
}


void TlColVectorMatrix::set(index_type row, index_type col, double value)
{
    std::vector<ColVector>::iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), ColVector(col));

    if ((it != this->data_.end()) && (it->col == col)) {
        it->rows[row] = value;
    }
}


TlVector TlColVectorMatrix::getColVector(index_type col) const
{
    TlVector answer;
    std::vector<ColVector>::const_iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), ColVector(col));
    if (it != this->data_.end()) {
        assert(it->col == col);
        answer = it->rows;
    }

    return answer;
}


TlMatrixObject::index_type 
TlColVectorMatrix::getColVector(index_type col,
                                double *pBuf,
                                index_type maxRowSize) const
{
    index_type copySize = 0;
    std::vector<ColVector>::const_iterator it = 
        std::lower_bound(this->data_.begin(), this->data_.end(), ColVector(col));
    if (it != this->data_.end()) {
        assert(it->col == col);
        assert(it->rows.size() == this->getNumOfRows());
        copySize = std::min(this->getNumOfRows(), maxRowSize);
        std::copy(it->rows.begin(), it->rows.begin() + copySize,
                  pBuf);
    }
    return copySize;
}


int TlColVectorMatrix::getPEinChargeByCol(const index_type col) const
{
    assert((0 <= col) && (col < this->getNumOfCols()));
    return this->col_PE_table_[col];
}


TlMatrix TlColVectorMatrix::getTlMatrix() const 
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlMatrix answer(row, col);
    
    std::vector<ColVector>::const_iterator itEnd = this->data_.end();
    for (std::vector<ColVector>::const_iterator it = this->data_.begin(); it != itEnd; ++it) {
        const index_type c = it->col;
        for (index_type r = 0; r < row; ++r) {
            answer.set(r, c, it->rows[r]);
        }
    }

    rComm.allReduce_SUM(answer);
    return answer;
}


TlDistributeMatrix TlColVectorMatrix::getTlDistributeMatrix() const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type row = this->getNumOfRows();
    const index_type col = this->getNumOfCols();
    TlDistributeMatrix answer(row, col);

    int blocks = 0;
    std::vector<index_type> cols;
    std::vector<double> vtr;
    for (int pe = 0; pe < numOfProcs; ++pe) {
        if (pe == myRank) {
            blocks = this->data_.size();
            vtr.resize(row * blocks);
            for (int i = 0; i < blocks; ++i) {
                cols[i] = this->data_[i].col;
                for (index_type r = 0; r < row; ++r) {
                    vtr[row * i + r] = this->data_[i].rows[r];
                }
            }
        }

        rComm.broadcast(blocks, pe);
        rComm.broadcast(cols, pe);
        rComm.broadcast(vtr, pe);

        // set
        for (int i = 0; i < blocks; ++i) {
            const index_type c = cols[i];
            for (index_type r = 0; r < row; ++r) {
                answer.set(r, c, vtr[row * i + r]);
            }
        }
    }

    return answer;
}


void TlColVectorMatrix::save(const std::string& basename,
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
    int numOfMyCols = this->col_PE_table_.size();
    ofs.write(reinterpret_cast<const char*>(&PE), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&row), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&col), sizeof(index_type));

    // data
    ofs.write(reinterpret_cast<const char*>(&numOfMyCols), sizeof(int));
    for (int i = 0; i < numOfMyCols; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i].col)), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i].rows[0])), sizeof(double) * row);
    }
    
    ofs.close();
}


void TlColVectorMatrix::load(const std::string& basename,
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
    int numOfMyCols = 0;
    ifs.read((char*)&numOfMyCols, sizeof(numOfMyCols));
    assert(numOfMyCols == this->col_PE_table_.size());
    for (int i = 0; i < numOfMyCols; ++i) {
        index_type col_index = 0;
        ifs.read((char*)&col_index, sizeof(index_type));
        assert(this->data_[i].col == col_index);
        assert(this->data_[i].rows.size() == row);
        ifs.read((char*)&(this->data_[i].rows[0]), sizeof(double)*row);
    }

    ifs.close();
}

