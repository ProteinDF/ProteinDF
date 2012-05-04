#include <cassert>
#include <cmath>
#include <limits>

#include "TlSparseMatrix.h"

const int TlSparseMatrix::INT_BITS = sizeof(int) * 8;
const int TlSparseMatrix::MAX_INT = std::numeric_limits<unsigned int>::max();

TlSparseMatrix::TlSparseMatrix(int row, int col)
        : m_nRows(row), m_nCols(col)
{
}

TlSparseMatrix::TlSparseMatrix(const TlSparseMatrix& rhs)
        : m_nRows(rhs.m_nRows), m_nCols(rhs.m_nCols), m_aMatrix(rhs.m_aMatrix)
{
}

TlSparseMatrix::~TlSparseMatrix()
{
}


TlMatrixObject::index_type TlSparseMatrix::getNumOfRows() const
{
    return this->m_nRows;
}


TlMatrixObject::index_type TlSparseMatrix::getNumOfCols() const
{
    return this->m_nCols;
}


TlMatrixObject::size_type TlSparseMatrix::getSize() const
{
    return this->m_aMatrix.size();
}


std::size_t TlSparseMatrix::getMemSize() const
{
    std::size_t answer = sizeof(TlSparseMatrix);
    answer += (sizeof(int) * 2 + sizeof(double)) * this->getSize();

    return answer;
}


void TlSparseMatrix::clear()
{
    this->m_aMatrix.clear();
}


void TlSparseMatrix::erase(int row, int col)
{
    //const long index = (row << (sizeof(int) * 8)) + col;

    //iterator p = this->m_aMatrix.find(TlMatrixIndexPair(row, col));
    iterator p = this->m_aMatrix.find(this->index(row, col));
    if (p != this->m_aMatrix.end()) {
        this->erase(p);
    }
}


void TlSparseMatrix::erase(TlSparseMatrix::iterator p)
{
    this->m_aMatrix.erase(p);
}


void TlSparseMatrix::resize(const int row, const int col)
{
    assert(0 < row);
    assert(0 < col);

    if ((row > this->getNumOfRows()) || (col > this->getNumOfCols())) {
        for (iterator p = this->m_aMatrix.begin(); p != this->m_aMatrix.end(); ++p) {
            int r;
            int c;
            this->index(p->first, &r, &c);

            if ((r >= row) || (c > col)) {
                this->m_aMatrix.erase(p);
            }
        }
    }

    this->m_nRows = row;
    this->m_nCols = col;
}


double TlSparseMatrix::pop(int* pRow, int* pCol)
{
    assert(pRow != NULL);
    assert(pCol != NULL);

    double answer = 0.0;

    if (this->m_aMatrix.size() > 0) {
        //iterator p = this->m_aMatrix.end();
        //--p;
        iterator p = this->m_aMatrix.begin();
        this->index(p->first, pRow, pCol);
        answer = p->second;
        this->erase(p);
    }

    return answer;
}

void TlSparseMatrix::merge(const TlSparseMatrix& rhs)
{
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    TlSparseMatrix::const_iterator pEnd = rhs.end();
    for (TlSparseMatrix::const_iterator p = rhs.begin(); p != pEnd; ++p) {
        this->m_aMatrix[p->first] = p->second;
    }
}

TlSparseMatrix& TlSparseMatrix::operator=(const TlSparseMatrix& rhs)
{
    this->m_nRows = rhs.m_nRows;
    this->m_nCols = rhs.m_nCols;
    this->m_aMatrix.clear();
    this->m_aMatrix = rhs.m_aMatrix;

    return (*this);
}

TlSparseMatrix& TlSparseMatrix::operator*=(const double& rhs)
{
    for (iterator p = this->begin(); p != this->end(); p++) {
        p->second *= rhs;
    }

    return (*this);
}

TlSparseMatrix& TlSparseMatrix::operator/=(const double& rhs)
{
    assert(std::fabs(rhs) > 1E-16);

    const double v = 1.0 / rhs;

    return this->operator*=(v);
}

TlVector TlSparseMatrix::getRowVector(const index_type row) const
{
    const index_type numOfCols = this->getNumOfCols();
    TlVector ans(numOfCols);

    if (this->m_aMatrix.size() > (std::size_t)(numOfCols / 2)) {
        for (index_type i = 0; i < numOfCols; ++i) {
            ans.set(i, this->get(row, i));
        }
    } else {
        const_iterator pEnd = this->end();
        for (const_iterator p = this->begin(); p != pEnd; ++p) {
            int r, c;
            this->index(p->first, &r, &c);
            if (r == row) {
                ans[c] = p->second;
            }
        }
    }

    return ans;
}

TlVector TlSparseMatrix::getColVector(const int col) const
{
    const int numOfRows = this->getNumOfRows();
    TlVector ans(numOfRows);

    if (this->m_aMatrix.size() > (std::size_t)(numOfRows / 2)) {
        for (index_type i = 0; i < numOfRows; ++i) {
            ans.set(i, this->get(i, col));
        }
    } else {
        const_iterator pEnd = this->end();
        for (const_iterator p = this->begin(); p != pEnd; ++p) {
            int r, c;
            this->index(p->first, &r, &c);
            if (c == col) {
                ans[r] = p->second;
            }
        }
    }

    return ans;
}

const TlSparseMatrix& TlSparseMatrix::dot(const TlSparseMatrix& X)
{
    iterator p = this->begin();
    iterator pEnd = this->end();
    const_iterator q = X.begin();
    const_iterator qEnd = X.end();

    while ((p != pEnd) && (q != qEnd)) {
        const unsigned long p_index = p->first;
        const unsigned long q_index = q->first;
        if (p_index < q_index) {
            this->m_aMatrix.erase(p++);
        } else if (p_index > q_index) {
            ++q;
        } else {
            // p_index == q_index
            this->m_aMatrix[p_index] *= X.m_aMatrix[p_index];
            ++p;
            ++q;
        }
    }

    return (*this);
}

double TlSparseMatrix::sum() const
{
    double ans = 0.0;

    const_iterator pEnd = this->end();
    for (const_iterator p = this->begin(); p != pEnd; ++p) {
        ans += p->second;
    }

    return ans;
}

std::vector<int> TlSparseMatrix::getRowIndexList() const
{
    std::vector<int> rowIndexList;
    rowIndexList.reserve(this->m_aMatrix.size());

    const_iterator pEnd = this->m_aMatrix.end();
    for (const_iterator p = this->m_aMatrix.begin(); p != pEnd; ++p) {
        int row = 0;
        int col = 0;
        this->index(p->first, &row, &col);
        rowIndexList.push_back(row);
    }

    std::sort(rowIndexList.begin(), rowIndexList.end());
    std::unique(rowIndexList.begin(), rowIndexList.end());
    std::vector<int>(rowIndexList).swap(rowIndexList); // see effective STL 17.

    return rowIndexList;
}

std::vector<int> TlSparseMatrix::getColIndexList() const
{
    std::vector<int> colIndexList;
    colIndexList.reserve(this->m_aMatrix.size());

    const_iterator pEnd = this->m_aMatrix.end();
    for (const_iterator p = this->m_aMatrix.begin(); p != pEnd; ++p) {
        int row = 0;
        int col = 0;
        this->index(p->first, &row, &col);
        colIndexList.push_back(col);
    }

    std::sort(colIndexList.begin(), colIndexList.end());
    std::unique(colIndexList.begin(), colIndexList.end());
    std::vector<int>(colIndexList).swap(colIndexList); // see effective STL 17.

    return colIndexList;
}


bool TlSparseMatrix::save(const std::string& filePath) const
{
    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

    bool answer = this->save(ofs);

    ofs.close();
    return answer;
}


bool TlSparseMatrix::save(std::ofstream& ofs) const
{
    const int type = 128; // means SparseMatrix
    ofs.write(reinterpret_cast<const char*>(&type), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));
    unsigned long counts = this->m_aMatrix.size();
    ofs.write(reinterpret_cast<const char*>(&counts), sizeof(unsigned long));

    index_type row = 0;
    index_type col = 0;
    SparseMatrixData::const_iterator itEnd = this->m_aMatrix.end();
    for (SparseMatrixData::const_iterator it = this->m_aMatrix.begin(); it != itEnd; ++it) {
        this->index(it->first, &row, &col);
        ofs.write(reinterpret_cast<const char*>(&row), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&col), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&it->second), sizeof(double));
    }

    return true;
}


bool TlSparseMatrix::load(const std::string& filePath)
{
    std::ifstream ifs;

    ifs.open(filePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlMatrix::load(): could not open file. " << filePath << std::endl;
#endif // DEBUG
        return false;
    }

    bool answer = this->load(ifs);
    ifs.close();

    if (answer != true) {
        std::cerr << "TlMatrix::load() is not supported: " << filePath << std::endl;
        return false;
    }

    return true;
}


bool TlSparseMatrix::load(std::ifstream& ifs)
{
    // read header
    int type = 0;
    int numOfRows = 0;
    int numOfCols = 0;
    unsigned long counts = 0;
    ifs.read((char*)&type, sizeof(int));
    ifs.read((char*)&numOfRows, sizeof(int));
    ifs.read((char*)&numOfCols, sizeof(int));
    ifs.read((char*)&counts, sizeof(unsigned long));
    assert(type == 128);
    this->resize(numOfRows, numOfCols);

    int row = 0;
    int col = 0;
    double value = 0.0;
    for (unsigned long i = 0; i < counts; ++i) {
        ifs.read((char*)&row, sizeof(int));
        ifs.read((char*)&col, sizeof(int));
        ifs.read((char*)&value, sizeof(double));
        this->set(row, col, value);
    }

    return true;
}
