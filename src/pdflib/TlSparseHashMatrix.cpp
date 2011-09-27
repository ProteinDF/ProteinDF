#include <cassert>
#include <limits>

#include "TlSparseHashMatrix.h"

const int TlSparseHashMatrix::INT_BITS = sizeof(int) * 8;
const int TlSparseHashMatrix::MAX_INT = std::numeric_limits<unsigned int>::max();

TlSparseHashMatrix::TlSparseHashMatrix(int rows, int cols)
        : numOfRows(rows), numOfCols(cols)
{
}

TlSparseHashMatrix::TlSparseHashMatrix(const TlSparseHashMatrix& rhs)
        : numOfRows(rhs.numOfRows), numOfCols(rhs.numOfCols), container(rhs.container)
{
}

TlSparseHashMatrix::~TlSparseHashMatrix()
{
}

int TlSparseHashMatrix::getNumOfRows() const
{
    return this->numOfRows;
}

int TlSparseHashMatrix::getNumOfCols() const
{
    return this->numOfCols;
}

void TlSparseHashMatrix::resize(const int row, const int col)
{
    assert(0 < row);
    assert(0 < col);

    if ((row > this->getNumOfRows()) || (col > this->getNumOfCols())) {
        ContainerType::iterator pEnd = this->container.end();
        for (ContainerType::iterator p = this->container.begin(); p != pEnd; ++p) {
            int r;
            int c;
            this->index(p->first, &r, &c);

            if ((r >= row) || (c > col)) {
                this->container.erase(p);
            }
        }
    }

    this->numOfRows = row;
    this->numOfCols = col;
}

int TlSparseHashMatrix::getSize() const
{
    return this->container.size();
}

TlSparseHashMatrix& TlSparseHashMatrix::operator=(const TlSparseHashMatrix& rhs)
{
    this->numOfRows = rhs.numOfRows;
    this->numOfCols = rhs.numOfCols;
    this->container.clear();
    this->container = rhs.container;

    return (*this);
}

double& TlSparseHashMatrix::operator()(const int row, const int col)
{
    const unsigned long i = this->index(row, col);

    return (this->container[i]);
}

