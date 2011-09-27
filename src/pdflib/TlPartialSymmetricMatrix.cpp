#include <iostream>
#include <cassert>
#include "TlPartialSymmetricMatrix.h"
#include "TlUtils.h"

//#define WARNING_OUT_OF_RANGE


TlPartialSymmetricMatrix::TlPartialSymmetricMatrix(int numOfSize, int startRow, int startCol, int range)
    : TlPartialMatrix(numOfSize, numOfSize, startRow, startCol, range, range)
{
}


TlPartialSymmetricMatrix::TlPartialSymmetricMatrix(index_type numOfSize, index_type startRow, index_type startCol,
                                                   index_type rowRange, index_type colRange)
    : TlPartialMatrix(numOfSize, numOfSize, startRow, startCol, rowRange, colRange)
{
}


TlPartialSymmetricMatrix::TlPartialSymmetricMatrix(const TlPartialSymmetricMatrix& rhs) : TlPartialMatrix(rhs)
{
}


TlPartialSymmetricMatrix::~TlPartialSymmetricMatrix()
{
}


TlPartialSymmetricMatrix& TlPartialSymmetricMatrix::operator=(const TlPartialSymmetricMatrix& rhs)
{
    if (&rhs != this) {
        TlPartialMatrix::operator=(rhs);
    }

    return *this;
}


void TlPartialSymmetricMatrix::set(const index_type globalRow, const index_type globalCol, const double value)
{
    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        this->pData_[index] = value;
    }
}


void TlPartialSymmetricMatrix::add(const index_type globalRow, const index_type globalCol, const double value)
{
    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        this->pData_[index] += value;
    }
}


double TlPartialSymmetricMatrix::get(const index_type globalRow, const index_type globalCol) const
{
    double answer = 0.0;

    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        answer = this->pData_[index];
    }
#ifdef WARNING_OUT_OF_RANGE
    else {
        std::cerr << "out of range: TlPartialSymmetricMatrix::get() row="
                  << globalRow << ", col=" << globalCol
                  << std::endl;
    }
#endif // WARNING_OUT_OF_RANGE

    return answer;
}


double TlPartialSymmetricMatrix::getLocal(const index_type localRow, const index_type localCol) const
{
    assert(0 <= localRow);
    assert(localRow < this->localNumOfRows_);
    assert(0 <= localCol);
    assert(localCol < this->localNumOfCols_);

    return this->pData_[localRow * this->localNumOfCols_ + localCol];
}


TlVector TlPartialSymmetricMatrix::getRowVector(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));

    const index_type numOfCols = this->getNumOfCols();
    TlVector answer(numOfCols);

    const index_type startRow = this->getStartRow();
    const index_type endRow = startRow + this->getRowRange();
    const index_type startCol = this->getStartCol();
    const index_type endCol = startCol + this->getColRange();

    // 探査領域を狭める
    if ((startRow <= row) && (row < endRow)) {
        for (index_type i = startCol; i < endCol; ++i) {
            answer.set(i, this->get(row, i));
        }
    } else if ((startCol <= row) && (row < endCol)) {
        for (index_type i = startRow; i < endRow; ++i) {
            answer.set(i, this->get(i, row));
        }
    }
    
    return answer;
}


TlVector TlPartialSymmetricMatrix::getColVector(const index_type col) const
{
    return this->getRowVector(col);
}
