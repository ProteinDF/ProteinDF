#include <cmath>
#include <limits>
#include "TlMatrixObject.h"

void TlMatrixObject::addByList(const index_type* pIndexPairs,
                               const double* pValues,
                               const std::size_t size)
{
    for (std::size_t i = 0; i < size; ++i) {
        const index_type globalRow = pIndexPairs[i*2   ];
        const index_type globalCol = pIndexPairs[i*2 +1];
        const double value = pValues[i];
        
        this->add(globalRow, globalCol, value);
    }
}

TlVector TlMatrixObject::getRowVector(const int row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));

    const int numOfCols = this->getNumOfCols();
    TlVector answer(numOfCols);

    for (int i = 0; i < numOfCols; ++i) {
        answer[i] = this->get(row, i);
    }

    return answer;
}


TlVector TlMatrixObject::getColVector(const int col) const
{
    assert((0 <= col) && (col < this->getNumOfCols()));

    const int numOfRows = this->getNumOfRows();
    TlVector answer(numOfRows);

    for (int i = 0; i < numOfRows; ++i) {
        answer[i] = this->get(i, col);
    }

    return answer;
}


double TlMatrixObject::getMaxAbsoluteElement(int* pOutRow,
                                             int* pOutCol) const
{
    const int maxRow = this->getNumOfRows();
    const int maxCol = this->getNumOfCols();
    int outRow = -1;
    int outCol = -1;
    double value = 0.0;
    for (int r = 0; r < maxRow; ++r) {
        for (int c = 0; c < maxCol; ++c) {
            const double tmp = std::fabs(this->get(r, c));
            if (value < tmp) {
                outRow = r;
                outCol = c;
                value = tmp;
            }
        }
    }

    if (pOutRow != NULL) {
        *pOutRow = outRow;
    }
    if (pOutCol != NULL) {
        *pOutCol = outCol;
    }

    return value;
}


// TlMatrixObject& TlMatrixObject::operator+=(const TlMatrixObject& rhs)
// {
//     const index_type numOfRows = this->getNumOfRows();
//     const index_type numOfCols = this->getNumOfCols();
//     assert(numOfRows == rhs.getNumOfRows());
//     assert(numOfCols == rhs.getNumOfCols());

//     for (index_type r = 0; r < numOfRows; ++r) {
//         for (index_type c = 0; c < numOfCols; ++c) {
//             this->add(r, c, rhs.get(r, c));
//         }
//     }

//     return *this;
// }


