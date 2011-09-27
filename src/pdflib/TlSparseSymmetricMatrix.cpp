#include "TlSparseSymmetricMatrix.h"

TlSparseSymmetricMatrix::TlSparseSymmetricMatrix(const int size) : TlSparseMatrix(size, size)
{
}

TlSparseSymmetricMatrix::TlSparseSymmetricMatrix(const TlSparseSymmetricMatrix& rhs) : TlSparseMatrix(rhs)
{
}

TlSparseSymmetricMatrix::~TlSparseSymmetricMatrix()
{
}

void TlSparseSymmetricMatrix::resize(const int size)
{
    TlSparseMatrix::resize(size, size);
}


TlVector TlSparseSymmetricMatrix::getRowVector(const index_type row) const
{
    const int numOfCols = this->getNumOfCols();
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
            } else if (c == row) {
                ans[r] = p->second;
            }
        }
    }

    return ans;
}


TlVector TlSparseSymmetricMatrix::getColVector(const index_type col) const
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
            } else if (r == col) {
                ans[c] = p->second;
            }
        }
    }

    return ans;
}

const TlSparseSymmetricMatrix& TlSparseSymmetricMatrix::dot(const TlSparseSymmetricMatrix& X)
{
    iterator p = this->begin();
    iterator pEnd = this->end();
    const_iterator q = X.begin();
    const_iterator qEnd = X.end();

    while ((p != pEnd) && (q != qEnd)) {
        const unsigned long p_index = p->first;
        const unsigned long q_index = q->first;
        if (p_index < q_index) {
            ++p;
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


