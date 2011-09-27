#include "TlSparseHashSymmetricMatrix.h"

TlSparseHashSymmetricMatrix::TlSparseHashSymmetricMatrix(const int size) : TlSparseHashMatrix(size, size)
{
}

TlSparseHashSymmetricMatrix::TlSparseHashSymmetricMatrix(const TlSparseHashSymmetricMatrix& rhs)
        : TlSparseHashMatrix(rhs)
{
}

TlSparseHashSymmetricMatrix::~TlSparseHashSymmetricMatrix()
{
}

void TlSparseHashSymmetricMatrix::resize(const int size)
{
    TlSparseHashMatrix::resize(size, size);
}

