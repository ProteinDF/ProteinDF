#include "TlCatMatrix.h"

TlCatMatrix::TlCatMatrix() : catMatrix(NULL)
{
}

TlCatMatrix::~TlCatMatrix()
{
    if (this->catMatrix != NULL) {
        delete this->catMatrix;
        this->catMatrix = NULL;
    }
}

void TlCatMatrix::addMatrix(const TlMatrix& obj)
{
    const int baseMatrixRows = this->catMatrix->getNumOfRows();
    const int baseMatrixCols = this->catMatrix->getNumOfRows();

    const int

}
