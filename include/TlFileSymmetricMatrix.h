#ifndef TLFILESYMMETRICMATRIX_H
#define TLFILESYMMETRICMATRIX_H

#include "TlFileMatrix.h"
#include "TlPartialSymmetricMatrix.h"

class TlFileSymmetricMatrix : public TlFileMatrix {
public:
    explicit TlFileSymmetricMatrix(const std::string& filePath, int dim = 0, size_t cacheSize = DEFAULT_CACHE_SIZE);
    virtual ~TlFileSymmetricMatrix();

public:
    void add(int row, int col, double value);
    void add(const TlPartialSymmetricMatrix& psm);

    virtual TlFileSymmetricMatrix& operator*=(double coef);
    virtual TlFileSymmetricMatrix& operator/=(double coef) {
        return this->operator*=(1.0 / coef);
    }

    TlPartialSymmetricMatrix getPartialMatrix(const int startRow, const int startCol, const int range) const;

protected:
    virtual void open();

protected:
    virtual bool readHeader();

    virtual size_t index(int row, int col) const {
        if (row < col) {
            std::swap(row, col);
        }

        return (std::size_t(col) + (std::size_t(row +1) * std::size_t(row)) / std::size_t(2));
    }

    virtual size_t maxIndex() const {
        return (std::size_t(this->getNumOfRows()) * std::size_t(this->getNumOfCols() +1) / std::size_t(2));
    }
};

#endif // TLFILESYMMETRICMATRIX_H
