#ifndef TLPARTIALMATRIX_H
#define TLPARTIALMATRIX_H

#include <vector>
#include "TlMatrixObject.h"
#include "TlUtils.h"

class TlPartialMatrix : public TlMatrixObject {
public:
    explicit TlPartialMatrix(int globalNumOfRows =0, int globalNumOfCols =0,
                             int startRow =0, int startCol =0,
                             int rowRange =0, int colRange =0);
    TlPartialMatrix(const TlPartialMatrix& rhs);
    virtual ~TlPartialMatrix();

public:
    TlPartialMatrix& operator=(const TlPartialMatrix& rhs);

public:
    index_type getNumOfRows() const;
    index_type getNumOfCols() const;

    index_type getStartRow() const;
    index_type getStartCol() const;
    
    index_type getRowRange() const;
    index_type getColRange() const;

    virtual std::size_t getMemSize() const;
    
    double getMaxAbsoluteElement(index_type* outRow =NULL, index_type* outCol =NULL) const;

public:
    virtual void set(index_type globalRow, index_type globalCol, double value);
    virtual void add(index_type globalRow, index_type globalCol, double value);
    virtual double get(index_type globalRow, index_type globalCol) const;
    virtual double getLocal(index_type localRow, index_type localCol) const;

    virtual TlVector getRowVector(index_type row) const;
    virtual TlVector getColVector(index_type col) const;
    
protected:
    virtual inline size_type index(index_type globalRow, index_type globalCol) const;

protected:
    virtual bool load(const std::string& path) {
        return false;
    }

    virtual bool save(const std::string& path) const {
        return false;
    }
    
protected:
    index_type numOfRows_;
    index_type numOfCols_;
    index_type startRow_;
    index_type startCol_;
    index_type localNumOfRows_;
    index_type localNumOfCols_;

    double* pData_;

    friend class TlCommunicate;
};


inline TlMatrixObject::index_type TlPartialMatrix::getRowRange() const
{
    return this->localNumOfRows_;
}


inline TlMatrixObject::index_type TlPartialMatrix::getColRange() const
{
    return this->localNumOfCols_;
}


inline TlMatrixObject::size_type TlPartialMatrix::index(const index_type globalRow,
                                                        const index_type globalCol) const
{
    size_type answer = -1;

    const index_type r = globalRow - this->startRow_;
    const index_type c = globalCol - this->startCol_;
    if ((0 <= r) && (r < this->getRowRange()) &&
        (0 <= c) && (c < this->getColRange())) {
        answer = r * this->localNumOfCols_ + c;
    }

    return answer;
}


#endif // TLPARTIALMATRIX_H

