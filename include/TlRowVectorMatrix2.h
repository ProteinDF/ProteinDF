#ifndef TLROWVECTORMATRIX2_H
#define TLROWVECTORMATRIX2_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"

class TlRowVectorMatrix2 {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    explicit TlRowVectorMatrix2(index_type row = 1, index_type col = 1,
                                int allProcs = 1, int rank = 0);

    ~TlRowVectorMatrix2();
        
public:

    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->numOfRows_;
    };

    index_type getNumOfCols() const {
        return this->numOfCols_;
    };

    /// 列数のcapacityを設定する
    void reserve_cols(index_type col);

    void set(index_type row, index_type col, double value);
        
    TlVector getRowVector(index_type row) const;
    index_type getRowVector(index_type row, double *pBuf, index_type maxColSize) const;
        
    int getPEinChargeByRow(index_type row) const;
        
    TlMatrix getTlMatrix() const;

private:
    TlRowVectorMatrix2(const TlRowVectorMatrix2& rhs);
    TlRowVectorMatrix2& operator=(const TlRowVectorMatrix2& rhs);

private:
    index_type numOfRows_;
    index_type numOfCols_;
    index_type reserveCols_; // 列数のメモリ確保量

    int allProcs_;
    int rank_;

    index_type numOfLocalRows_;
    std::vector<double* > data_;
};

#endif // TLROWVECTORMATRIX2_H
