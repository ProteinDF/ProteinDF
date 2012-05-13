#ifndef TLCOLVECTORMATRIX2_H
#define TLCOLVECTORMATRIX2_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"

class TlColVectorMatrix2 {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    explicit TlColVectorMatrix2(index_type row = 1, index_type col = 1,
                                int allProcs = 1, int rank = 0);
    TlColVectorMatrix2(const TlColVectorMatrix2& rhs);
    ~TlColVectorMatrix2();
        
public:
    TlColVectorMatrix2& operator=(const TlColVectorMatrix2& rhs);

    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->numOfRows_;
    };

    index_type getNumOfCols() const {
        return this->numOfCols_;
    };

    /// 行数のcapacityを設定する
    void reserve_rows(index_type row);

    void set(index_type row, index_type col, double value);
        
    TlVector getColVector(index_type col) const;
    index_type getColVector(index_type col, double *pBuf, index_type maxRowSize) const;
        
    int getPEinChargeByCol(index_type col) const;
        
public:
    void save(const std::string& basename) const;
    void load(const std::string& basename);

private:
    index_type numOfRows_;
    index_type numOfCols_;
    index_type reserveRows_; // 行数のメモリ確保量

    int allProcs_;
    int rank_;

    std::vector<std::vector<double> > data_;
};

#endif // TLCOLVECTORMATRIX2_H
