#ifndef TLCOLVECTORMATRIX_H
#define TLCOLVECTORMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"
#include "TlDistributeMatrix.h"

class TlColVectorMatrix {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    TlColVectorMatrix(index_type row = 1, index_type col = 1);
    ~TlColVectorMatrix();
        
public:
    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->globalRows_;
    };

    index_type getNumOfCols() const {
        return this->globalCols_;
    };

    /// 行数のcapacityを設定する
    void reserve_rows(index_type col);

    void set(index_type row, index_type col, double value);
        
    TlVector getColVector(index_type col) const;
    index_type getColVector(index_type row, double *pBuf, index_type maxRowSize) const;
        
    int getPEinChargeByCol(index_type col) const;
        
    TlMatrix getTlMatrix() const;
    TlDistributeMatrix getTlDistributeMatrix() const;

public:
    void save(const std::string& basename,
              const int PE = 0) const;
    void load(const std::string& basename,
              const int PE = 0);

private:
    struct ColVector {
    public:
        explicit ColVector(index_type r =1, index_type c =0)
            : rows(r), col(c) {
        };
            
        bool operator<(const ColVector& rhs) const {
            return (this->col < rhs.col);
        };
            
    public:
        index_type col;
        std::vector<double> rows;
    };
        
private:
    index_type globalRows_;
    index_type globalCols_;
    index_type reserveRows_; // 行数のメモリ確保量
    std::vector<ColVector> data_;

    std::vector<int> col_PE_table_;
};

#endif // TLCOLVECTORMATRIX_H
