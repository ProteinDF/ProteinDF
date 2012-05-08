#ifndef TLROWVECTORMATRIX_H
#define TLROWVECTORMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"
#include "TlDistributeMatrix.h"
#include "TlColVectorMatrix.h"

class TlRowVectorMatrix {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    TlRowVectorMatrix(index_type row = 1, index_type col = 1);
    ~TlRowVectorMatrix();
        
public:
    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->globalRows_;
    };

    /// 列数のcapacityを設定する
    void reserve_cols(index_type col);

    index_type getNumOfCols() const {
        return this->globalCols_;
    };

    void set(index_type row, index_type col, double value);
        
    TlVector getRowVector(index_type row) const;
    index_type getRowVector(index_type row, double *pBuf, index_type maxColSize) const;
        
    int getPEinChargeByRow(index_type row) const;
        
    TlMatrix getTlMatrix() const;
    TlDistributeMatrix getTlDistributeMatrix() const;
    TlColVectorMatrix getTlColVectorMatrix() const;

public:
    void save(const std::string& basename,
              const int PE = 0) const;
    void load(const std::string& basename,
              const int PE = 0);

private:
    struct RowVector {
    public:
        explicit RowVector(index_type r =0, index_type c =1) 
            : row(r), cols(c) {
        };
            
        bool operator<(const RowVector& rhs) const {
            return (this->row < rhs.row);
        };
            
    public:
        index_type row;
        std::vector<double> cols;
    };
        
private:
    index_type globalRows_;
    index_type globalCols_;
    index_type reserveCols_; // 列数のメモリ確保量
    std::vector<RowVector> data_;

    std::vector<int> row_PE_table_;
};

#endif // TLROWVECTORMATRIX_H
