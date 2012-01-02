#ifndef TLFILEMATRIX_H
#define TLFILEMATRIX_H

#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cassert>

#include "TlMatrixObject.h"
#include "TlMatrix.h"
#include "TlVector.h"

// default cache size (16 MB)
#define DEFAULT_CACHE_SIZE (16 * 1024 * 1024)

class TlFileMatrix : public TlMatrixObject {
    friend class TlCommunicate;

protected:
    // for cache
    struct CacheUnit {
public:
        CacheUnit(unsigned long g) : group(g), isUpdate(false) {
        }

public:
        unsigned long group;
        bool isUpdate;
        std::vector<double> data;

        friend bool operator==(const CacheUnit& x, const CacheUnit& y) {
            return (x.group == y.group);
        }
    };

    struct CacheUnitComp : public std::unary_function<const CacheUnit&, bool> {
public:
        CacheUnitComp(unsigned long group) : group_(group) {
        }

        bool operator()(const CacheUnit& cu) {
            return (group_ == cu.group);
        }

private:
        unsigned long group_;
    };

public:
    explicit TlFileMatrix(const std::string& filePath, int row = 0, int col = 0, size_t cacheSize = DEFAULT_CACHE_SIZE);
    virtual ~TlFileMatrix();

protected:
    TlFileMatrix(const std::string& filePath, int row, int col, bool initialize, size_t cacheSize = DEFAULT_CACHE_SIZE);

public:
    int getNumOfRows() const {
        return this->numOfRows_;
    }

    int getNumOfCols() const {
        return this->numOfCols_;
    }

    virtual std::size_t getMemSize() const;
    
    virtual void set(int row, int col, double value);
    virtual double get(int row, int col) const;
    virtual void add(int row, int col, double value);

    virtual TlFileMatrix& operator*=(double coef);
    virtual TlFileMatrix& operator/=(double coef) {
        return this->operator*=(1.0 / coef);
    }

    /// 指定した行の要素から構成されるベクトルを返す
    ///
    /// @param[in] nRow 指定する行
    virtual TlVector getRowVector(int nRow) const;

    /// 指定した列の要素から構成されるベクトルを返す
    ///
    /// @param[in] nCol 指定する列
    virtual TlVector getColumnVector(int nCol) const;

    /// ブロック行列を返す
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] row_distance 取得する行数
    /// @param[in] col_distance 取得する列数
    /// @return row_distance × col_distance 次元のブロック行列
    virtual TlMatrix getBlockMatrix(index_type row, index_type col,
                                    index_type rowDistance, index_type colDistance) const;

    /// 行列要素を指定された位置に上書きする
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] matrix 行列要素
    virtual void setBlockMatrix(index_type row, index_type col,
                                const TlMatrix& matrix);
protected:
    virtual bool readHeader();

    virtual std::size_t index(const int row, const int col) const {
        assert(0 <= row);
        assert(row < this->numOfRows_);
        assert(0 <= col);
        assert(col < this->numOfCols_);

        return (std::size_t(row) * std::size_t(this->numOfCols_) + std::size_t(col));
    }

    virtual std::size_t maxIndex() const {
        return (std::size_t(this->getNumOfRows()) * std::size_t(this->getNumOfCols()));
    }

    double* getCachedData(const int row, const int col);
    double getCachedData(const int row, const int col) const;
    void updateCache(size_t index) const;

protected:
    virtual void open();
    void writeDisk(const CacheUnit& cu) const;

protected:
    virtual bool load(const std::string& path) {
        return false;
    }

    virtual bool save(const std::string& path) const {
        return false;
    }
    
protected:
    std::string filePath_;
    int numOfRows_;
    int numOfCols_;

    mutable std::fstream fs_;
    std::fstream::pos_type startPos_;
    std::fstream::pos_type endPos_;

    mutable std::list<CacheUnit> cache_;
    mutable size_t cacheCount_; // == cache_.size()
    size_t cacheSize_;
};

#endif // TLFILEMATRIX_H
