#ifndef TLSPARSEMATRIX_H
#define TLSPARSEMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <cassert>
#include <functional>
#include <iostream>

// -------------------------------------
#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#elifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elifdef HAVE_GOOGLE_SPARSE_HASH_MAP
#include <google/sparse_hash_map>
#else
#include <map>
#endif
// -------------------------------------

#include "TlMatrixObject.h"
#include "TlUtils.h"
#include "TlVector.h"

/// 疎行列クラス
class TlSparseMatrix : public TlMatrixObject {
public:
    /// 行列オブジェクトを作成する
    ///
    /// @param[in] row 作成する行列の行数
    /// @param[in] col 作成する行列の列数
    explicit TlSparseMatrix(index_type row =1, index_type col =1);

    /// コピーコンストラクタ
    TlSparseMatrix(const TlSparseMatrix& rhs);

    /// デストラクタ
    virtual ~TlSparseMatrix();

public:
    typedef std::size_t KeyType;

#ifdef HAVE_UNORDERED_MAP
    typedef std::unordered_map<KeyType, double> SparseMatrixData;
#elifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<KeyType, double> SparseMatrixData;
#elifdef HAVE_GOOGLE_SPARSE_HASH_MAP
    typedef google::sparse_hash_map<KeyType, double> SparseMatrixData;
#else
    typedef std::map<KeyType, double> SparseMatrixData;
    #define TSM_DATATYPE_BINTREE 1
#endif

    typedef SparseMatrixData::const_iterator const_iterator;
    typedef SparseMatrixData::iterator iterator;

    // operator
public:
    /// 最初のイテレータを返す
    const_iterator begin() const {
        return m_aMatrix.begin();
    }

    /// 最初のイテレータを返す
    iterator begin() {
        return m_aMatrix.begin();
    }

    /// 最後のイテレータを返す
    const_iterator end() const {
        return m_aMatrix.end();
    }

    /// 最後のイテレータを返す
    iterator end() {
        return m_aMatrix.end();
    }

    /// オブジェクトの内容を破棄する
    virtual void clear();

    virtual void zeroClear() {
        this->clear();
    }
    
    /// 要素を削除する
    virtual void erase(index_type row, index_type col);
    virtual void erase(iterator p);

    /// 行数を返す
    ///
    /// @return 行数
    virtual index_type getNumOfRows() const;

    /// 列数を返す
    ///
    /// @return 列数
    virtual index_type getNumOfCols() const;

    /// 要素の数を返す
    ///
    /// @return 要素の数
    virtual size_type getSize() const;

    virtual std::size_t getMemSize() const;
    
    /// 行列のサイズを変更する
    ///
    /// 行列を大きくする場合、追加される要素は0で初期化される。
    /// 行列を小さくする場合、切り詰められる要素は破棄される。
    /// @param[in] row 行数
    /// @param[in] col 列数
    virtual void resize(index_type row, index_type col);

    double pop(index_type* pRow, index_type* pCol);

    /// 要素を返す(読み取り専用)
    ///
    /// 内部では、行列要素を(2次元配列ではなく)
    /// 1次元配列として保持しているので、
    /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
    ///
    /// @param[in] row 行数
    /// @param[in] col 列数
    /// @return 要素
    virtual double operator()(index_type row, index_type col) const;
    virtual double get(index_type row, index_type col) const;

    /** 要素を返す(代入可能)
     *
     *  内部では、行列要素を(2次元配列ではなく)
     *  1次元配列として保持しているので、
     *  他のメンバ関数内でもこのメンバ関数が呼ばれる。
     *
     *  @param[in] row 行数
     *  @param[in] col 列数
     *  @return 要素
     */
    virtual double& operator()(index_type row, index_type col);
    virtual void set(const index_type row, const index_type col, const double value);
    virtual void set(const std::pair<unsigned long, double>& obj);

    virtual void add(const index_type row, const index_type col, const double value);
    virtual void add(const std::pair<unsigned long, double>& obj);

    /** 指定された要素が存在すればtrueを返す
     *
     *  @retval true 要素が存在する
     *  @retval false 要素が存在しない
     */
    virtual bool hasKey(index_type row, index_type col) {
        //return (this->m_aMatrix.find(TlMatrixIndexPair(row, col)) != this->m_aMatrix.end());
        return (this->m_aMatrix.find(this->index(row, col)) != this->m_aMatrix.end());
    }

    virtual void merge(const TlSparseMatrix& rhs);

    /// 代入演算子
    TlSparseMatrix& operator =(const TlSparseMatrix& rhs);

    /// 行列を定数倍する
    ///
    /// @param[in] rhs 定数倍の値
    /// @return 計算後のオブジェクト
    TlSparseMatrix& operator*=(const double& rhs);

    /// 行列を定数で割る
    ///
    /// @param[in] rhs 割る定数の値
    /// @return 計算後のオブジェクト
    TlSparseMatrix& operator/=(const double& rhs);

    /// 指定した行の要素から構成されるベクトルを返す
    ///
    /// @param[in] row 指定する行
    virtual TlVector getRowVector(index_type row) const;

    /// 指定した列の要素から構成されるベクトルを返す
    ///
    /// @param[in] col 指定する列
    virtual TlVector getColVector(index_type col) const;

    /// Hadamard product
    const TlSparseMatrix& dot(const TlSparseMatrix& X);

    double sum() const;

    std::vector<int> getRowIndexList() const;
    std::vector<int> getColIndexList() const;

public:
    /// オブジェクトの内容をテキスト出力する
    ///
    /// @param[in] out 出力用の<< 演算子を持つオブジェクト
    template <typename T> void print(T& out) const;

public:
    unsigned long index(const index_type row, const index_type col) const;

    void index(const KeyType i, index_type* pRow, index_type* pCol) const;

public:
    virtual bool load(const std::string& path);
    virtual bool save(const std::string& path) const;

protected:
    bool load(std::ifstream& ifs);
    bool save(std::ofstream& ofs) const;
    
protected:
    index_type m_nRows;                  /// 行数
    index_type m_nCols;                  /// 列数
    mutable SparseMatrixData m_aMatrix;   /// 行列要素

    static const int INT_BITS;
    static const int MAX_INT;

    friend class TlCommunicate;
};

////////////////////////////////////////////////////////////////////////
// inline functions
inline unsigned long TlSparseMatrix::index(const index_type row, const index_type col) const
{
    const unsigned int r = static_cast<unsigned int>(row);
    const unsigned int c = static_cast<unsigned int>(col);

    return ((static_cast<unsigned long>(r) << INT_BITS) + c);
}


inline void TlSparseMatrix::index(const KeyType i, index_type* pRow, index_type* pCol) const
{
    assert(pRow != NULL);
    assert(pCol != NULL);

    const unsigned long r = i >> INT_BITS;
    const unsigned long c = i & TlSparseMatrix::MAX_INT;

    *pRow = static_cast<index_type>(r);
    *pCol = static_cast<index_type>(c);
}


inline void TlSparseMatrix::set(const index_type row, const index_type col, const double value)
{
    const unsigned long i = this->index(row, col);
    this->set(std::pair<unsigned long, double>(i, value));
}


inline void TlSparseMatrix::set(const std::pair<unsigned long, double>& obj)
{
    const unsigned long index = obj.first;

#pragma omp critical(TlSparseMatrix__update)
    {
#ifdef TSM_DATATYPE_BINTREE
        {
            iterator p = this->m_aMatrix.lower_bound(index);
            if ((p != this->m_aMatrix.end()) && (index == p->first)) {
                p->second = obj.second;
            } else {
                this->m_aMatrix.insert(p, obj);
            }
        }
#else
        {
            this->m_aMatrix[index] = obj.second;
        }
#endif // TSM_DATATYPE_BINTREE
    }
}


inline void TlSparseMatrix::add(const index_type row, const index_type col, const double value)
{
    const unsigned long index = this->index(row, col);
    this->add(std::pair<unsigned long, double>(index, value));
}


inline void TlSparseMatrix::add(const std::pair<unsigned long, double>& obj)
{
    const unsigned long index = obj.first;

#pragma omp critical(TlSparseMatrix__update)
    {
#ifdef TSM_DATATYPE_BINTREE
        {
            iterator p = this->m_aMatrix.lower_bound(index);
            if ((p != this->m_aMatrix.end()) && (index == p->first)) {
                p->second += obj.second;
            } else {
                this->m_aMatrix.insert(p, obj);
            }
        }
#else
        {
            this->m_aMatrix[index] += obj.second;
        }
#endif // TSM_DATATYPE_BINTREE
    }
}


//
inline double TlSparseMatrix::operator()(const index_type row, const index_type col) const
{
    return this->get(row, col);
}


inline double TlSparseMatrix::get(const index_type row, const index_type col) const
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));

    double answer = 0.0;
    const_iterator p = this->m_aMatrix.find(this->index(row, col));
    if (p != this->m_aMatrix.end()) {
        answer = p->second;
    }

    return answer;
}


inline double& TlSparseMatrix::operator()(const index_type row, const index_type col)
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));

    return this->m_aMatrix[this->index(row, col)];
}


template <typename T> void TlSparseMatrix::print(T& out) const
{
    const int nNumOfRows = this->getNumOfRows();
    const int nNumOfCols = this->getNumOfCols();

    for (int ord = 0; ord < nNumOfCols; ord += 10) {
        out << "       ";
        for (int j = ord; j < ord+10 && j < nNumOfCols; j++) {
            out << TlUtils::format("   %5d th", j+1);
        }
        out << "\n ----";

        for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); j++) {
            out << "-----------";
        }
        out << "----\n";

        for (int i = 0; i < nNumOfRows; i++) {
            out << TlUtils::format(" %5d  ", i+1);

            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); j++) {
                out << TlUtils::format(" %10.6lf", (*this)(i, j));
            }
            out << "\n";
        }
        out << "\n\n";
    }
    out.flush();

}

#endif // TLMATRIX_SPARSE_H

