#ifndef TLVECTOR_H
#define TLVECTOR_H

#include <cassert>
#include <valarray>

#include "TlVectorObject.h"
#include "TlSparseVector.h"
#include "TlUtils.h"

class TlMatrix;
class TlSymmetricMatrix;
class TlDistributeMatrix;
class TlDistributeSymmetricMatrix;
class TlCommunicate;

/// double専用のベクトルクラス
class TlVector : public TlVectorObject {
    friend class TlCommunicate;

public:
    //typedef std::valarray<double> BufferType;

public:
    /// コンストラクタ
    ///
    /// @param[in] size ベクトルの要素数
    explicit TlVector(size_type size =0);

    /// コンストラクタ
    ///
    /// @param[in] rhs ベクトル
    TlVector(const double* p, size_type size);

    /// コンストラクタ(変換)
    TlVector(const std::vector<double>& rhs);

    /// コピーコンストラクタ
    TlVector(const TlVector& rhs);

    /// デストラクタ
    ~TlVector();

public:
    /// 要素数を返す
    ///
    /// @return 要素数
    size_type getSize() const;

    /// サイズを変更する
    ///
    /// 現在の要素数よりも大きいサイズを指定された場合は、要素数が拡張される.
    /// 現在の要素数よりも小さいサイズを指定された場合は、切り詰められる.
    /// @param[in] size 要素数
    void resize(size_type size);

    /// 末尾に値を追加する
    ///
    /// @param[in] value 末尾に追加する値
    void push_back(double value);

    /// copy
    //void copy(BufferType* pBuf) const;

    /// 要素のうち、要素の絶対値の最大を返す.
    ///
    /// @return 要素の絶対値の最大値
    double getMaxAbsoluteElement() const;

    TlVector& dot(const TlVector& rhs);
    
    /// 全要素の和を返す
    ///
    /// @return 全要素の和
    double sum() const;

    /// 要素を値の大きい順にソートする
    void sortByGrater();

    /// 大きな要素のインデックスを返す
    /// @param startIndex[in]: 探査を開始するインデックス
    std::vector<TlVectorObject::size_type>::const_iterator
    argmax(const std::vector<TlVectorObject::size_type>::const_iterator& begin,
           const std::vector<TlVectorObject::size_type>::const_iterator& end) const;

public:
    static bool isLoadable(const std::string& sFilePath);
    static bool isLoadable(std::ifstream& ifs);
    static bool getHeaderInfo(std::ifstream& ifs, int* pnSize);

    /// ファイルからデータを読み込む
    ///
    /// @param[in] sFilePath ファイルのパス名
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    virtual bool load(const std::string& sFilePath);

    /// std::ifstreamオブジェクトからデータを読み込む
    ///
    /// @param[in,out] ifs std::ifstreamオブジェクト
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    virtual bool load(std::ifstream& ifs);

    virtual bool loadText(const std::string& filePath);

    /// ファイルへデータを出力する
    ///
    /// @param[in] sFilePath ファイルのパス名
    /// @retval true 読み込み成功
    /// @retval false 読み込み失敗
    virtual bool save(const std::string& sFilePath) const;

    /// std::ofstreamオブジェクトへデータを出力する
    ///
    /// @param[out] ofs std::ofstreamオブジェクト
    /// @retval true 出力成功
    /// @retval false 出力失敗
    virtual bool save(std::ofstream& ofs) const;

    virtual void outputText(std::ostream& os) const;

public:
    /// オブジェクトの内容をテキスト出力する
    ///
    /// @param[in] out 出力用の<< 演算子を持つオブジェクト
    template<typename T> void print(T& out) const;

public:
    /// 内容を0(zero)で埋める
    void zeroClear();

    double get(size_type index) const;

    void set(size_type index, double value);

    virtual void add(const size_type index, const double value);

    /// 要素を返す(読み取り専用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    double operator[](const size_type index) const;

    /// 要素を返す(代入用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    double& operator[](const size_type index);

    /// 代入演算子
    TlVector& operator =(const TlVector& rhs);

    /// 要素数は同じである必要がある
    TlVector& operator+=(const TlVector& rhs);

    TlVector& operator+=(const TlSparseVector& rhs);

    /// 要素数は同じである必要がある
    TlVector& operator-=(const TlVector& rhs);

    /// 全要素に対し乗算する
    ///
    /// @param[in] rhs 乗算する値
    TlVector& operator*=(const double& rhs);

    /// 全要素に対し除算する
    ///
    /// @param[in] rhs 除算する値
    TlVector& operator/=(const double& rhs);

    /// ベクトル同士の加算
    friend TlVector operator+(const TlVector& X, const TlVector& Y);

    /// ベクトル同士の減算
    friend TlVector operator-(const TlVector& X, const TlVector& Y);

    /// ベクトルの定数倍
    friend TlVector operator*(const TlVector& X, const double& Y);

    /// ベクトルの除算
    friend TlVector operator*(const double& X,   const TlVector& Y) {
        return (Y * X);
    }

    /// 行列と縦ベクトルの積
    friend TlVector operator*(const TlMatrix& X, const TlVector& Y);

    /// 横ベクトルと行列の積
    friend TlVector operator*(const TlVector& Y, const TlMatrix& X);

    /// ベクトルの内積
    friend double operator*(const TlVector& X, const TlVector& Y);

    /// ベクトルの除算
    ///
    /// @param[in] X TlVectorオブジェクト
    /// @param[in] Y 除算する値
    friend TlVector operator/(const TlVector& X, const double& Y) {
        assert(std::fabs(Y) < 1.0E-20);
        return (X *(1.0/Y));
    }

protected:
    void initialize(bool isZeroClear = true);

    virtual void destroy();
    
protected:
    size_type size_;

    //BufferType m_aVector; /// 要素格納用
    double* data_; /// 要素格納用

#ifdef HAVE_LAPACK
    friend TlVector multiplicationByLapack(const TlMatrix& A, const TlVector& X);
    friend TlVector multiplicationByLapack(const TlVector& X, const TlMatrix& A);
    friend TlVector multiplicationByLapack(const TlSymmetricMatrix& A, const TlVector& X);
    friend bool diagonalByLapack(const TlSymmetricMatrix& inMatrix, TlVector* outEigVal, TlMatrix* outEigVec);
#endif // HAVE_LAPACK

#ifdef HAVE_SCALAPACK
    friend bool diagonalByScaLapack_QR(const TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);
    friend bool diagonalByScaLapack_DC(const TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);
#endif // HAVE_SCALAPACK
};

////////////////////////////////////////////////////////////////////////
// inline functions
inline TlVector::size_type TlVector::getSize() const
{
    return this->size_;
}


inline void TlVector::zeroClear()
{
    const size_type size = this->getSize();
    std::fill(this->data_, this->data_ + size, 0.0);
}


inline double TlVector::get(size_type index) const
{
    assert(index < this->getSize());
    return this->data_[index];
}


inline void TlVector::set(size_type index, double value)
{
    assert(index < this->getSize());

#pragma omp critical(TlVector__set)
    {
        this->data_[index] = value;
    }
}


inline void TlVector::add(size_type index, double value)
{
    assert(index < this->getSize());

#pragma omp atomic
    this->data_[index] += value;
}


inline double TlVector::operator[](const size_type index) const
{
    assert(index < this->getSize());

    return (this->data_[index]);
}


inline double& TlVector::operator[](const size_type index)
{
    assert(index < this->getSize());

    return (this->data_[index]);
}


template <typename T>
void TlVector::print(T& out) const
{
    const size_type size = this->getSize();
    for (size_type ord = 0; ord < size; ord += 10) {
        out << "\n";
        for (size_type j = ord; (j < ord +10) && (j < size); ++j) {
            out << TlUtils::format("   %5d th", j+1);
        }
        out << "\n";

        for (size_type j = ord; (j < ord +10) && (j < size); ++j) {
            out << "-----------";
        }
        out << "----\n\n";

        for (size_type j = ord; (j < ord +10) && (j < size); ++j) {
            out << TlUtils::format(" %10.6lf", this->get(j));
        }
        out << "\n\n";
    }

    out.flush();
}

#endif // TLVECTOR_H
