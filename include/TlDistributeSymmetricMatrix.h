#ifndef TLDISTRIBUTESYMMETRICMATRIX_H
#define TLDISTRIBUTESYMMETRICMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <list>
#include "TlDistributeMatrix.h"
#include "TlDistributeVector.h"
#include "TlSparseSymmetricMatrix.h"

class TlDistributeSymmetricMatrix : public TlDistributeMatrix {
public:
    enum DIAGONAL_METHOD {
        QR,
        DIVIDE_AND_CONQUER
    };

public:
    explicit TlDistributeSymmetricMatrix(int dim =1);
    TlDistributeSymmetricMatrix(const TlDistributeSymmetricMatrix& rhs);
    TlDistributeSymmetricMatrix(const TlDistributeMatrix& rhs);
    explicit TlDistributeSymmetricMatrix(const TlDistributeVector& rhs,
                                         int dim);

    virtual ~TlDistributeSymmetricMatrix();

public:
    void resize(int nSize);

public:
    // need to call by all processes.
    TlDistributeSymmetricMatrix& operator=(const TlDistributeSymmetricMatrix& rhs);
    TlDistributeSymmetricMatrix& operator=(const TlDistributeMatrix& rhs);

    /// 要素の絶対値の最大値を返す
    ///
    /// @param[out] outRow 該当要素の行
    /// @param[out] outCol 該当要素の列
    /// @return 要素の絶対値の最大値
    virtual double getMaxAbsoluteElement(int* pOutRow =NULL, int* pOutCol =NULL) const;

    // need to call by all processes.
    virtual double get(index_type row, index_type col) const;
    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);

    /** ローカルから該当する要素があれば値を返す
     *
     *  範囲外であれば0.0を返す
     *  すべてのプロセスが同時に呼ぶ必要はない
     */
    virtual double getLocal(index_type row, index_type col) const;

    /// 全行列を各プロセスに均等分割された疎行列を返す
    TlSparseSymmetricMatrix getPartialMatrix(double threshold = 1.0E-16) const;

    /// 指定された要素を持つ疎行列を返す
    void getPartialMatrix(TlSparseSymmetricMatrix& ioMatrix) const;
    bool getSparseMatrixX(TlSparseSymmetricMatrix* pMatrix, bool isFinalize = false) const;
    bool getPartialMatrixX(TlPartialSymmetricMatrix* pMatrix, bool isFinalize = false) const;

    TlVector getPartialMatrix(int* pStartRow, int* pEndRow, int* pStartCol, int* pEndCol) const;

    /// 各ノードが与えた疎行列を大域行列に加算する。
    /// 
    /// 全ノードがこの関数を呼び出す必要がある。
    void mergeSparseMatrix(const TlSparseSymmetricMatrix& M);

    /// 各ノードが与えた部分行列を大域行列に加算する。
    /// 
    /// 全ノードがこの関数を呼び出す必要がある。
    void mergePartialMatrix(const TlPartialSymmetricMatrix& M);
    void mergePartialMatrix(const std::list<TlPartialSymmetricMatrix>& M);

    virtual void mergeSparseMatrixAsync(const TlSparseMatrix* pMatrix,
                                        bool isFinalize = false);

    virtual double operator()(int row, int col) const;
    virtual double& operator()(int row, int col);

    /// 固有値を求める
    ///
    /// @param[out] pEigVal 固有値が格納されたベクトル
    /// @param[out] pEigVec 固有値ベクトルが格納された行列
    /// @retval true 固有値が求められた
    /// @retval false エラーが発生した
    virtual bool diagonal(TlVector* pEigVal, TlDistributeMatrix* pEigVec,
                          DIAGONAL_METHOD method = DIVIDE_AND_CONQUER);

    virtual bool inverse();

    //const TlDistributeSymmetricMatrix& dot(const TlDistributeSymmetricMatrix& X);
    virtual double sum() const;

    TlDistributeMatrix choleskyFactorization(const double threshold = 1.0E-16) const;

public:
    /// 指定された入力ストリームが読み込み可能かどうかを返す
    ///
    /// @param[in,out] ifs 入力ファイルストリーム
    /// @retval true 読み取り可能
    /// @retval false 読み取り不可能
    static bool isLoadable(std::ifstream& ifs);
    static bool isLoadable(const std::string& rFilePath);

    virtual bool load(const std::string& sFilePath);
    virtual bool load(std::ifstream& ifs);
    virtual bool save(const std::string& sFilePath) const;
    //virtual bool save(std::ofstream& ofs) const;
//   bool saveText(const std::string& sFilePath) const;
//   bool saveText(std::ofstream& ofs) const;

public:
    // need to call by all processes.
    template <typename T>
    void print(T& out) const;

protected:
    virtual void getPartialMatrixX_registerTask(TlPartialMatrix* pMatrix) const;
    
protected:
    bool load_RLHD(std::ifstream& ifs);
    bool load_CLHD(std::ifstream& ifs);

protected:
    virtual bool loadLocal(const std::string& filePath);
    virtual bool saveLocal(const std::string& filePath) const;
    
    // friend
protected:
#ifdef HAVE_SCALAPACK
    friend TlDistributeMatrix operator*(const TlDistributeSymmetricMatrix& X, const TlDistributeSymmetricMatrix& Y);
    friend TlDistributeMatrix operator*(const TlDistributeSymmetricMatrix& X, const TlDistributeMatrix& Y);
    friend TlDistributeMatrix operator*(const TlDistributeMatrix& X, const TlDistributeSymmetricMatrix& Y);
    friend TlDistributeVector operator*(const TlDistributeSymmetricMatrix& A, const TlDistributeVector& X);

    friend TlDistributeSymmetricMatrix operator*(double X, const TlDistributeSymmetricMatrix& Y);
    friend TlDistributeSymmetricMatrix operator*(const TlDistributeSymmetricMatrix& X, double Y);

    friend TlDistributeMatrix multiplicationByScaLapack(const TlDistributeMatrix& X, const TlDistributeSymmetricMatrix& Y);
    friend TlDistributeMatrix multiplicationByScaLapack(const TlDistributeSymmetricMatrix& X, const TlDistributeMatrix& Y);

    /// 対称行列の積を求める
    //friend TlMatrix multiplicationByLapack(const TlMatrix_Symmetric& X, const TlMatrix& Y);
    //friend TlMatrix multiplicationByLapack(const TlMatrix& X, const TlMatrix_Symmetric& Y);

    /// 対称行列の固有値を求める(QR法)
    ///
    /// @param[in] inMatrix 対称行列
    /// @param[out] outEigVal 固有値が格納されたベクトル
    /// @param[out] outEigVec 固有値ベクトルが格納された行列
    /// @retval true 固有値が求められた
    /// @retval false エラーが発生した
    friend bool diagonalByScaLapack_QR(TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);

    /// 対称行列の固有値を求める(Divide&Conquer)
    ///
    /// @param[in] inMatrix 対称行列
    /// @param[out] outEigVal 固有値が格納されたベクトル
    /// @param[out] outEigVec 固有値ベクトルが格納された行列
    /// @retval true 固有値が求められた
    /// @retval false エラーが発生した
    friend bool diagonalByScaLapack_DC(TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);

    friend bool inverseByScaLapack(TlDistributeSymmetricMatrix& inoutMatrix);
#endif // HAVE_SCALAPACK
};

template <typename T>
void TlDistributeSymmetricMatrix::print(T& out) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int nNumOfDim = this->getNumOfRows(); // == this->getNumOfCols()

    if (rComm.isMaster() == true) {
        out << "\n\n";
    }
    for (int ord = 0; ord < nNumOfDim; ord += 10) {
        if (rComm.isMaster() == true) {
            out << "       ";
            for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
                out << TlUtils::format("   %5d th", j+1);
            }
            out << "\n" << " ----";

            for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
                out << "-----------";
            }
            out << "----\n";
        }

        for (int i = 0; i < nNumOfDim; ++i) {
            if (rComm.isMaster() == true) {
                out << TlUtils::format(" %5d  ", i+1);
            }

            for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
                if (j > i) {
                    if (rComm.isMaster() == true) {
                        out << "    ----   ";
                    }
                } else {
                    const double dValue = this->get(i, j);
                    if (rComm.isMaster() == true) {
                        out << TlUtils::format(" %10.6lf", dValue);
                    }
                }
            }
            if (rComm.isMaster() == true) {
                out << "\n";
            }
        }
        if (rComm.isMaster() == true) {
            out << "\n\n";
        }
    }

    if (rComm.isMaster() == true) {
        out.flush();
    }
}

#endif // TLDISTRIBUTESYMMETRICMATRIX_H
