#ifndef DFXMATRIX_H
#define DFXMATRIX_H

#include "DfObject.h"
#include "TlMatrix.h"
#include "TlVector.h"

/// X行列を求めるクラス
/// S行列から、固有値, 固有ベクトルを求め、X 行列および−1 X 行列の計算を行い、
/// X, Xinv(Xの逆行列)を出力する
class DfXMatrix : public DfObject {
public:
    DfXMatrix(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix();

public:
    virtual void buildX();

    template<typename SymmetricMatrixType, typename MatrixType>
    void canonicalOrthogonalize(const SymmetricMatrixType& S,
                                MatrixType* pX, MatrixType* pXinv);

protected:
    template<typename MatrixType>
    void check_X(const MatrixType& X,
                 const MatrixType& Xinv,
                 const std::string& savePathPrefix);

protected:
    /// 一次従属性の判定値
    double threshold_trancation;

    /// デバッグ用に行列を保存する
    bool debug_save_mat_;

    /// X行列が性質を満たすかどうかをテストする(true)
    bool debug_check_X_;
};


template<typename SymmetricMatrixType, typename MatrixType>
void DfXMatrix::canonicalOrthogonalize(const SymmetricMatrixType& S,
                                       MatrixType* pX, MatrixType* pXinv)
{
    const index_type dim = S.getNumOfRows();
    index_type rest = 0;

    TlVector sqrt_s; // Sの固有値の平方根
    MatrixType U; // Sの固有ベクトル
    {
        this->loggerTime("diagonalization of S matrix");
        TlVector EigVal;
        MatrixType EigVec;
        S.diagonal(&EigVal, &EigVec);
        assert(EigVal.getSize() == dim);

        this->loggerTime("truncation of linear dependent");
        {
            const double threshold = this->threshold_trancation;
            int cutoffCount = 0;
            for (index_type k = 0; k < dim; ++k) {
                if (EigVal.get(k) < threshold) {
                    ++cutoffCount;
                } else {
                    // 固有値は小さい方から格納されているはずなので、探査終了
                    break;
                }
            }
            rest = dim - cutoffCount;
        }

        this->loggerTime(" generation of U matrix");
        const index_type cutoffBasis = dim - rest;

        {
            MatrixType trans(dim, rest);
#pragma omp parallel for
            for (index_type i = 0; i < rest; ++i) {
                trans.set(cutoffBasis + i, i, 1.0);
            }
            U = EigVec * trans;
        }
        {
            sqrt_s.resize(rest);
#pragma omp parallel for
            for (index_type k = 0; k < rest; ++k) {
                const index_type index = cutoffBasis + k;
                sqrt_s.set(k, std::sqrt(EigVal.get(index)));
            }
        }
    }
    if (this->debug_save_mat_) {
        U.save("U.mat");
        sqrt_s.save("sqrt_s.vct");
    }

    if (pX != NULL) {
        this->loggerTime("generate X matrix");
        SymmetricMatrixType S12(rest);
        for (index_type i = 0; i < rest; ++i) {
            S12.set(i, i, (1.0 / sqrt_s.get(i)));
        }

        *pX = U * S12;
    }

    if (pXinv != NULL) {
        this->loggerTime("generate X^-1 matrix");

        SymmetricMatrixType S12(rest);
        for (index_type i = 0; i < rest; ++i) {
            S12.set(i, i, sqrt_s.get(i));
        }

        U.transpose();
        *pXinv = S12 * U;
    }

    this->loggerTime(" finalize");
    
    if (this->debug_check_X_) {
        this->check_X(*pX, *pXinv,
                      TlUtils::format("%s/S_", this->m_sWorkDirPath.c_str()));
    }
}

template<typename MatrixType>
void DfXMatrix::check_X(const MatrixType& X,
                        const MatrixType& Xinv,
                        const std::string& savePathPrefix)
{
    this->log_.info("check X");

    {
        const std::string path = TlUtils::format("XXinv.mat", savePathPrefix.c_str());
        this->log_.info(TlUtils::format("calc X * Xinv to save %s.", path.c_str()));
        const MatrixType XinvX = X * Xinv;
        XinvX.save(path);
    }

    {
        const std::string path = TlUtils::format("XinvX.mat", savePathPrefix.c_str());
        this->log_.info(TlUtils::format("calc invX * X to save %s.", path.c_str()));
        const MatrixType invXX = Xinv * X;
        invXX.save(path);
    }
}

#endif // DFXMATRIX_H

