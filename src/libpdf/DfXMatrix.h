// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DFXMATRIX_H
#define DFXMATRIX_H

#include <cassert>
#include <cmath>

#include "DfObject.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

/// X行列を求めるクラス
/// S行列から、固有値, 固有ベクトルを求め、X 行列および−1 X 行列の計算を行い、
/// X, Xinv(Xの逆行列)を出力する
class DfXMatrix : public DfObject {
public:
    DfXMatrix(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix();

public:
    virtual void buildX();

    virtual void canonicalOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                        TlDenseGeneralMatrix_Lapack* pX,
                                        TlDenseGeneralMatrix_Lapack* pXinv,
                                        const std::string& eigvalFilePath = "");

    virtual void lowdinOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                     TlDenseGeneralMatrix_Lapack* pX,
                                     TlDenseGeneralMatrix_Lapack* pXinv,
                                     const std::string& eigvalFilePath = "");

protected:
    template <typename SymmetricMatrixType, typename MatrixType, typename VectorTYpe>
    void buildX_templ();

    template <typename SymmetricMatrixType, typename MatrixType, typename VectorType>
    void canonicalOrthogonalizeTmpl(const SymmetricMatrixType& S,
                                    MatrixType* pX, MatrixType* pXinv,
                                    const std::string& eigvalFilePath = "");

    template <typename SymmetricMatrixType, typename MatrixType>
    void lowdinOrthogonalizeTmpl(const SymmetricMatrixType& S, MatrixType* pX,
                                 MatrixType* pXinv,
                                 const std::string& eigvalFilePath = "");

    template <typename MatrixType>
    void check_X(const MatrixType& X, const MatrixType& Xinv,
                 const std::string& savePathPrefix);

protected:
    /// 一次従属性の判定値(for canonical method)
    double threshold_trancation_canonical_;

    /// 一次従属性の判定値(for lowdin method)
    double threshold_trancation_lowdin_;

    /// X行列作成時に求められた固有値を保存するファイル名
    // std::string XEigvalFilePath_;

    /// 固有値を保持するかどうか
    bool debugSaveEigval_;

    /// デバッグ用に行列を保存する
    bool debugSaveMatrix_;

    /// X行列が性質を満たすかどうかをテストする(true)
    bool debugCheckX_;
};

template <typename SymmetricMatrixType, typename MatrixType, typename VectorType>
void DfXMatrix::buildX_templ() {
    SymmetricMatrixType S = this->getSpqMatrix<SymmetricMatrixType>();
    MatrixType X;
    MatrixType Xinv;

    std::string eigvalFilePath = "";
    if (this->debugSaveEigval_) {
        eigvalFilePath = DfObject::getXEigvalVtrPath();
    }
    this->canonicalOrthogonalizeTmpl<SymmetricMatrixType, MatrixType, VectorType>(S, &X, &Xinv, eigvalFilePath);

    DfObject::saveXMatrix(X);
    DfObject::saveXInvMatrix(Xinv);
    (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();
}

template <typename SymmetricMatrixType, typename MatrixType, typename VectorType>
void DfXMatrix::canonicalOrthogonalizeTmpl(const SymmetricMatrixType& S,
                                           MatrixType* pX, MatrixType* pXinv,
                                           const std::string& eigvalFilePath) {
    this->log_.info("orthogonalize by canonical method");
    this->log_.info(
        TlUtils::format("S: %d x %d", S.getNumOfRows(), S.getNumOfCols()));

    const TlMatrixObject::index_type dim = S.getNumOfRows();
    TlMatrixObject::index_type rest = 0;

    TlDenseVector_Lapack sqrt_s;  // Sの固有値の平方根
    MatrixType U;                 // Sの固有ベクトル
    {
        this->loggerTime("diagonalization of S matrix");
        VectorType EigVal;
        MatrixType EigVec;
        S.eig(&EigVal, &EigVec);
        assert(EigVal.getSize() == dim);

        if (!eigvalFilePath.empty()) {
            this->log_.info(
                TlUtils::format("save eigvals to %s", eigvalFilePath.c_str()));
            EigVal.save(eigvalFilePath);
        }

        this->loggerTime("truncation of linear dependent");
        {
            const double threshold = this->threshold_trancation_canonical_;
            this->log_.info(TlUtils::format("threshold: %f", threshold));
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
        const TlMatrixObject::index_type cutoffBasis = dim - rest;

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
    if (this->debugSaveMatrix_) {
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
        for (TlMatrixObject::index_type i = 0; i < rest; ++i) {
            S12.set(i, i, sqrt_s.get(i));
        }

        this->loggerTime("transpose U matrix");
        U.transposeInPlace();

        *pXinv = S12 * U;
    }

    this->loggerTime(" finalize");

    if (this->debugCheckX_) {
        this->check_X(*pX, *pXinv,
                      TlUtils::format("%s/S_", this->m_sWorkDirPath.c_str()));
    }
}

template <typename SymmetricMatrixType, typename MatrixType>
void DfXMatrix::lowdinOrthogonalizeTmpl(const SymmetricMatrixType& S,
                                        MatrixType* pX, MatrixType* pXinv,
                                        const std::string& eigvalFilePath) {
    this->log_.info("orthogonalize by lowdin method");
    const index_type dim = S.getNumOfRows();
    index_type rest = 0;

    TlDenseVector_Lapack sqrt_s;  // Sの固有値の平方根
    MatrixType U;                 // Sの固有ベクトル
    {
        this->loggerTime("diagonalization of S matrix");
        TlDenseVector_Lapack EigVal;
        MatrixType EigVec;
        S.eig(&EigVal, &EigVec);
        assert(EigVal.getSize() == dim);

        if (!eigvalFilePath.empty()) {
            this->log_.info(
                TlUtils::format("save eigvals to %s", eigvalFilePath.c_str()));
            EigVal.save(eigvalFilePath);
        }

        this->loggerTime("truncation of linear dependent");
        {
            const double threshold = this->threshold_trancation_lowdin_;
            this->log_.info(TlUtils::format("threshold: %f", threshold));
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
    if (this->debugSaveMatrix_) {
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

        MatrixType Ut = U;
        Ut.transposeInPlace();

        *pX *= Ut;
    }

    if (pXinv != NULL) {
        this->loggerTime("generate X^-1 matrix");
        *pXinv = *pX;
        pXinv->inverse();
    }

    this->loggerTime(" finalize");

    if (this->debugCheckX_) {
        this->check_X(*pX, *pXinv,
                      TlUtils::format("%s/S_", this->m_sWorkDirPath.c_str()));
    }
}

// X^(-1) が X の逆行列になっているかチェック
template <typename MatrixType>
void DfXMatrix::check_X(const MatrixType& X, const MatrixType& Xinv,
                        const std::string& savePathPrefix) {
    this->log_.info("check X");
    {
        const std::string path =
            TlUtils::format("XXinv.mat", savePathPrefix.c_str());
        this->log_.info(
            TlUtils::format("calc X * Xinv to save %s.", path.c_str()));
        const MatrixType XinvX = X * Xinv;
        XinvX.save(path);
    }

    {
        const std::string path =
            TlUtils::format("XinvX.mat", savePathPrefix.c_str());
        this->log_.info(
            TlUtils::format("calc invX * X to save %s.", path.c_str()));
        const MatrixType invXX = Xinv * X;
        invXX.save(path);
    }
}

#endif  // DFXMATRIX_H
