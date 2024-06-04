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

#ifndef DFINITIALGUESS_H
#define DFINITIALGUESS_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <fstream>

#include "CnError.h"
#include "DfObject.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_mmap.h"
#include "tl_matrix_utils.h"

class DfDmatrix;
class TlDenseVector_Lapack;

/// 初期電子密度を求めるクラス
class DfInitialGuess : public DfObject {
protected:
    enum CalcState {
        GUESS_UNPROCESSED = 0,
        GUESS_DONE = 1
    };

    // ------------------------------------------------------------------------
    // constructor & destructor
    // ------------------------------------------------------------------------
public:
    DfInitialGuess(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess();

    // ------------------------------------------------------------------------
    // exec
    // ------------------------------------------------------------------------
public:
    void exec();

    // ------------------------------------------------------------------------
    // state
    // ------------------------------------------------------------------------
protected:
    virtual unsigned int loadCalcState() const;
    virtual void saveCalcState(unsigned int cs);

    // ------------------------------------------------------------------------
    // create occupation
    // ------------------------------------------------------------------------
protected:
    void createOccupation();
    virtual void createOccupation(const RUN_TYPE runType);

    virtual void createOccupationByFile(const RUN_TYPE runType);

    std::vector<int> getLevel(const std::string& level);

    // ------------------------------------------------------------------------
    // create initial guess (rho)
    // ------------------------------------------------------------------------
protected:
    // virtual void createRho();

    // ------------------------------------------------------------------------
    // create initial guess (Huckel/core)
    // ------------------------------------------------------------------------
protected:
    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();

    // ------------------------------------------------------------------------
    // create initial guess (Harris)
    // ------------------------------------------------------------------------
protected:
    virtual void createInitialGuessUsingHarris();

    // ------------------------------------------------------------------------
    // create initial guess (density matrix)
    // ------------------------------------------------------------------------
protected:
    void createInitialGuessUsingDensityMatrix();

    template <typename SymmetricMatrixType, typename DfPopulationType>
    void createInitialGuessUsingDensityMatrix_tmpl();

    template <typename SymmetricMatrixType, typename DfPopulationType>
    void createInitialGuessUsingDensityMatrix(const RUN_TYPE runType);

    // ------------------------------------------------------------------------
    // create initial guess (LCAO)
    // ------------------------------------------------------------------------
protected:
    void createInitialGuessUsingLCAO();
    virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);

    /// LCAO行列を取得する
    void createLcaoByFile(const RUN_TYPE runType);

    static std::string getLcaoPath_txt(const RUN_TYPE runType);
    static std::string getLcaoPath_bin(const RUN_TYPE runType);

    void createLcaoByBinFile(const RUN_TYPE runType);
    void createLcaoByTxtFile(const RUN_TYPE runType);

    /// 行列C0を保存する
    // template <typename MatrixType>
    // void saveC0(const RUN_TYPE runType, const MatrixType& C0);

    /// Cprime0 を作成し、保存する
    // template <typename MatrixType>
    // void buildCprime0(const RUN_TYPE runType, const MatrixType& C);

    /// 密度行列を作成する
    ///
    /// DfDmatrixクラスを用いる
    void makeDensityMatrix();

    /// DfDmatrixクラスオブジェクトを作成して返す
    virtual DfDmatrix* getDfDmatrixObject(TlSerializeData* param);

    /// SCF 1回転用の密度行列を作成する
    void copyDensityMatrix();

    /// 初期値用電子密度行列を返す
    template <typename SymmetricMatrixType>
    SymmetricMatrixType getInitialDensityMatrix(const RUN_TYPE runType);

    template <typename SymmetricMatrixType, class DfPopulationType>
    SymmetricMatrixType normalizeDensityMatrix(const RUN_TYPE runType, const SymmetricMatrixType& P);

    // ------------------------------------------------------------------------
    // member variables
    // ------------------------------------------------------------------------
protected:
    bool isNormalizeDensityMatrix_;
};

// template <typename MatrixType>
// MatrixType DfInitialGuess::getLCAO(const RUN_TYPE runType)
// {
//     MatrixType lcaoMatrix;
//     const std::string binFile = TlUtils::format("./guess.lcao.%s.mat",
//                                                 this->m_sRunTypeSuffix[runType].c_str());
//     const std::string txtFile = TlUtils::format("./guess.lcao.%s.txt",
//                                                 this->m_sRunTypeSuffix[runType].c_str());
//
//     this->log_.info("get LCAO");
//     if (TlFile::isExist(binFile)) {
//         // LCAO file is prepared by binary file.
//         this->log_.info(TlUtils::format("LCAO: loading: %s",
//         binFile.c_str())); lcaoMatrix.load(binFile);
//     } else if (TlFile::isExist(txtFile)) {
//         this->log_.info(TlUtils::format("LCAO: loading: %s",
//         txtFile.c_str()));
//         // LCAO file is prepared by text file.
//         std::ifstream fi;
//         fi.open(txtFile.c_str(), std::ios::in);
//         if ((fi.rdstate() & std::ifstream::failbit) != 0) {
//             CnErr.abort(TlUtils::format("cannot open file %s.\n",
//             txtFile.c_str()));
//         }
//
//         std::string dummy_line;
//         fi >> dummy_line;
//
//         int row_dimension, col_dimension;
//         fi >> row_dimension >> col_dimension;
//         if (row_dimension != this->m_nNumOfAOs) {
//             CnErr.abort("DfInitialGuess", "", "prepare_occupation_and_or_mo",
//                         "inputted guess lcao has illegal dimension");
//         }
//         lcaoMatrix.resize(row_dimension, col_dimension);
//
//         const int maxRows = row_dimension;
//         const int maxCols = col_dimension;
//         for (int i = 0; i < maxRows; ++i) {
//             for (int j = 0; j < maxCols; ++j) {
//                 fi >> lcaoMatrix(i, j);
//             }
//         }
//     } else {
//         this->log_.warn(TlUtils::format("file not found.: %s",
//         binFile.c_str()));
//     }
//
//     return lcaoMatrix;
// }

template <typename SymmetricMatrixType, typename DfPopulationType>
void DfInitialGuess::createInitialGuessUsingDensityMatrix_tmpl() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            this->createInitialGuessUsingDensityMatrix<SymmetricMatrixType, DfPopulationType>(RUN_RKS);
        } break;

        case METHOD_UKS: {
            this->createInitialGuessUsingDensityMatrix<SymmetricMatrixType, DfPopulationType>(RUN_UKS_ALPHA);
            this->createInitialGuessUsingDensityMatrix<SymmetricMatrixType, DfPopulationType>(RUN_UKS_BETA);
        } break;

        case METHOD_ROKS: {
            this->createInitialGuessUsingDensityMatrix<SymmetricMatrixType, DfPopulationType>(RUN_ROKS_CLOSED);
            this->createInitialGuessUsingDensityMatrix<SymmetricMatrixType, DfPopulationType>(RUN_ROKS_OPEN);
        } break;

        default:
            abort();
    }

    this->createOccupation();
}

template <typename SymmetricMatrixType, typename DfPopulationType>
void DfInitialGuess::createInitialGuessUsingDensityMatrix(const RUN_TYPE runType) {
    // read guess lcao
    SymmetricMatrixType P = this->getInitialDensityMatrix<SymmetricMatrixType>(runType);
    if (this->isNormalizeDensityMatrix_) {
        P = this->normalizeDensityMatrix<SymmetricMatrixType, DfPopulationType>(runType, P);
    }
    this->savePOutMatrix(runType, 0, P);

    // spin density matrix
    const int iteration = 0;
    if (runType == RUN_RKS) {
        this->saveSpinDensityMatrix(runType, iteration, 0.5 * P);
    } else {
        this->saveSpinDensityMatrix(runType, iteration, P);
    }

    // make occupation data
    this->createOccupation();
}

template <typename SymmetricMatrixType>
SymmetricMatrixType DfInitialGuess::getInitialDensityMatrix(const RUN_TYPE runType) {
    SymmetricMatrixType lcaoMatrix;
    const std::string binFile = TlUtils::format("./guess.density.%s.mat", this->m_sRunTypeSuffix[runType].c_str());

    if (TlFile::isExistFile(binFile)) {
        lcaoMatrix.load(binFile);
    } else {
        this->log_.warn(TlUtils::format("file not found.: %s", binFile.c_str()));
    }

    return lcaoMatrix;
}

template <typename SymmetricMatrixType, class DfPopulationType>
SymmetricMatrixType DfInitialGuess::normalizeDensityMatrix(const RUN_TYPE runType, const SymmetricMatrixType& P) {
    this->log_.info(" normalize initial density matrix. ");
    DfPopulationType dfPop(this->pPdfParam_);
    const double trPS = dfPop.getSumOfElectrons(P);

    double numOfElectrons = 0.0;

    switch (runType) {
        case RUN_RKS:
            numOfElectrons = this->m_nNumOfElectrons;
            break;

        case RUN_UKS_ALPHA:
            numOfElectrons = this->m_nNumOfAlphaElectrons;
            break;

        case RUN_UKS_BETA:
            numOfElectrons = this->m_nNumOfBetaElectrons;
            break;

        case RUN_ROKS_CLOSED:
            numOfElectrons = this->numOfClosedShellElectrons_;
            break;

        case RUN_ROKS_OPEN:
            numOfElectrons = this->numOfOpenShellElectrons_;
            break;

        default:
            this->log_.critical(
                TlUtils::format("program error: %s %s", __FILE__, __LINE__));
            abort();
            break;
    }

    const double coef = numOfElectrons / trPS;

    this->log_.info(TlUtils::format(" trPS = % 8.3f", trPS));
    this->log_.info(TlUtils::format(" elec = % 8.3f", numOfElectrons));
    this->log_.info(TlUtils::format(" coef = % 8.3f", coef));

    return coef * P;
}

#endif  // DFINITIALGUESS_H
