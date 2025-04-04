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

#ifndef DFOBJECT_H
#define DFOBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <cassert>
#include <iostream>
#include <string>

#include "CnError.h"
#include "TlLogging.h"
#include "TlMatrixCache.h"
#include "TlSerializeData.h"
#include "tl_dense_vector_lapack.h"
#include "tl_matrix_object.h"
#include "tl_vector_utils.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#endif  // HAVE_LAPACK

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

/// Dfクラスの親クラス
class DfObject {
public:
    typedef signed int index_type;  // = TlMatrixObject::index_type
    typedef signed long size_type;  // = TlMatrixObject::size_type

public:
    enum METHOD_TYPE {
        METHOD_UNDEFINED = 0,
        METHOD_RKS = 1,
        METHOD_UKS,
        METHOD_ROKS
    };

    enum RUN_TYPE {
        RUN_UNDEFINED = 0,
        RUN_RKS = 1,
        RUN_UKS_ALPHA,
        RUN_UKS_BETA,
        RUN_ROKS,
        RUN_ROKS_CLOSED,
        RUN_ROKS_OPEN,
        RUN_ROKS_ALPHA,  // for XC term
        RUN_ROKS_BETA,   // for XC term
        RUN_MAXINDEX
    };

    enum GUESS_TYPE {
        GUESS_UNKNOWN = 0,
        GUESS_RHO = 1,
        GUESS_FILE_RHO,
        GUESS_CORE,
        GUESS_HUCKEL,
        GUESS_HARRIS,
        GUESS_DENSITY,
        GUESS_LCAO
    };

    struct IJShellPair {
    public:
        explicit IJShellPair(index_type i = 0, index_type j = 0)
            : nIShell(i), nJShell(j){};

        IJShellPair(const IJShellPair& rhs)
            : nIShell(rhs.nIShell), nJShell(rhs.nJShell){};

    public:
        index_type nIShell;
        index_type nJShell;
    };

public:
    DfObject(TlSerializeData* pPdfParam);
    virtual ~DfObject();

    // --------------------------------------------------------------------------
    // properties
    // --------------------------------------------------------------------------
public:
    int getNumOfAtoms() const;
    int iteration() const;

public:
    double getTotalEnergy(int iteration) const;
    double getTotalEnergy_elec(int iteration) const;

    // --------------------------------------------------------------------------
    // PATH
    // --------------------------------------------------------------------------
public:
    std::string getSpqMatrixPath();
    std::string getFpqMatrixPath(RUN_TYPE runType, int iteration) const;
    std::string getCMatrixPath(RUN_TYPE runType, int iteration,
                               const std::string& fragment = "") const;

    // density matrix
    [[deprecated]] std::string getPpqMatrixPath(RUN_TYPE runType, int iteration) const;

    std::string getPInMatrixPath(RUN_TYPE runType, int iteration) const;
    std::string getPOutMatrixPath(RUN_TYPE runType, int iteration) const;

    //
    std::string getEigenvaluesPath(const RUN_TYPE runType,
                                   const int iteration) const;

protected:
    /// 行列・ベクトルファイルのパスを返す
    ///
    /// "baseFileName"をキーとして、pdfpararm["file"]["base_base_name"]に登録されている
    /// フォーマット文字列に応じて"iteration"を置き換えてファイル名を作成する。
    std::string makeFilePath(const std::string& baseFileName,
                             const std::string& iteration = "",
                             const std::string& dir = "",
                             const std::string& suffix = "") const;
    std::string getHpqMatrixPath();
    std::string getHpq2MatrixPath();
    std::string getSabMatrixPath();
    std::string getSab2MatrixPath();
    std::string getSgdMatrixPath();
    std::string getSabInvMatrixPath();
    std::string getSgdInvMatrixPath();
    std::string getI2pqVtrPath();
    std::string getI2prVtrPath();
    std::string getI2pqVtrXCPath();
    std::string getLjkMatrixPath(const std::string& dir = "", bool isTmp = false);
    std::string getLkMatrixPath(const std::string& dir = "");
    std::string getLxcMatrixPath(const std::string& dir = "");
    std::string getLjkErrorsVtrPath();
    std::string getLkErrorsVtrPath();
    std::string getLxcErrorsVtrPath();
    std::string getXMatrixPath();
    std::string getXInvMatrixPath();
    std::string getXEigvalVtrPath();
    std::string getNalphaPath();
    std::string getOccupationPath(RUN_TYPE runType);
    std::string getGridDataFilePath() const;
    std::string getGridMatrixPath(const int iteration) const;

    // density matrix
    std::string getDiffDensityMatrixPath(RUN_TYPE runType, int iteration) const;
    std::string getSpinDensityMatrixPath(RUN_TYPE runType, int iteration) const;
    std::string getP1pqMatrixPath(int iteration);
    std::string getP2pqMatrixPath(int iteration);

    std::string getHFxMatrixPath(RUN_TYPE runType, int iteration);
    std::string getFxcMatrixPath(RUN_TYPE runType, int iteration);
    std::string getExcMatrixPath(RUN_TYPE runType, int iteration);
    std::string getFxcPureMatrixPath(RUN_TYPE runType, int iteration);
    std::string getJMatrixPath(int iteration);
    // std::string getKMatrixPath(int iteration);
    std::string getFprimeMatrixPath(RUN_TYPE runType, int iteration,
                                    const std::string& fragment = "");
    std::string getCprimeMatrixPath(RUN_TYPE runType, int iteration,
                                    const std::string& fragment = "");

    // DIIS
    std::string getDiisResidualMatrixPath(const RUN_TYPE runType, int iteration) const;

    // GridFree
    std::string getGfSMatrixPath() const;
    std::string getGfStildeMatrixPath() const;
    std::string getGfOmegaMatrixPath() const;
    std::string getGfVMatrixPath() const;
    std::string getGfVEigvalVtrPath() const;

    std::string getDipoleVelocityIntegralsXPath() const;
    std::string getDipoleVelocityIntegralsYPath() const;
    std::string getDipoleVelocityIntegralsZPath() const;

    std::string getRhoPath(RUN_TYPE nRunType, int nIteration) const;
    std::string getMyuPath(RUN_TYPE nRunType, int nIteration) const;
    std::string getNyuPath(RUN_TYPE nRunType, int nIteration) const;
    std::string getTalphaPath(RUN_TYPE runType, int iteration) const;

    // Population
    std::string getPopGrossOrbPath(RUN_TYPE runType, int iteration) const;
    std::string getPopGrossAtomPath(RUN_TYPE runType, int iteration) const;
    std::string getPopMullikenPath(RUN_TYPE runType, int iteration) const;

    // LO
    std::string getCloMatrixPath(RUN_TYPE runType, int itr) const;

protected:
    template <class SymmetricMatrixType>
    void saveHpqMatrix(const SymmetricMatrixType& Hpq);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getHpqMatrix();

    template <class SymmetricMatrixType>
    void saveHpq2Matrix(const SymmetricMatrixType& Hpq);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getHpq2Matrix();

    template <class SymmetricMatrixType>
    void saveSpqMatrix(const SymmetricMatrixType& Spq);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSpqMatrix();

    template <class MatrixType>
    void saveXInvMatrix(const MatrixType& Xinv);

    template <class MatrixType>
    MatrixType getXInvMatrix();

    template <class SymmetricMatrixType>
    void saveSabMatrix(const SymmetricMatrixType& Sab);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSabMatrix();

    template <class SymmetricMatrixType>
    void saveSab2Matrix(const SymmetricMatrixType& Sab2);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSab2Matrix();

    template <class SymmetricMatrixType>
    void saveSgdMatrix(const SymmetricMatrixType& Sab);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSgdMatrix();

    template <class SymmetricMatrixType>
    void saveSabInvMatrix(const SymmetricMatrixType& Sab);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSabInvMatrix();

    template <class SymmetricMatrixType>
    void saveSgdInvMatrix(const SymmetricMatrixType& Sgd);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSgdInvMatrix();

    template <class MatrixType>
    MatrixType getLjkMatrix(const std::string& dir = "");
    template <class MatrixType>
    void saveLjkMatrix(const MatrixType& Ljk, const std::string& dir = "");

    template <class MatrixType>
    MatrixType getLkMatrix(const std::string& dir = "");
    template <class MatrixType>
    void saveLkMatrix(const MatrixType& Lk, const std::string& dir = "");

    template <class MatrixType>
    MatrixType getLxcMatrix(const std::string& dir = "");
    template <class MatrixType>
    void saveLxcMatrix(const MatrixType& Lxc, const std::string& dir = "");

    template <class VectorType>
    VectorType loadLjkErrorsVector();
    template <class VectorType>
    void saveLjkErrorsVector(const VectorType& LjkErrors);

    template <class VectorType>
    VectorType loadLkErrorsVector();
    template <class VectorType>
    void saveLkErrorsVector(const VectorType& LkErrors);

    template <class MatrixType>
    void saveGridMatrix(const int iteration, const MatrixType& gridMatrix);

    template <class MatrixType>
    MatrixType getGridMatrix(const int iteration);

    template <class SymmetricMatrixType>
    void saveDiffDensityMatrix(const RUN_TYPE runType, const int iteration,
                               const SymmetricMatrixType& deltaP);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getDiffDensityMatrix(const RUN_TYPE runType,
                                             const int iteration) const;

    template <class SymmetricMatrixType>
    void saveSpinDensityMatrix(const RUN_TYPE runType, const int iteration,
                               const SymmetricMatrixType& P);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getSpinDensityMatrix(const RUN_TYPE runType,
                                             const int iteration);

    // GridFree S matrix -------------------------------------------------------
    template <class SymmetricMatrixType>
    void saveGfSMatrix(const SymmetricMatrixType& gfS);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getGfSMatrix();

    // GridFree S~ matrix
    // -------------------------------------------------------
    template <class MatrixType>
    void saveGfStildeMatrix(const MatrixType& gfStilde);

    template <class MatrixType>
    MatrixType getGfStildeMatrix();

    // GridFree omega matrix ---------------------------------------------------
    template <class MatrixType>
    void saveGfOmegaMatrix(const MatrixType& GfOmega);

    template <class MatrixType>
    MatrixType getGfOmegaMatrix();

    // GridFree V matrxi -------------------------------------------------------
    template <class MatrixType>
    void saveGfVMatrix(const MatrixType& gfV);

    template <class MatrixType>
    MatrixType getGfVMatrix();

    // dipole velocity integral ------------------------------------------------
    template <class MatrixType>
    void saveDipoleVelocityIntegralsXMatrix(const MatrixType& dSdx);

    template <class MatrixType>
    MatrixType getDipoleVelocityIntegralsXMatrix();

    template <class MatrixType>
    void saveDipoleVelocityIntegralsYMatrix(const MatrixType& dSdx);

    template <class MatrixType>
    MatrixType getDipoleVelocityIntegralsYMatrix();

    template <class MatrixType>
    void saveDipoleVelocityIntegralsZMatrix(const MatrixType& dSdx);

    template <class MatrixType>
    MatrixType getDipoleVelocityIntegralsZMatrix();

    // -------------------------------------------------------------------------

    template <class VectorType>
    void saveRho(const RUN_TYPE runType, const int iteration,
                 const VectorType& rho);

    template <class VectorType>
    VectorType getRho(RUN_TYPE runType, int iteration) const;

    template <class SymmetricMatrixType>
    void saveHFxMatrix(RUN_TYPE runType, int iteration,
                       const SymmetricMatrixType& HFx);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getHFxMatrix(RUN_TYPE runType, int iteration);

    template <class SymmetricMatrixType>
    void saveFxcMatrix(RUN_TYPE runType, int iteration,
                       const SymmetricMatrixType& Fxc);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getFxcMatrix(RUN_TYPE runType, int iteration);

    template <class SymmetricMatrixType>
    void saveExcMatrix(RUN_TYPE runType, int iteration,
                       const SymmetricMatrixType& Fxc);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getExcMatrix(RUN_TYPE runType, int iteration);

    template <class SymmetricMatrixType>
    void saveFxcPureMatrix(RUN_TYPE runType, int iteration,
                           const SymmetricMatrixType& FxcPure);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getFxcPureMatrix(RUN_TYPE runType, int iteration);

    template <class SymmetricMatrixType>
    void saveJMatrix(const int iteration, const SymmetricMatrixType& J);

    template <class SymmetricMatrixType>
    SymmetricMatrixType getJMatrix(int iteration);

    /** Fpq行列を保存する
     *
     */
    template <class SymmetricMatrixType>
    void saveFpqMatrix(RUN_TYPE runType, int iteration,
                       const SymmetricMatrixType& Fpq);

    /** Fpq行列を返す
     *
     */
    template <class SymmetricMatrixType>
    SymmetricMatrixType getFpqMatrix(RUN_TYPE runType, int iteration) const;

    template <class MatrixType>
    void saveXMatrix(const MatrixType& X);
    /** X行列を返す
     */
    template <class MatrixType>
    MatrixType getXMatrix();

    /** F'行列を保存する
     *
     * @param runType the type of calculation (RKS, UKS_ALPHA, UKS_BETA, ROKS)
     * @param iteration SCF iteration
     * @param fragment the name of fragment, for PDF-QCLO method
     */
    template <class SymmetricMatrixType>
    void saveFprimeMatrix(RUN_TYPE runType, int iteration,
                          const SymmetricMatrixType& Fprime,
                          const std::string& fragment = "");
    /** F'行列を返す
     *
     * @param runType the type of calculation (RKS, UKS_ALPHA, UKS_BETA, ROKS)
     * @param iteration SCF iteration
     * @param fragment the name of fragment, for PDF-QCLO method
     */
    template <class MatrixType>
    MatrixType getFprimeMatrix(RUN_TYPE runType, int iteration,
                               const std::string& fragment = "");

    template <class MatrixType>
    void saveCprimeMatrix(RUN_TYPE runType, int iteration,
                          const std::string& fragment,
                          const MatrixType& Cprime);

    template <class MatrixType>
    MatrixType getCprimeMatrix(RUN_TYPE runType, int iteration,
                               const std::string& fragment = "");

    template <class MatrixType>
    void saveCMatrix(RUN_TYPE nRunType, int nIteration, const MatrixType& C);

    /** C(LCAO)行列を返す
     *
     * @param runType the type of calculation (RKS, UKS_ALPHA, UKS_BETA, ROKS)
     * @param iteration SCF iteration
     * @param fragment the name of fragment, for PDF-QCLO method
     */
    template <class MatrixType>
    MatrixType getCMatrix(RUN_TYPE runType, int iteration,
                          const std::string& fragment = "");

    template <class SymmetricMatrixType>
    [[deprecated]] void savePpqMatrix(const RUN_TYPE runType, const int iteration,
                                      const SymmetricMatrixType& Ppq);
    template <class SymmetricMatrixType>
    [[deprecated]] SymmetricMatrixType getPpqMatrix(RUN_TYPE runType, int iteration) const;

    template <class SymmetricMatrixType>
    void savePInMatrix(const RUN_TYPE runType, const int iteration,
                       const SymmetricMatrixType& Ppq);
    template <class SymmetricMatrixType>
    SymmetricMatrixType getPInMatrix(RUN_TYPE runType, int iteration) const;

    template <class SymmetricMatrixType>
    void savePOutMatrix(const RUN_TYPE runType, const int iteration,
                        const SymmetricMatrixType& Ppq);
    template <class SymmetricMatrixType>
    SymmetricMatrixType getPOutMatrix(RUN_TYPE runType, int iteration) const;

    // template<class SymmetricMatrixType>
    // void savePCMatrix(const int iteration,
    //                   const SymmetricMatrixType& PC);
    // template <class SymmetricMatrixType>
    // SymmetricMatrixType getPCMatrix(int iteration);
    // template<class SymmetricMatrixType>
    // void savePOMatrix(const int iteration,
    //                   const SymmetricMatrixType& PO);
    // template <class SymmetricMatrixType>
    // SymmetricMatrixType getPOMatrix(int iteration);

    template <class VectorType>
    void saveNalpha(const VectorType& Na);
    template <class VectorType>
    VectorType getNalpha();

    template <class VectorType>
    void saveMyu(const RUN_TYPE runType, const int iteration,
                 const VectorType& myu);
    template <class VectorType>
    VectorType getMyu(RUN_TYPE nRunType, int nIteration);

    template <class VectorType>
    void saveNyu(const RUN_TYPE runType, const int iteration,
                 const VectorType& nyu);
    template <class VectorType>
    VectorType getNyu(RUN_TYPE runType, int iteration);

    // for LO
    template <class MatrixType>
    void saveCloMatrix(const RUN_TYPE runType, const int itr,
                       const MatrixType& Clo);
    template <class MatrixType>
    MatrixType getCloMatrix(RUN_TYPE runType, int itr);

    /// return occupation vector
    template <typename Vector>
    Vector getOccVtr(const RUN_TYPE runType);

    // --------------------------------------------------------------------------
    // methods
    // --------------------------------------------------------------------------
protected:
    virtual void setParam(const TlSerializeData& data);
    void updateLinearAlgebraPackageParam(const std::string& keyword);

protected:
    // void clearCache(const std::string& path);

    // --------------------------------------------------------------------------
    // logger
    // --------------------------------------------------------------------------
protected:
    virtual void logger(const std::string& str) const;
    void loggerTime(const std::string& str) const;
    void loggerStartTitle(const std::string& stepName,
                          const char lineChar = '-') const;
    void loggerEndTitle(const std::string& stepName = "",
                        const char lineChar = '-') const;

    // --------------------------------------------------------------------------
    // constants
    // --------------------------------------------------------------------------
protected:
    enum J_Engine_Type { J_ENGINE_CONVENTIONAL,
                         J_ENGINE_RI_J,
                         J_ENGINE_CD };

    enum K_Engine_Type {
        K_ENGINE_CONVENTIONAL,
        K_ENGINE_RI_K,
        K_ENGINE_CD,
        K_ENGINE_FASTCDK
    };

    enum XC_Engine_Type {
        XC_ENGINE_GRID,
        XC_ENGINE_GRIDFREE,
        XC_ENGINE_GRIDFREE_CD
    };

    enum LinearAlgebraPackageType {
        LAP_LAPACK,
        LAP_EIGEN,
        LAP_VIENNACL,
        LAP_SCALAPACK
    };

    // --------------------------------------------------------------------------
    // parameters
    // --------------------------------------------------------------------------
protected:
    static const std::string m_sWorkDirPath;  // fl_Work directory name
    static const std::string m_sTempDirPath;  // fl_Work directory name // before fl_Temp
    static const std::string m_sRunTypeSuffix[RUN_MAXINDEX];

protected:
    TlSerializeData* pPdfParam_;

    TlLogging& log_;

    // system parameters -------------------------------------------------------
    /// プロセスあたりの最大メモリ容量(byte)
    // std::size_t procMaxMemSize_;

    /// OpenMPスレッド数
    int numOfThreads_;

    // bool isWorkOnDisk_;
    std::string localTempPath_;
    bool useMmap_;

    bool isRestart_;

    bool isUseNewEngine_;

    METHOD_TYPE m_nMethodType;
    int m_nIteration;

    int m_nNumOfAtoms;
    int m_nNumOfDummyAtoms;
    // int numOfRealAtoms_;
    index_type m_nNumOfAOs;
    index_type m_nNumOfMOs;
    size_type m_nNumOfAux;
    size_type numOfAuxXC_;

    int m_nNumOfElectrons;
    int m_nNumOfAlphaElectrons;
    int m_nNumOfBetaElectrons;
    int numOfClosedShellElectrons_;
    int numOfOpenShellElectrons_;

    GUESS_TYPE initialGuessType_;

    int chargeExtrapolateNumber_;
    bool m_bDiskUtilization;
    bool m_bMemorySave;

    /// update法を使用する(true)
    bool isUpdateMethod_;

    // J
    J_Engine_Type J_engine_;

    // K
    K_Engine_Type K_engine_;

    // XC
    XC_Engine_Type XC_engine_;

    std::string m_sXCFunctional;
    bool isDFT_;          /// XC項にpure DFTを含む場合はtrue
    bool m_bIsXCFitting;  /// true => XC項をRI法で計算, false => XC項を直接計算
    bool m_bIsUpdateXC;   /// XC項をupdate法で計算する(true)
    bool
        enableGrimmeDispersion_;  /// Grimmeの経験的分散力補正を計算するかどうか

    bool
        isDedicatedBasisForGridFree_;  /// Grid-Free法を使用する場合、専用の基底関数を使用する

    // bool isRI_J_; /// RI_J法を用いる(true)
    bool isRI_K_;  /// RI-K法を用いる(true)

    //
    LinearAlgebraPackageType linearAlgebraPackage_;

    /// ScaLAPACKを使用する(true)かどうかを表すフラグ (deprecated)
    bool m_bUsingSCALAPACK;

    /// ScaLAPACKのブロックサイズを指定する
    int scalapackBlockSize_;

    /// 分散行列をローカルディスクへ保存する(true)かどうかを表すフラグ
    /// falseの場合は分散行列でも1つのファイルに保存される。
    /// default value is false
    // [TODO]
    // bool isSaveDistributedMatrixToLocalDisk_;

    // 一時ファイル保存に利用できるローカルディスクのパス
    // std::string localDiskPath_;

    // for HPC ==== ====================================================
    static int rank_;  /// the number of process
    int grainSize_;
    bool isMasterSlave_;  /// Master-Slave方式のときtrue

    // for debug =======================================================
    bool isFileWarning;            /// fileが無い場合のメッセージを表示する
    bool isSaveJMatrix_;           /// J行列を保存する
    bool enableExperimentalCode_;  /// 実験コードを有効にする

    // matrix memory cache
    // bool isUseCache_;
    static int objectCount_;  // DfObjectの存在個数(派生クラスを含む)
    static TlMatrixCache matrixCache_;
};

// template ============================================================
template <class SymmetricMatrixType>
void DfObject::saveHpqMatrix(const SymmetricMatrixType& Hpq) {
    const std::string path = this->getHpqMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Hpq, true);
    // } else {
    Hpq.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getHpqMatrix() {
    SymmetricMatrixType Hpq;
    const std::string path = this->getHpqMatrixPath();
    // Hpq = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Hpq.resize((this->m_nNumOfAOs));
    Hpq.load(path);

    return Hpq;
}

template <class SymmetricMatrixType>
void DfObject::saveHpq2Matrix(const SymmetricMatrixType& Hpq2) {
    const std::string path = this->getHpq2MatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Hpq2, true);
    // } else {
    Hpq2.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getHpq2Matrix() {
    SymmetricMatrixType Hpq2;
    const std::string path = this->getHpq2MatrixPath();
    // Hpq2 = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Hpq2.resize((this->m_nNumOfAOs));
    Hpq2.load(path);

    return Hpq2;
}

template <class SymmetricMatrixType>
void DfObject::saveSpqMatrix(const SymmetricMatrixType& Spq) {
    const std::string path = this->getSpqMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Spq, true);
    // } else {
    Spq.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSpqMatrix() {
    SymmetricMatrixType Spq;
    const std::string path = this->getSpqMatrixPath();
    // Spq = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Spq.resize(this->m_nNumOfAOs);
    Spq.load(path);
    return Spq;
}

template <class MatrixType>
MatrixType DfObject::getLjkMatrix(const std::string& dir) {
    const std::string path = this->getLjkMatrixPath(dir);
    MatrixType Ljk;
    Ljk.load(path, 0);

    return Ljk;
}

template <class MatrixType>
void DfObject::saveLjkMatrix(const MatrixType& Ljk, const std::string& dir) {
    const std::string path = this->getLjkMatrixPath(dir);
    Ljk.save(path, 0);
}

template <class MatrixType>
MatrixType DfObject::getLkMatrix(const std::string& dir) {
    const std::string path = this->getLkMatrixPath(dir);
    MatrixType Lk;
    Lk.load(path, 0);

    return Lk;
}

template <class MatrixType>
void DfObject::saveLkMatrix(const MatrixType& Lk, const std::string& dir) {
    const std::string path = this->getLkMatrixPath();
    Lk.save(path, 0);
}

template <class MatrixType>
MatrixType DfObject::getLxcMatrix(const std::string& dir) {
    const std::string path = this->getLxcMatrixPath();
    MatrixType Lxc;
    Lxc.load(path, 0);

    return Lxc;
}

template <class MatrixType>
void DfObject::saveLxcMatrix(const MatrixType& Lxc, const std::string& dir) {
    const std::string path = this->getLxcMatrixPath(dir);
    Lxc.save(path);
}

template <class VectorType>
VectorType DfObject::loadLjkErrorsVector() {
    const std::string path = this->getLjkErrorsVtrPath();
    VectorType Lerrors;
    Lerrors.load(path);

    return Lerrors;
}

template <class VectorType>
void DfObject::saveLjkErrorsVector(const VectorType& Lerrors) {
    const std::string path = this->getLjkErrorsVtrPath();
    Lerrors.save(path);
}

template <class VectorType>
VectorType DfObject::loadLkErrorsVector() {
    const std::string path = this->getLkErrorsVtrPath();
    VectorType Lerrors;
    Lerrors.load(path);

    return Lerrors;
}

template <class VectorType>
void DfObject::saveLkErrorsVector(const VectorType& Lerrors) {
    const std::string path = this->getLkErrorsVtrPath();
    Lerrors.save(path);
}

template <class MatrixType>
void DfObject::saveXInvMatrix(const MatrixType& XInv) {
    const std::string path = this->getXInvMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, XInv, true);
    // } else {
    XInv.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getXInvMatrix() {
    MatrixType XInv;
    const std::string path = this->getXInvMatrixPath();

    // XInv = this->matrixCache_.get<MatrixType>(path);
    if (TlFile::isExistFile(path)) {
        XInv.load(path);
    } else {
        CnErr.abort("XInv matrix file is not found: " + path);
    }

    return XInv;
}

template <class SymmetricMatrixType>
void DfObject::saveSabMatrix(const SymmetricMatrixType& Sab) {
    const std::string path = this->getSabMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Sab, true);
    // } else {
    Sab.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSabMatrix() {
    SymmetricMatrixType Sab;
    const std::string path = this->getSabMatrixPath();
    // Sab = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Sab.resize(this->m_nNumOfAux);
    Sab.load(path);

    return Sab;
}

template <class SymmetricMatrixType>
void DfObject::saveSab2Matrix(const SymmetricMatrixType& Sab2) {
    const std::string path = this->getSab2MatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Sab2, true);
    // } else {
    Sab2.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSab2Matrix() {
    SymmetricMatrixType Sab2;
    const std::string path = this->getSab2MatrixPath();
    // Sab2 = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Sab2.resize(this->m_nNumOfAux);
    Sab2.load(path);

    return Sab2;
}

template <class SymmetricMatrixType>
void DfObject::saveSgdMatrix(const SymmetricMatrixType& Sgd) {
    const std::string path = this->getSgdMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Sgd, true);
    // } else {
    Sgd.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSgdMatrix() {
    SymmetricMatrixType Sgd;
    const std::string path = this->getSgdMatrixPath();
    // Sgd = this->matrixCache_.get<SymmetricMatrixType>(path);
    // Sgd.resize(this->m_nNumOfAux);
    Sgd.load(path);

    return Sgd;
}

template <class SymmetricMatrixType>
void DfObject::saveSabInvMatrix(const SymmetricMatrixType& SabInv) {
    const std::string path = this->getSabInvMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, SabInv, true);
    // } else {
    SabInv.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSabInvMatrix() {
    SymmetricMatrixType SabInv;
    const std::string path = this->getSabInvMatrixPath();
    // SabInv = this->matrixCache_.get<SymmetricMatrixType>(path);
    // SabInv.resize(this->m_nNumOfAux);
    SabInv.load(path);

    return SabInv;
}

template <class SymmetricMatrixType>
void DfObject::saveSgdInvMatrix(const SymmetricMatrixType& SgdInv) {
    const std::string path = this->getSgdInvMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, SgdInv, true);
    // } else {
    SgdInv.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSgdInvMatrix() {
    SymmetricMatrixType SgdInv;
    const std::string path = this->getSgdInvMatrixPath();
    // SgdInv = this->matrixCache_.get<SymmetricMatrixType>(path);
    // SgdInv.resize(this->m_nNumOfAux);
    SgdInv.load(path);

    return SgdInv;
}

template <class MatrixType>
void DfObject::saveGridMatrix(const int iteration,
                              const MatrixType& gridMatrix) {
    const std::string path = this->getGridMatrixPath(iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, gridMatrix, true);
    // } else {
    gridMatrix.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getGridMatrix(const int iteration) {
    MatrixType gridMatrix;
    const std::string path = this->getGridMatrixPath(iteration);
    // gridMatrix = this->matrixCache_.get<MatrixType>(path);
    gridMatrix.load(path);

    return gridMatrix;
}

// diff density matrix
template <class SymmetricMatrixType>
void DfObject::saveDiffDensityMatrix(const RUN_TYPE runType,
                                     const int iteration,
                                     const SymmetricMatrixType& deltaP) {
    const std::string path = this->getDiffDensityMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, deltaP, true);
    // } else {
    deltaP.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getDiffDensityMatrix(const RUN_TYPE runType,
                                                   const int iteration) const {
    SymmetricMatrixType deltaP;
    const std::string path = this->getDiffDensityMatrixPath(runType, iteration);

    // deltaP = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        deltaP.load(path);
    } else {
        CnErr.abort("DiffDensity matrix file is not found: " + path);
    }

    return deltaP;
}

// spin density matrix
template <class SymmetricMatrixType>
void DfObject::saveSpinDensityMatrix(const RUN_TYPE runType,
                                     const int iteration,
                                     const SymmetricMatrixType& P) {
    const std::string path = this->getSpinDensityMatrixPath(runType, iteration);
    P.save(path);
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getSpinDensityMatrix(const RUN_TYPE runType,
                                                   const int iteration) {
    SymmetricMatrixType P;
    const std::string path = this->getSpinDensityMatrixPath(runType, iteration);

    if (TlFile::isExistFile(path)) {
        P.load(path);
    } else {
        CnErr.abort("SpinDensity matrix file is not found: " + path);
    }

    return P;
}

// rho
template <class VectorType>
void DfObject::saveRho(const RUN_TYPE runType, const int iteration,
                       const VectorType& rho) {
    const std::string path = this->getRhoPath(runType, iteration);
    rho.save(path);
}

template <class VectorType>
VectorType DfObject::getRho(const RUN_TYPE runType, const int iteration) const {
    VectorType rho;
    const std::string path = this->getRhoPath(runType, iteration);

    if (TlFile::isExistFile(path)) {
        rho.load(path);
    } else {
        CnErr.abort("Rho file is not found: " + path);
    }
    // rho.resize(this->m_nNumOfAux);

    return rho;
}

template <class SymmetricMatrixType>
void DfObject::saveHFxMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& HFx) {
    const std::string path = this->getHFxMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, HFx, true);
    // } else {
    HFx.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getHFxMatrix(const RUN_TYPE runType,
                                           const int iteration) {
    SymmetricMatrixType HFx;
    const std::string path = this->getHFxMatrixPath(runType, iteration);

    // HFx = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        HFx.load(path);
    } else {
        CnErr.abort("HFx matrix file is not found: " + path);
    }
    // HFx.resize(this->m_nNumOfAOs);

    return HFx;
}

template <class SymmetricMatrixType>
void DfObject::saveFxcMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& Fxc) {
    const std::string path = this->getFxcMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Fxc, true);
    // } else {
    Fxc.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getFxcMatrix(const RUN_TYPE runType,
                                           const int iteration) {
    SymmetricMatrixType Fxc;
    const std::string path = this->getFxcMatrixPath(runType, iteration);

    // Fxc = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Fxc.load(path);
    } else {
        CnErr.abort("Fxc matrix file is not found: " + path);
    }
    // Fxc.resize(this->m_nNumOfAOs);

    return Fxc;
}

template <class SymmetricMatrixType>
void DfObject::saveExcMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& Exc) {
    const std::string path = this->getExcMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Exc, true);
    // } else {
    Exc.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getExcMatrix(const RUN_TYPE runType,
                                           const int iteration) {
    SymmetricMatrixType Exc;
    const std::string path = this->getExcMatrixPath(runType, iteration);

    // Exc = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Exc.load(path);
    } else {
        CnErr.abort("Exc matrix file is not found: " + path);
    }

    return Exc;
}

template <class SymmetricMatrixType>
void DfObject::saveFxcPureMatrix(const RUN_TYPE runType, const int iteration,
                                 const SymmetricMatrixType& FxcPure) {
    const std::string path = this->getFxcPureMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, FxcPure, true);
    // } else {
    FxcPure.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getFxcPureMatrix(const RUN_TYPE runType,
                                               const int iteration) {
    SymmetricMatrixType FxcPure;
    const std::string path = this->getFxcPureMatrixPath(runType, iteration);
    // FxcPure = this->matrixCache_.get<SymmetricMatrixType>(path);
    // FxcPure.resize(this->m_nNumOfAOs);
    FxcPure.load(path);

    return FxcPure;
}

template <class SymmetricMatrixType>
void DfObject::saveJMatrix(const int iteration, const SymmetricMatrixType& J) {
    const std::string path = this->getJMatrixPath(iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, J, true);
    // } else {
    J.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getJMatrix(const int iteration) {
    SymmetricMatrixType J;
    const std::string path = this->getJMatrixPath(iteration);

    // J = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        J.load(path);
    } else {
        CnErr.abort("J matrix file is not found: " + path);
    }
    // J.resize(this->m_nNumOfAOs);

    return J;
}

// template<class SymmetricMatrixType>
// void DfObject::saveKMatrix(const RUN_TYPE runType, const int iteration,
//                            const SymmetricMatrixType& K)
// {
//     const std::string path = this->getKMatrixPath(runType, iteration);
//     if (this->isUseCache_ == true) {
//         this->matrixCache_.set(path, K, true);
//     } else {
//         K.save(path);
//     }
// }

// template<class SymmetricMatrixType>
// SymmetricMatrixType DfObject::getKMatrix(const RUN_TYPE runType,
//                                          const int iteration)
// {
//     SymmetricMatrixType K;
//     const std::string path = this->getKMatrixPath(runType, iteration);
//     K = this->matrixCache_.get<SymmetricMatrixType>(path);
//     K.resize(this->m_nNumOfAOs);
//     return K;
// }

template <class SymmetricMatrixType>
void DfObject::saveFpqMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& Fpq) {
    const std::string path = this->getFpqMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Fpq, true);
    // } else {
    Fpq.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getFpqMatrix(const RUN_TYPE runType,
                                           const int iteration) const {
    SymmetricMatrixType Fpq;
    const std::string path = this->getFpqMatrixPath(runType, iteration);

    // Fpq = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Fpq.load(path);
    } else {
        CnErr.abort("Fpq matrix file is not found: " + path);
    }
    // Fpq.resize(this->m_nNumOfAOs);

    return Fpq;
}

template <class MatrixType>
void DfObject::saveXMatrix(const MatrixType& X) {
    const std::string path = this->getXMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, X, true);
    // } else {
    X.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getXMatrix() {
    MatrixType X;
    const std::string path = this->getXMatrixPath();

    // X = this->matrixCache_.get<MatrixType>(path);
    if (TlFile::isExistFile(path)) {
        X.load(path);
    } else {
        CnErr.abort("X matrix file is not found: " + path);
    }

    return X;
}

template <class SymmetricMatrixType>
void DfObject::saveFprimeMatrix(RUN_TYPE runType, int iteration,
                                const SymmetricMatrixType& Fprime,
                                const std::string& fragment) {
    const std::string path =
        this->getFprimeMatrixPath(runType, iteration, fragment);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Fprime, true);
    // } else {
    Fprime.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getFprimeMatrix(const RUN_TYPE runType,
                                              const int iteration,
                                              const std::string& fragment) {
    SymmetricMatrixType Fprime;
    const std::string path = this->getFprimeMatrixPath(runType, iteration, fragment);

    // Fprime = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Fprime.load(path);
    } else {
        CnErr.abort("Fprime matrix file is not found: " + path);
    }

    return Fprime;
}

template <class MatrixType>
void DfObject::saveCprimeMatrix(const RUN_TYPE runType, const int iteration,
                                const std::string& fragment,
                                const MatrixType& Cprime) {
    const std::string path = this->getCprimeMatrixPath(runType, iteration, fragment);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Cprime, true);
    // } else {
    Cprime.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getCprimeMatrix(const RUN_TYPE runType,
                                     const int iteration,
                                     const std::string& fragment) {
    MatrixType Cprime;
    const std::string path = this->getCprimeMatrixPath(runType, iteration, fragment);

    // Cprime = this->matrixCache_.get<MatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Cprime.load(path);
    } else {
        CnErr.abort("Cprime matrix file is not found: " + path);
    }

    return Cprime;
}

template <class MatrixType>
void DfObject::saveCMatrix(const RUN_TYPE runType, const int iteration,
                           const MatrixType& C) {
    const std::string path = this->getCMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, C, true);
    // } else {
    C.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getCMatrix(const RUN_TYPE runType, const int iteration,
                                const std::string& fragment) {
    MatrixType C;
    const std::string path = this->getCMatrixPath(runType, iteration, "");

    // C = this->matrixCache_.get<MatrixType>(path);
    if (TlFile::isExistFile(path)) {
        C.load(path);
    } else {
        CnErr.abort("C matrix file is not found: " + path);
    }

    return C;
}

template <class SymmetricMatrixType>
void DfObject::savePpqMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& Ppq) {
    const std::string path = this->getPpqMatrixPath(runType, iteration);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Ppq, true);
    // } else {
    Ppq.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getPpqMatrix(const RUN_TYPE runType,
                                           const int iteration) const {
    SymmetricMatrixType Ppq;
    const std::string path = this->getPpqMatrixPath(runType, iteration);

    // Ppq = this->matrixCache_.get<SymmetricMatrixType>(path);
    if (TlFile::isExistFile(path)) {
        Ppq.load(path);
    } else {
        CnErr.abort("Ppq matrix file is not found: " + path);
    }

    return Ppq;
}

// density matrix (input)
template <class SymmetricMatrixType>
void DfObject::savePInMatrix(const RUN_TYPE runType, const int iteration,
                             const SymmetricMatrixType& P) {
    const std::string path = this->getPInMatrixPath(runType, iteration);
    P.save(path);
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getPInMatrix(const RUN_TYPE runType,
                                           const int iteration) const {
    SymmetricMatrixType P;
    const std::string path = this->getPInMatrixPath(runType, iteration);

    if (TlFile::isExistFile(path)) {
        P.load(path);
    } else {
        CnErr.abort("PIn matrix file is not found: " + path);
    }

    return P;
}

// density matrix (output)
template <class SymmetricMatrixType>
void DfObject::savePOutMatrix(const RUN_TYPE runType, const int iteration,
                              const SymmetricMatrixType& P) {
    const std::string path = this->getPOutMatrixPath(runType, iteration);
    P.save(path);
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getPOutMatrix(const RUN_TYPE runType,
                                            const int iteration) const {
    SymmetricMatrixType P;
    const std::string path = this->getPOutMatrixPath(runType, iteration);
    if (TlFile::isExistFile(path)) {
        P.load(path);
    } else {
        CnErr.abort("POut matrix file is not found: " + path);
    }

    return P;
}

// template<class SymmetricMatrixType>
// void DfObject::savePCMatrix(const int iteration,
//                             const SymmetricMatrixType& PC)
// {
//     const std::string path = this->getP1pqMatrixPath(iteration);
//     if (this->isUseCache_ == true) {
//         this->matrixCache_.set(path, PC, true);
//     } else {
//         PC.save(path);
//     }
// }

// template<class SymmetricMatrixType>
// SymmetricMatrixType DfObject::getPCMatrix(const int iteration)
// {
//     SymmetricMatrixType PC;
//     const std::string path = this->getP1pqMatrixPath(iteration);
//     PC = this->matrixCache_.get<SymmetricMatrixType>(path);
//     PC.resize(this->m_nNumOfAOs);
//     return PC;
// }

// template<class SymmetricMatrixType>
// void DfObject::savePOMatrix(const int iteration,
//                             const SymmetricMatrixType& PO)
// {
//     const std::string path = this->getP2pqMatrixPath(iteration);
//     if (this->isUseCache_ == true) {
//         this->matrixCache_.set(path, PO, true);
//     } else {
//         PO.save(path);
//     }
// }

// template<class SymmetricMatrixType>
// SymmetricMatrixType DfObject::getPOMatrix(const int iteration)
// {
//     SymmetricMatrixType PO;
//     const std::string path = this->getP2pqMatrixPath(iteration);
//     PO = this->matrixCache_.get<SymmetricMatrixType>(path);
//     PO.resize(this->m_nNumOfAOs);
//     return PO;
// }

template <class SymmetricMatrixType>
void DfObject::saveGfSMatrix(const SymmetricMatrixType& gfS) {
    const std::string path = this->getGfSMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, gfS, true);
    // } else {
    gfS.save(path);
    // }
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfObject::getGfSMatrix() {
    SymmetricMatrixType gfS;
    const std::string path = this->getGfSMatrixPath();
    // gfS = this->matrixCache_.get<SymmetricMatrixType>(path);
    gfS.load(path);

    return gfS;
}

template <class MatrixType>
void DfObject::saveGfStildeMatrix(const MatrixType& gfStilde) {
    const std::string path = this->getGfStildeMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, gfStilde, true);
    // } else {
    gfStilde.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getGfStildeMatrix() {
    MatrixType gfStilde;
    const std::string path = this->getGfStildeMatrixPath();
    // gfStilde = this->matrixCache_.get<MatrixType>(path);
    gfStilde.load(path);

    return gfStilde;
}

template <class MatrixType>
void DfObject::saveGfOmegaMatrix(const MatrixType& gfOmega) {
    const std::string path = this->getGfOmegaMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, gfOmega, true);
    // } else {
    gfOmega.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getGfOmegaMatrix() {
    MatrixType gfOmega;
    const std::string path = this->getGfOmegaMatrixPath();
    // gfOmega = this->matrixCache_.get<MatrixType>(path);
    gfOmega.load(path);

    return gfOmega;
}

template <class MatrixType>
void DfObject::saveGfVMatrix(const MatrixType& gfV) {
    const std::string path = this->getGfVMatrixPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, gfV, true);
    // } else {
    gfV.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getGfVMatrix() {
    MatrixType gfV;
    const std::string path = this->getGfVMatrixPath();
    // gfV = this->matrixCache_.get<MatrixType>(path);
    gfV.load(path);

    return gfV;
}

template <class MatrixType>
void DfObject::saveDipoleVelocityIntegralsXMatrix(const MatrixType& dSdx) {
    const std::string path = this->getDipoleVelocityIntegralsXPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, dSdx, true);
    // } else {
    dSdx.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getDipoleVelocityIntegralsXMatrix() {
    MatrixType dSdx;
    const std::string path = this->getDipoleVelocityIntegralsXPath();
    // dSdx = this->matrixCache_.get<MatrixType>(path);
    dSdx.load(path);

    return dSdx;
}

template <class MatrixType>
void DfObject::saveDipoleVelocityIntegralsYMatrix(const MatrixType& dSdy) {
    const std::string path = this->getDipoleVelocityIntegralsYPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, dSdy, true);
    // } else {
    dSdy.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getDipoleVelocityIntegralsYMatrix() {
    MatrixType dSdy;
    const std::string path = this->getDipoleVelocityIntegralsYPath();
    // dSdy = this->matrixCache_.get<MatrixType>(path);
    dSdy.load(path);

    return dSdy;
}

template <class MatrixType>
void DfObject::saveDipoleVelocityIntegralsZMatrix(const MatrixType& dSdz) {
    const std::string path = this->getDipoleVelocityIntegralsZPath();
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, dSdz, true);
    // } else {
    dSdz.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getDipoleVelocityIntegralsZMatrix() {
    MatrixType dSdz;
    const std::string path = this->getDipoleVelocityIntegralsZPath();
    // dSdz = this->matrixCache_.get<MatrixType>(path);
    dSdz.load(path);

    return dSdz;
}

template <class VectorType>
void DfObject::saveNalpha(const VectorType& Na) {
    const std::string path = this->getNalphaPath();
    Na.save(path);
}

template <class VectorType>
VectorType DfObject::getNalpha() {
    VectorType Na(this->m_nNumOfAux);
    const std::string path = this->getNalphaPath();
    if (TlVectorUtils::isLoadable(path) == true) {
        Na.load(path);
    }
    return Na;
}

template <class VectorType>
void DfObject::saveMyu(const RUN_TYPE runType, const int iteration,
                       const VectorType& myu) {
    const std::string path = this->getMyuPath(runType, iteration);
    myu.save(path);
}

template <class VectorType>
VectorType DfObject::getMyu(const RUN_TYPE runType, const int iteration) {
    VectorType myu(this->numOfAuxXC_);
    const std::string path = this->getMyuPath(runType, iteration);
    if (TlVectorUtils::isLoadable(path) == true) {
        myu.load(path);
    }
    return myu;
}

template <class VectorType>
void DfObject::saveNyu(const RUN_TYPE runType, const int iteration,
                       const VectorType& nyu) {
    const std::string path = this->getNyuPath(runType, iteration);
    nyu.save(path);
}

template <class VectorType>
VectorType DfObject::getNyu(const RUN_TYPE runType, const int iteration) {
    VectorType nyu(this->numOfAuxXC_);
    const std::string path = this->getNyuPath(runType, iteration);
    if (TlVectorUtils::isLoadable(path) == true) {
        nyu.load(path);
    }
    return nyu;
}

template <class MatrixType>
void DfObject::saveCloMatrix(const RUN_TYPE runType, const int itr,
                             const MatrixType& Clo) {
    const std::string path = this->getCloMatrixPath(runType, itr);
    // if (this->isUseCache_ == true) {
    //   this->matrixCache_.set(path, Clo, true);
    // } else {
    Clo.save(path);
    // }
}

template <class MatrixType>
MatrixType DfObject::getCloMatrix(RUN_TYPE runType, int itr) {
    MatrixType Clo;
    const std::string path = this->getCloMatrixPath(runType, itr);
    // Clo = this->matrixCache_.get<MatrixType>(path);
    Clo.load(path);

    return Clo;
}

template <typename Vector>
Vector DfObject::getOccVtr(const RUN_TYPE runType) {
    const std::string fileName = this->getOccupationPath(runType);

    Vector occ;
    occ.load(fileName);
    assert(occ.getSize() == this->m_nNumOfMOs);

    return occ;
}

#endif  // DFOBJECT_H
