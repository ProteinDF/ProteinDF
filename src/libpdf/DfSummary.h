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

#ifndef DFSUMMARY_H
#define DFSUMMARY_H

#include <string>

#include "CnError.h"
#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "tl_dense_vector_lapack.h"

/** 計算結果の出力を行うクラス
 */
class DfSummary : public DfObject {
 public:
  DfSummary(TlSerializeData* pPdfParam);
  virtual ~DfSummary();

  void exec();

 protected:
  template <class MatrixType>
  void exec();

  /** 各軌道のエネルギーを出力する
   */
  void printEigen(DfObject::RUN_TYPE runType);

  /** LCAO展開係数を出力する
   */
  template <class MatrixType>
  void printMO(const MatrixType& C);

  /** LCAO展開係数、電子占有数データをファイルへ出力する
   */
  template <class MatrixType>
  void saveGuessFile(const DfObject::RUN_TYPE runType, const MatrixType& C,
                     const std::string& suffix);

  /** 密度行列展開係数、交換相関ポテンシャル展開係数を出力する
   */
  void printAux(const TlDenseVector_Lapack& rho,
                const TlDenseVector_Lapack& myu,
                const TlDenseVector_Lapack& nyu);

  /** rho population を出力する
   */
  void printRhoPop(const TlDenseVector_Lapack& rho);

 private:
};

template <class MatrixType>
void DfSummary::exec() {
  this->logger(TlUtils::format("SCF %d th iteration\n", this->m_nIteration));

  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      this->logger("eigenvalues of KS matrix");
      this->printEigen(DfObject::RUN_RKS);
      this->logger("\n");

      this->logger("LCAO coefficent\n");
      const MatrixType C =
          this->getCMatrix<MatrixType>(DfObject::RUN_RKS, this->m_nIteration);
      this->printMO(C);
      this->logger("\n");

      // this->logger("output guess.lcao.rks and guess.occ.rks");
      this->saveGuessFile(DfObject::RUN_RKS, C, "rks");

      const TlDenseVector_Lapack rho = this->getRho<TlDenseVector_Lapack>(
          DfObject::RUN_RKS, this->m_nIteration);

      TlDenseVector_Lapack myu;
      if (this->m_bIsXCFitting == true) {
        myu.load(this->getMyuPath(DfObject::RUN_RKS, this->m_nIteration));
      }

      TlDenseVector_Lapack nyu;
      if (this->m_sXCFunctional == "xalpha") {
        nyu.load(this->getNyuPath(DfObject::RUN_RKS, this->m_nIteration));
      }

      this->logger("auxiliary coefficient\n");
      this->printAux(rho, myu, nyu);
      this->logger("\n");

      this->loggerStartTitle("Rho Population\n");
      this->printRhoPop(rho);
      this->logger("\n");
    } break;

    case METHOD_ROKS: {
      this->logger("eigenvalues of KS matrix");
      this->printEigen(DfObject::RUN_ROKS);
      this->logger("\n");

      this->logger("LCAO coefficent\n");
      const MatrixType C =
          this->getCMatrix<MatrixType>(DfObject::RUN_ROKS, this->m_nIteration);
      this->printMO(C);
      this->logger("\n");

      // this->loggerStartTitle("output guess.lcao.roks and guess.occ.roks");
      this->saveGuessFile(DfObject::RUN_ROKS, C, "roks");

      const TlDenseVector_Lapack rho = this->getRho<TlDenseVector_Lapack>(
          DfObject::RUN_ROKS, this->m_nIteration);

      TlDenseVector_Lapack myu;
      if (this->m_bIsXCFitting == true) {
        myu.load(this->getMyuPath(DfObject::RUN_ROKS, this->m_nIteration));
      }

      TlDenseVector_Lapack nyu;
      if (this->m_sXCFunctional == "xalpha") {
        nyu.load(this->getNyuPath(DfObject::RUN_ROKS, this->m_nIteration));
      }

      this->logger("auxiliary coefficient\n");
      this->printAux(rho, myu, nyu);
      this->logger("\n");

      this->logger("Rho Population\n");
      this->printRhoPop(rho);
      this->logger("\n");
    } break;

    case METHOD_UKS: {
      this->logger("eigenvalues of KS matrix(alpha)\n");
      this->printEigen(DfObject::RUN_UKS_ALPHA);
      this->logger("\n");
      this->logger("eigenvalues of KS matrix(beta)\n");
      this->printEigen(DfObject::RUN_UKS_BETA);
      this->logger("\n");

      this->loggerStartTitle("LCAO coefficent(alpha)\n");
      const MatrixType CA = this->getCMatrix<MatrixType>(
          DfObject::RUN_UKS_ALPHA, this->m_nIteration);
      this->printMO(CA);
      this->logger("\n");

      this->loggerStartTitle("LCAO coefficent(beta)\n");
      const MatrixType CB = this->getCMatrix<MatrixType>(DfObject::RUN_UKS_BETA,
                                                         this->m_nIteration);
      this->printMO(CB);
      this->logger("\n");

      // this->loggerStartTitle("output guess.lcao.uks-alpha and
      // guess.occ.uks-alpha");
      this->saveGuessFile(DfObject::RUN_UKS_ALPHA, CA, "uks-alpha");

      // this->loggerStartTitle("output guess.lcao.uks-beta and
      // guess.occ.uks-beta");
      this->saveGuessFile(DfObject::RUN_UKS_BETA, CB, "uks-beta");

      const TlDenseVector_Lapack rhoA = this->getRho<TlDenseVector_Lapack>(
          DfObject::RUN_UKS_ALPHA, this->m_nIteration);

      TlDenseVector_Lapack myuA;
      if (this->m_bIsXCFitting == true) {
        myuA.load(
            this->getMyuPath(DfObject::RUN_UKS_ALPHA, this->m_nIteration));
      }

      TlDenseVector_Lapack nyuA;
      if (this->m_sXCFunctional == "xalpha") {
        nyuA.load(
            this->getNyuPath(DfObject::RUN_UKS_ALPHA, this->m_nIteration));
      }

      const TlDenseVector_Lapack rhoB = this->getRho<TlDenseVector_Lapack>(
          DfObject::RUN_UKS_BETA, this->m_nIteration);

      TlDenseVector_Lapack myuB;
      if (this->m_bIsXCFitting == true) {
        myuB.load(this->getMyuPath(DfObject::RUN_UKS_BETA, this->m_nIteration));
      }

      TlDenseVector_Lapack nyuB;
      if (this->m_sXCFunctional == "xalpha") {
        nyuB.load(this->getNyuPath(DfObject::RUN_UKS_BETA, this->m_nIteration));
      }

      this->logger("auxiliary coefficient (alpha)\n");
      this->printAux(rhoA, myuA, nyuA);
      this->logger("\n");
      this->logger("auxiliary coefficient (beta)\n");
      this->printAux(rhoB, myuB, nyuB);
      this->logger("\n");

      this->loggerStartTitle("Rho Population (alpha)");
      this->printRhoPop(rhoA);
      this->logger("\n");
      this->loggerStartTitle("Rho Population (beta)");
      this->printRhoPop(rhoB);
      this->logger("\n");
    } break;

    default:
      CnErr.abort("program error. (DfSummary)\n");
      break;
  }
}

template <class MatrixType>
void DfSummary::printMO(const MatrixType& C) {
  const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                              (*this->pPdfParam_)["basis_set"]);

  // print out LCAO coefficent
  const int numOfAOs = this->m_nNumOfAOs;
  const int numOfMOs = this->m_nNumOfMOs;
  for (int orderMO = 0; orderMO < numOfMOs; orderMO += 10) {
    this->logger("                         ");
    for (int j = orderMO; ((j < orderMO + 10) && (j < numOfMOs)); ++j) {
      this->logger(TlUtils::format("   %5ld th", j + 1));
    }
    this->logger("\n ---------------------");

    for (int j = orderMO; ((j < orderMO + 10) && (j < numOfMOs)); ++j) {
      this->logger("-----------");
    }
    this->logger("----\n");

    for (int i = 0; i < numOfAOs; ++i) {
      this->logger(TlUtils::format(
          " %6d %-2s %6s", i + 1, orbInfo.getAtomName(i).c_str(),
          TlOrbitalInfoObject::basisTypeNameTbl_[orbInfo.getBasisType(i)]));

      for (int j = orderMO; (j < orderMO + 10) && (j < numOfMOs); ++j) {
        this->logger(TlUtils::format(" %10.6lf", C.get(i, j)));
      }
      this->logger("\n");
    }
    this->logger("\n\n");
  }
  this->logger("\n\n");
}

template <class MatrixType>
void DfSummary::saveGuessFile(const DfObject::RUN_TYPE runType,
                              const MatrixType& C, const std::string& suffix) {
  C.saveText(std::string("result.guess.lcao." + suffix));

  // occupation file
  TlDenseVector_Lapack occ;
  occ.load(this->getOccupationPath(runType));

  std::ofstream ofs;
  ofs.open(std::string("result.guess.occ." + suffix).c_str(), std::ios::out);
  ofs << occ;
  ofs.close();
}

#endif  // DFSUMMARY_H
