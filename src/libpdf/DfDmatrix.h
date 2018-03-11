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

#ifndef DFDMATRIX_H
#define DFDMATRIX_H

#include <string>
#include <vector>

#include "CnError.h"
#include "DfObject.h"
#include "TlSymmetricMatrix.h"
#include "TlTime.h"

class TlVector;
class TlMatrix;

/// 分子軌道の占有数の決定と密度行列の作成を行うクラス
///
/// 収束させるための処理として、軌道の重なりの方法、射影演算子法を用いる
class DfDmatrix : public DfObject {
 public:
  DfDmatrix(TlSerializeData* pPdfParam);
  virtual ~DfDmatrix();

 public:
  void DfDmatrixMain();

 protected:
  virtual void main(DfObject::RUN_TYPE runType);

  // virtual TlVector getOccupation(DfObject::RUN_TYPE runType);

  template <typename MatrixType>
  TlVector getOccupationUsingOverlap(DfObject::RUN_TYPE runType);

  template <typename MatrixType, typename SymmetricMatrixType>
  TlVector getOccupationUsingProjection(DfObject::RUN_TYPE runType);

  virtual void checkOccupation(const TlVector& prevOcc,
                               const TlVector& currOcc);
  virtual void printOccupation(const TlVector& occ);

  template <typename MatrixType, typename SymmetricMatrixType>
  void generateDensityMatrix(DfObject::RUN_TYPE runType,
                             const TlVector& currOcc);

  template <typename MatrixType, typename SymmetricMatrixType>
  SymmetricMatrixType calcDensMatrix(const MatrixType& inputC,
                                     const TlVector& f, double base);

  void printTwoVectors(const TlVector& a, const TlVector& b,
                       const std::string& title, int pnumcol);

  // for simple
  // TlVector createOccupation(DfObject::RUN_TYPE runType);
  // std::vector<int> getLevel(std::string sLevel);

 protected:
  enum ORBITAL_CORRESPONDENCE_METHOD { OCM_NONE, OCM_OVERLAP, OCM_PROJECTION };

 protected:
  ORBITAL_CORRESPONDENCE_METHOD orbitalCorrespondenceMethod_;

  //     /// 軌道関連づけを行う(true)かどうか
  //     bool isOrbitalCorrespondence_;

  //     /// 1 回目のiteration から上記の方法を用いるかのキーワード
  //     std::string orbital_overlap_first;

  //     /// 軌道重なりの方法もしくは軌道射影法の指定を行うキーワード
  //     std::string orbital_overlap_method;

  //     /// 軌道重なりの方法をiteration 何回目まで行うかを指定するキーワード
  //     int mo_overlap_iter;
};

//
// templates
//

// 軌道の重なりの対応のルーチン ====================================
template <typename MatrixType>
TlVector DfDmatrix::getOccupationUsingOverlap(DfObject::RUN_TYPE runType) {
  const int nNumOfMOs = this->m_nNumOfMOs;

  this->log_.info(" MO overlap method is started.\n");

  this->log_.info(" load previous C' matrix");
  MatrixType prevCprime(nNumOfMOs, nNumOfMOs);
  {
    // read current orbital in orthonormal basis
    MatrixType Cprime;
    Cprime = DfObject::getCprimeMatrix<MatrixType>(runType, this->m_nIteration);

    // read previous orbital in orthonormal basis
    if ((this->m_nIteration != 1) ||
        (this->initialGuessType_ == GUESS_LCAO ||
         this->initialGuessType_ == GUESS_HUCKEL)) {
      prevCprime = DfObject::getCprimeMatrix<MatrixType>(
          runType, this->m_nIteration - 1);
    } else if (this->m_nIteration == 1) {
      prevCprime = Cprime;
    }

    // calculate and print MO overlap between previous and current orbitals
    // C <-- C_pre^dagger * C_cur, in program  D <-- B^dagger * A
    prevCprime.transpose();
    prevCprime *= Cprime;
  }

  TlVector prevOcc = this->getOccVtr(runType);

  // construct occupation number of current orbital with MO overlap matrix
  // 旧MO(pre)がどの新MO(crr)との重なりが一番大きいかを探す
  this->log_.info(" check overlap");
  TlVector currOcc(nNumOfMOs);
  {
    std::vector<bool> g(nNumOfMOs, false);
    bool bListHeaderOutput = false;
    for (int pre = 0; pre < nNumOfMOs; ++pre) {
      const double dPrevOcc = prevOcc[pre];
      if ((std::fabs(dPrevOcc - 1.00) < 1.0e-10) ||  // 占有数が 1 or 2 の時
          (std::fabs(dPrevOcc - 2.00) < 1.0e-10)) {
        const TlVector prevCprimeRowVector = prevCprime.getRowVector(pre);
        double xval = 0.0;
        int xord = -1;
        for (int crr = 0; crr < nNumOfMOs; ++crr) {
          const double val = std::fabs(prevCprimeRowVector[crr]);
          if ((g[crr] == false) && (val > xval)) {
            xval = val;  // fabs(prevCprime(pre, crr));
            xord = crr;
          }
        }

        if (xord == -1) {
          this->log_.info(
              TlUtils::format(" crr MO %d th is not corresponded!\n", pre));
          // CnErr.abort("DfDmatrix", "", "", "MO Overlap is not corresponded");
        }

        if (xval < 0.3) {
          if (bListHeaderOutput == false) {
            this->log_.info(" ##### MO Overlap is less than 0.3 #####\n");
            bListHeaderOutput = true;
          }
          this->log_.info(
              TlUtils::format("pre MO %6d th -> crr MO %6d th %12.8lf\n",
                              (pre + 1), (xord + 1), xval));
        }
        currOcc[xord] = dPrevOcc;
        g[xord] = true;
      }
    }
  }

  this->log_.info(" check occupation vectors");
  this->checkOccupation(prevOcc, currOcc);

  //     if (this->m_nIteration != 1) {
  //         if (orbital_overlap != "on") {
  //             this->logger(" orbital overlap correspondence is not carried
  //             out by the inputted keyword\n"); this->logger(" temporal
  //             occupation from orbital correspondence\n");
  //             this->printOccupation(currOcc);

  //             currOcc = prevOcc;
  //         }
  //     } else {
  //         if (orbital_overlap_first != "on") {
  //             this->logger(" orbital overlap correspondence is not carried
  //             out by the inputted keyword\n"); this->logger(" temporal
  //             occupation from orbital correspondence\n");
  //             this->printOccupation(currOcc);

  //             // 下から電子を詰めるようにする
  //             this->logger(" the electrons are occupied from the 1 st
  //             orbital\n");

  //             currOcc = prevOcc;

  //             currOcc.sortByGrater();
  //         }
  //     }

  this->log_.info(" finish");
  return currOcc;
}

// 射影演算子法 ====================================================
template <typename MatrixType, typename SymmetricMatrixType>
TlVector DfDmatrix::getOccupationUsingProjection(
    const DfObject::RUN_TYPE runType) {
  this->log_.info("orbital_overlap_method is mo-projection");

  const index_type numOfMOs = this->m_nNumOfMOs;

  // read molecular orbital occupation of (n-1) SCF iteration
  // TlVector prevOcc = this->getOccupation(runType);
  // {
  //     int num_mo_closed  = 0;
  //     int num_mo_open    = 0;
  //     int num_mo_virtual = 0;
  //     for (index_type k = 0; k < numOfMOs; ++k) {
  //         if (std::fabs(prevOcc[k] -2.00) < 1.0e-10) {
  //             ++num_mo_closed;
  //         } else if (std::fabs(prevOcc[k] -1.00) < 1.0e-10) {
  //             ++num_mo_open;
  //         } else if (std::fabs(prevOcc[k]) < 1.0e-10) {
  //             ++num_mo_virtual;
  //         }
  //     }
  //     this->log_.info(TlUtils::format(" closed  orbital = %5ld",
  //     num_mo_closed)); this->log_.info(TlUtils::format(" open    orbital =
  //     %5ld", num_mo_open)); this->log_.info(TlUtils::format(" virtual orbital
  //     = %5ld", num_mo_virtual));
  // }

  // read molecular orbital ^(n)
  MatrixType C;
  C = DfObject::getCMatrix<MatrixType>(runType, this->m_nIteration);

  MatrixType Ct = C;
  Ct.transpose();

  SymmetricMatrixType S;
  S = DfObject::getSpqMatrix<SymmetricMatrixType>();

  // calculation of projection diagonal
  TlVector pd(numOfMOs);

  // read density matrix of (n-1) SCF iteration
  const SymmetricMatrixType D = DfObject::getPpqMatrix<SymmetricMatrixType>(
      runType, this->m_nIteration - 1);

  const MatrixType SDS = S * D * S;

  // diagonal
  const MatrixType judge = Ct * SDS * C;
  for (index_type i = 0; i < numOfMOs; ++i) {
    pd[i] = judge.get(i, i);
  }

  for (index_type i = 0; i < numOfMOs; ++i) {
    if (pd[i] < 0.0) {
      pd[i] *= -1.0;
    }
  }

  // set order of MO (0から始まる軌道の番号)
  TlVector pd_ord(numOfMOs);
  for (index_type k = 0; k < numOfMOs; ++k) {
    pd_ord[k] = k;
  }

  // sort pd
  for (index_type i = 0; i < numOfMOs - 1; ++i) {
    for (index_type j = i; j < numOfMOs; ++j) {
      if (pd[i] < pd[j]) {
        std::swap(pd[i], pd[j]);
        std::swap(pd_ord[i], pd_ord[j]);
      }
    }
  }

  // output
  double occupation = 0.0;
  std::string title = "";
  switch (runType) {
    case RUN_RKS:
      occupation = 2.0;
      title = "projection diagonal of MO";
      break;

    case RUN_UKS_ALPHA:
      occupation = 1.0;
      title = "projection diagonal of MO(alpha)";
      break;

    case RUN_UKS_BETA:
      occupation = 1.0;
      title = "projection diagonal of MO(beta)";
      break;

    case RUN_ROKS_CLOSED:
      occupation = 2.0;
      title = "projection diagonal of MO(closed)";
      break;

    case RUN_ROKS_OPEN:
      occupation = 1.0;
      title = "projection diagonal of MO(open)";
      break;

    default:
      this->log_.critical("program error.");
      break;
  }
  this->printTwoVectors(pd_ord, pd, title, 10);

  // store electrons
  TlVector currOcc = TlVector(numOfMOs);
  for (index_type k = 0; k < numOfMOs; ++k) {
    currOcc[static_cast<index_type>(pd_ord[k])] = occupation;
  }

  // check
  {
    TlVector prevOcc = this->getOccVtr(runType);
    const double sumOfPrevOcc = prevOcc.sum();
    const double sumOfCurrOcc = currOcc.sum();

    this->log_.info(TlUtils::format("prev occ = %10.2f", sumOfPrevOcc));
    this->log_.info(TlUtils::format("prev occ = %10.2f", sumOfCurrOcc));
    if (std::fabs(sumOfPrevOcc - sumOfCurrOcc) > 1.0E-5) {
      this->log_.info("projection occupation was failed.");
      CnErr.abort();
    }
  }

  return currOcc;
}

template <typename MatrixType, typename SymmetricMatrixType>
void DfDmatrix::generateDensityMatrix(const DfObject::RUN_TYPE runType,
                                      const TlVector& currOcc) {
  switch (runType) {
    case RUN_RKS: {
      MatrixType C =
          DfObject::getCMatrix<MatrixType>(runType, this->m_nIteration);
      SymmetricMatrixType P =
          this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc,
                                                                2.0);
      P *= 2.0;

      this->savePpqMatrix(runType, this->m_nIteration, P);
    } break;

    case RUN_UKS_ALPHA:
    case RUN_UKS_BETA: {
      MatrixType C =
          DfObject::getCMatrix<MatrixType>(runType, this->m_nIteration);
      SymmetricMatrixType P =
          this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc,
                                                                1.0);
      this->savePpqMatrix(runType, this->m_nIteration, P);
    } break;

    case RUN_ROKS_CLOSED: {
      MatrixType C =
          DfObject::getCMatrix<MatrixType>(RUN_ROKS, this->m_nIteration);
      SymmetricMatrixType P =
          this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc,
                                                                2.0);
      P *= 2.0;
      this->savePpqMatrix(runType, this->m_nIteration, P);

    } break;

    case RUN_ROKS_OPEN: {
      MatrixType C =
          DfObject::getCMatrix<MatrixType>(RUN_ROKS, this->m_nIteration);
      SymmetricMatrixType P =
          this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc,
                                                                1.0);
      this->savePpqMatrix(runType, this->m_nIteration, P);
    } break;

    default:
      std::cerr << " DfDmatrix::generateDensityMatrix() program error."
                << __FILE__ << __LINE__ << std::endl;
      break;
  }
}

template <typename MatrixType, typename SymmetricMatrixType>
SymmetricMatrixType DfDmatrix::calcDensMatrix(const MatrixType& inputC,
                                              const TlVector& f, double base) {
  const std::size_t nNumOfAOs = inputC.getNumOfRows();
  const std::size_t nNumOfMOs = inputC.getNumOfCols();

  MatrixType C = inputC;
  {
    MatrixType E(nNumOfMOs, nNumOfAOs);

    std::size_t max_k = std::min<std::size_t>(nNumOfMOs, f.getSize());
    for (std::size_t k = 0; k < max_k; ++k) {
      if (std::fabs(f[k] - base) < 1.0e-10) {
        E(k, k) = 1.0;
      }
    }
    C = C * E;
  }

  MatrixType Ct = C;
  Ct.transpose();

  // DEBUG版ではここで対称性のチェックを行う
  const SymmetricMatrixType P = C * Ct;

  return P;
}

#endif  // DFDMATRIX_H
