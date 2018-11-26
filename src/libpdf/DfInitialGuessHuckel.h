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

#ifndef DFINITIALGUESSHUCKEL_H
#define DFINITIALGUESSHUCKEL_H

#include <string>

#include "DfInitialGuess.h"
#include "TlOrbitalInfo.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

/// 拡張Huckel法による初期値を作成する
class DfInitialGuessHuckel : public DfInitialGuess {
 public:
  DfInitialGuessHuckel(TlSerializeData* pPdfParam);
  virtual ~DfInitialGuessHuckel();

 public:
  virtual void createGuess();

 protected:
  void initialize();

  template <typename SymmetricMatrixType, typename MatrixType>
  void createGuess();

  template <typename SymmetricMatrixType>
  SymmetricMatrixType getHuckelMatrix();

  double getHii(const std::string& sAtomName, int nOrbitalType);

  template <typename SymmetricMatrixType, typename MatrixType>
  void generatePMatrix(const MatrixType& C, const TlDenseVector_Lapack& occ,
                       SymmetricMatrixType& P1, SymmetricMatrixType& P2);

 protected:
  static const double EPS;
  std::map<std::string, std::map<int, double> >
      m_Hii;  // [atom_number][orbital_type] = value
};

template <typename SymmetricMatrixType, typename MatrixType>
void DfInitialGuessHuckel::createGuess() {
  SymmetricMatrixType F;
  switch (initialGuessType_) {
    case GUESS_HUCKEL:
      F = this->getHuckelMatrix<SymmetricMatrixType>();
      break;
    case GUESS_CORE:
      F = DfObject::getHpqMatrix<SymmetricMatrixType>();
      break;
    default:
      this->log_.critical("unknown guess type");
      abort();
      break;
  }

  MatrixType X = DfObject::getXMatrix<MatrixType>();
  {
    MatrixType tX = X;
    tX.transposeInPlace();

    F = tX * F * X;
  }

  MatrixType C;
  {
    TlDenseVector_Lapack eigval;
    F.eig(&eigval, &C);
  }
  C = X * C;

  // P
  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      const TlDenseVector_Lapack occ = this->createOccupation(RUN_RKS);
      this->saveCMatrix(RUN_RKS, 0, C);
      this->makeDensityMatrix();
    } break;

    case METHOD_UKS: {
      const TlDenseVector_Lapack occA = this->createOccupation(RUN_UKS_ALPHA);
      const TlDenseVector_Lapack occB = this->createOccupation(RUN_UKS_BETA);
      this->saveCMatrix(RUN_UKS_ALPHA, 0, C);
      this->saveCMatrix(RUN_UKS_BETA, 0, C);
      this->makeDensityMatrix();
    } break;

    case METHOD_ROKS: {
      const TlDenseVector_Lapack occA = this->createOccupation(RUN_ROKS_CLOSED);
      const TlDenseVector_Lapack occB = this->createOccupation(RUN_ROKS_OPEN);
      this->saveCMatrix(RUN_ROKS, 0, C);
      this->makeDensityMatrix();
    } break;

    default:
      this->log_.critical("unknown method");
      abort();
      break;
  }
}

template <typename SymmetricMatrixType>
SymmetricMatrixType DfInitialGuessHuckel::getHuckelMatrix() {
  TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                        (*this->pPdfParam_)["basis_set"]);

  const int nNumOfOrbs = orbInfo.getNumOfOrbitals();

  SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
  assert(nNumOfOrbs == Hpq.getNumOfRows());

  SymmetricMatrixType Spq = DfObject::getSpqMatrix<SymmetricMatrixType>();
  assert(nNumOfOrbs == Spq.getNumOfRows());

  const double K = 1.75;

  // create Fock matrix using modified extended Huckel
  SymmetricMatrixType F(nNumOfOrbs);
  int nPrevAtomI = -1;
  int nPrevShellType = -1;
  for (int i = 0; i < nNumOfOrbs; i++) {
    const int nAtomI = orbInfo.getAtomIndex(i);
    const int nShellTypeI = orbInfo.getShellType(i);
    const double Hii = this->getHii(orbInfo.getAtomName(i), nShellTypeI);

    // i != j
    for (int j = 0; j < i; j++) {
      const int nAtomJ = orbInfo.getAtomIndex(j);
      if (nAtomI != nAtomJ) {
        const double Hjj =
            this->getHii(orbInfo.getAtomName(j), orbInfo.getShellType(j));
        F.set(i, j, 0.5 * K * (Hii + Hjj) * Spq.get(i, j));
      }
    }

    // i == j
    if ((nPrevAtomI != nAtomI) || (nPrevShellType != nShellTypeI)) {
      F.set(i, i, Hii);
      nPrevAtomI = nAtomI;
      nPrevShellType = nShellTypeI;
    }
  }

  // F.print(std::cout);
  return F;
}

template <typename SymmetricMatrixType, typename MatrixType>
void DfInitialGuessHuckel::generatePMatrix(const MatrixType& C,
                                           const TlDenseVector_Lapack& occ,
                                           SymmetricMatrixType& P1,
                                           SymmetricMatrixType& P2) {
  // generate density matrix
  const int nNumOfRows = C.getNumOfRows();
  const int nNumOfCols = C.getNumOfCols();

  P1 = SymmetricMatrixType(nNumOfRows);
  P2 = SymmetricMatrixType(nNumOfRows);

  for (int p = 0; p < nNumOfRows; p++) {
    for (int q = 0; q < nNumOfCols; q++) {
      double p1 = 0.0;
      double p2 = 0.0;

      for (int i = 0; i < occ.getSize(); i++) {
        const double occNum = occ.get(i);

        if (std::fabs(occNum - 2.0) < EPS) {
          p1 += C(p, i) * C(q, i);
        } else if (std::fabs(occNum - 1.0) < EPS) {
          p2 += C(p, i) * C(q, i);
        }
      }

      P1(p, q) = p1;
      P2(p, q) = p2;
    }
  }
}

#endif  // DFINITIALGUESSHUCKEL_H
