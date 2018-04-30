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

#include <iostream>
#include <limits>

#include "TlCombineDensityMatrix.h"

const double TlCombineDensityMatrix::EPSILON =
    std::numeric_limits<double>::epsilon();

TlCombineDensityMatrix::TlCombineDensityMatrix(double nearThreshold,
                                               bool verbose)
    : nearThreshold_(nearThreshold),
      isCheckAtomPosition_(true),
      isVerbose_(verbose) {}

TlCombineDensityMatrix::~TlCombineDensityMatrix() {}

void TlCombineDensityMatrix::make(const TlOrbitalInfoObject& orbInfo_ref,
                                  const TlMatrixObject& P_ref,
                                  const TlOrbitalInfoObject& orbInfo_target,
                                  TlMatrixObject* pP_target) {
  const int numOfAOs_target = orbInfo_target.getNumOfOrbitals();
  if (numOfAOs_target != pP_target->getNumOfRows()) {
    std::cerr << "the number of rows in target P matrix is not consistent with "
                 "orbital information."
              << std::endl;
  }
  if (numOfAOs_target != pP_target->getNumOfCols()) {
    std::cerr << "the number of columns int target P matrix is not consistent "
                 "with orbital information."
              << std::endl;
  }

  const std::vector<int> matchTable =
      this->getMatchTable(orbInfo_ref, orbInfo_target);
  for (int row = 0; row < numOfAOs_target; ++row) {
    const int row_ref = matchTable[row];
    if (row_ref == -1) {
      continue;
    }

    for (int col = row; col < numOfAOs_target; ++col) {
      const int col_ref = matchTable[col];
      if (col_ref == -1) {
        continue;
      }

      const double value = P_ref.get(row_ref, col_ref);
      pP_target->set(row, col, value);
    }
  }
}

bool TlCombineDensityMatrix::check(const TlMatrixObject& P) {
  bool isOK = true;
  const int numOfRows = P.getNumOfRows();
  const int numOfCols = P.getNumOfCols();
  const int size = std::min(numOfRows, numOfCols);
  for (int i = 0; i < size; ++i) {
    const double value = P.get(i, i);
    if (std::fabs(value) < TlCombineDensityMatrix::EPSILON) {
      std::cerr << TlUtils::format(
                       "[WARN] %d th orbital is empty. diagonal value = %e", i,
                       value)
                << std::endl;
      isOK = false;
    }
  }

  return isOK;
}

std::vector<int> TlCombineDensityMatrix::getMatchTable(
    const TlOrbitalInfoObject& refOrbInfo,
    const TlOrbitalInfoObject& targetOrbInfo) {
  const bool verbose = this->isVerbose_;
  const double nearThreshold = this->nearThreshold_;
  const double eps = TlCombineDensityMatrix::EPSILON;
  const bool isCheckAtomPosition = this->isCheckAtomPosition_;

  const int numOfAOs_ref = refOrbInfo.getNumOfOrbitals();
  const int numOfAOs_target = targetOrbInfo.getNumOfOrbitals();
  if (verbose == true) {
    std::cerr << "number of AOs in reference = " << numOfAOs_ref << std::endl;
    std::cerr << "number of AOs in target    = " << numOfAOs_target
              << std::endl;
  }

  std::vector<int> ans(numOfAOs_target, -1);
  for (int i = 0; i < numOfAOs_target; ++i) {
    const std::string iAtomName = targetOrbInfo.getAtomName(i);
    const TlPosition iPos = targetOrbInfo.getPosition(i);
    const int iContractions = targetOrbInfo.getCgtoContraction(i);
    const int iShellType = targetOrbInfo.getShellType(i);  // s, p, d
    const int iBasisType =
        targetOrbInfo.getBasisType(i);  // s, px, py, pz, dxy, ...

    bool isMatched = false;
    for (int j = 0; j < numOfAOs_ref; ++j) {
      const std::string jAtomName = refOrbInfo.getAtomName(j);
      const TlPosition jPos = refOrbInfo.getPosition(j);
      const int jContractions = refOrbInfo.getCgtoContraction(j);
      const int jShellType = refOrbInfo.getShellType(j);
      const int jBasisType = refOrbInfo.getBasisType(j);

      // filter ========================================================
      if ((isCheckAtomPosition == true) &&
          (iPos.distanceFrom(jPos) > nearThreshold)) {
        continue;
      }

      if ((iAtomName == jAtomName) && (iContractions == jContractions) &&
          (iShellType == jShellType) && (iBasisType == jBasisType)) {
        // contract- and primitive- GTO check
        bool isPGTO_OK = true;
        for (int c = 0; c < iContractions; ++c) {
          const double targetCoef = targetOrbInfo.getCoefficient(i, c);
          const double refCoef = refOrbInfo.getCoefficient(j, c);
          const double targetExp = targetOrbInfo.getExponent(i, c);
          const double refExp = refOrbInfo.getExponent(j, c);
          if ((std::fabs(targetCoef - refCoef) > eps) ||
              (std::fabs(targetExp - refExp) > eps)) {
            isPGTO_OK = false;
            break;
          }
        }

        if (isPGTO_OK == true) {
          // PGTOs parameters are also same.
          ans[i] = j;
          isMatched = true;

          // for debug
          if (verbose == true) {
            std::cerr << TlUtils::format(
                             "[%d] atom=%2s, shell=%d, type=%d => [%d] "
                             "atom=%2s, shell=%d, type=%d",
                             i, iAtomName.c_str(), iShellType, iBasisType, j,
                             jAtomName.c_str(), jShellType, jBasisType)
                      << std::endl;
          }

          break;  // break j-loop
        }
      }
    }

    if ((verbose == true) && (isMatched == false)) {
      std::cerr << TlUtils::format(
                       "warning! [%3d] atom=%2s, cont=%d, shell=%d, type=%d: "
                       "not found.",
                       i, iAtomName.c_str(), iContractions, iShellType,
                       iBasisType)
                << std::endl;
    }
  }

  return ans;
}
