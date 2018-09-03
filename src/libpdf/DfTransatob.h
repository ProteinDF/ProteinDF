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

#ifndef DFTRANSATOB_H
#define DFTRANSATOB_H

#include <cassert>
#include <string>

#include "CnError.h"
#include "DfObject.h"
#include "TlUtils.h"

// function : "C = X * C'"  (with Bender scheme)
//            RKS and UKS are suported
//
// input    : "X"
//            "C'"
//
// output   : "C"
//

/** X
 * 行列とC’行列を用いて、規格直交基底からCGTO基底に変換する変換行列を求め、ファイルに出力するクラス
 */
class DfTransatob : public DfObject {
 public:
  DfTransatob(TlSerializeData* pPdfParam);
  virtual ~DfTransatob();

  virtual void run();
  virtual void runQclo(const std::string& fragname, int norbcut);

 protected:
  template <typename GeneralMatrix>
  void run_method(const TlMatrixObject::index_type numOfMOs,
                  const std::string& fragname = "");

 protected:
  template <typename GeneralMatrix>
  void Cprime2C(const RUN_TYPE runtype,
                const TlMatrixObject::index_type numOfMOs,
                const std::string& fragname = "");
};

// template
template <typename GeneralMatrix>
void DfTransatob::Cprime2C(const DfObject::RUN_TYPE runType,
                           const TlMatrixObject::index_type numOfMOs,
                           const std::string& fragname) {
  // "read X matrix"
  GeneralMatrix X;
  X = this->getXMatrix<GeneralMatrix>();

  if (X.getNumOfRows() != this->m_nNumOfAOs || X.getNumOfCols() != numOfMOs) {
    this->log_.info(
        TlUtils::format("rowDim of X matrix = %d", X.getNumOfRows()));
    this->log_.info(
        TlUtils::format("colDim of X matrix = %d", X.getNumOfCols()));
    this->log_.info(TlUtils::format("number_ao_basis = %d", this->m_nNumOfAOs));
    this->log_.info(TlUtils::format("number_independant_basis = %d", numOfMOs));
    this->log_.info("DfTransatob dimension is not consistency, but continue\n");
  }

  // "read C' matrix"
  GeneralMatrix Cprime;
  Cprime = this->getCprimeMatrix<GeneralMatrix>(runType, this->m_nIteration,
                                                fragname);

  if ((Cprime.getNumOfRows() != numOfMOs) ||
      (Cprime.getNumOfCols() != numOfMOs)) {
    this->log_.info(TlUtils::format("row of C' = %d", Cprime.getNumOfRows()));
    this->log_.info(TlUtils::format("col of C' = %d", Cprime.getNumOfRows()));
    this->log_.info(TlUtils::format("number_mo_basis = %d", numOfMOs));
    this->log_.info("DfTransatob dimension is not consistency, but continue");
  }

  // calculate "C = X * C'"
  const GeneralMatrix C = X * Cprime;
  this->saveCMatrix(runType, this->m_nIteration, C);
}

template <typename GeneralMatrix>
void DfTransatob::run_method(const TlMatrixObject::index_type numOfMOs,
                             const std::string& fragment) {
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->Cprime2C<GeneralMatrix>(RUN_RKS, numOfMOs, fragment);  // RKS
      break;

    case METHOD_UKS:
      this->Cprime2C<GeneralMatrix>(RUN_UKS_ALPHA, numOfMOs, fragment);  // UKS alpha spin
      this->Cprime2C<GeneralMatrix>(RUN_UKS_BETA, numOfMOs, fragment);  // UKS beta spin
      break;

    case METHOD_ROKS:
      this->Cprime2C<GeneralMatrix>(RUN_ROKS, numOfMOs, fragment);
      break;

    default:
      CnErr.abort();
      break;
  }
}

#endif  // DFTRANSATOB_H
