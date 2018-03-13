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

#ifndef DFTRANSFMATRIX_H
#define DFTRANSFMATRIX_H

#include <cassert>
#include <iostream>
#include <string>

#include "DfObject.h"
#include "TlTime.h"
#include "TlUtils.h"

/// Fock行列を直交基底に変換し、F'に出力する
class DfTransFmatrix : public DfObject {
 public:
  DfTransFmatrix(TlSerializeData* pPdfParam, bool bExecDiis);
  virtual ~DfTransFmatrix();

 public:
  virtual void DfTrsFmatMain();
  virtual void DfTrsFmatQclo(const std::string& fragname, int norbcut);

 protected:
  template <typename MatrixType, typename SymmetricMatrixType>
  void main(RUN_TYPE runType, const std::string& fragname = "",
            bool bPdfQcloMode = false);

 protected:
  bool m_bExecDiis;
};

// template
// (注)
// 入力はiterationの番号になる
// Fprimeの出力はtmp_iterationの番号になる
template <typename MatrixType, typename SymmetricMatrixType>
void DfTransFmatrix::main(const RUN_TYPE runType, const std::string& fragname,
                          bool bPdfQcloMode) {
  // DIIS内挿Fpqを処理するか否か
  int tmp_iteration;
  if (this->m_bExecDiis == true) {
    // DIIS extrapolated F-matrix (a temporary file)
    // DIIS内挿Fpqを処理する
    tmp_iteration = this->m_nIteration;
    this->m_nIteration = -10;  // DfConverge2.cxxを関係している
  } else {
    // ordinary F-matrix, DAMPING
    // Log << "DIIS内挿Fpqを処理しない" << "\n";
    tmp_iteration = this->m_nIteration;
    // nothing
  }

  // "read F matrix"
  SymmetricMatrixType symmF(
      DfObject::getFpqMatrix<SymmetricMatrixType>(runType, this->m_nIteration));
  MatrixType F = MatrixType(symmF);
  // symmF.save("symmF.mtx");
  // F.save("F.mtx");

  // "read X matrix"
  MatrixType X;
  if (bPdfQcloMode == false) {
    // normal DF
    X = DfObject::getXMatrix<MatrixType>();
  } else {
    // PDF-QCLO method
    X = DfObject::getCMatrix<MatrixType>(runType, tmp_iteration - 1, fragname);
  }

  // caleculate F * X
  F *= X;
  // F.save("FX.mtx");

  // make X^dagger
  X.transposeInPlace();
  // X.save("trX.mtx");

  // calculate X^dagger * (F*X)
  X *= F;
  // X.save("XFX.mtx");

  // full matrixを対称行列に変換
  SymmetricMatrixType Fprime(X);

  // save
  if (bPdfQcloMode != true) {
    DfObject::saveFprimeMatrix(runType, tmp_iteration, Fprime);
  } else {
    std::string fname = "fl_Work/fl_Mtr_Fprime.matrix.";
    fname += fragname + ".";
    fname += (this->m_sRunTypeSuffix[runType] + TlUtils::xtos(tmp_iteration));
    Fprime.save(fname);
  }
}

#endif  // DFTRANSFMATRIX_H
