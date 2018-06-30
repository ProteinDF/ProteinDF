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

#ifndef DFJMATRIX_H
#define DFJMATRIX_H

#include "DfObject.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_vector_blas.h"

class DfJMatrix : public DfObject {
 public:
  DfJMatrix(TlSerializeData* pPdfParam);
  virtual ~DfJMatrix();

 public:
  virtual void buildJ();

 protected:
  virtual void getJ_RI();
  virtual void getJ_CD();
  virtual void getJ_conventional();

  virtual void getJ_RI_local(TlDenseSymmetricMatrix_BLAS_Old* pJ);
  virtual void getJ_CD_local(TlDenseSymmetricMatrix_BLAS_Old* pJ);
  virtual void getJ_conventional_local(TlDenseSymmetricMatrix_BLAS_Old* pJ);

 protected:
  virtual void saveJMatrix(const TlDenseSymmetricMatrix_BLAS_Old& J);
  virtual TlDenseSymmetricMatrix_BLAS_Old getJMatrix(const int iteration);

  virtual TlVector_BLAS getRho(const RUN_TYPE runType, const int iteration);

  template <class SymmetricMatrixType>
  SymmetricMatrixType getDiffDensityMatrix();

  template <class SymmetricMatrixType>
  SymmetricMatrixType getDensityMatrix();
};

template <class SymmetricMatrixType>
SymmetricMatrixType DfJMatrix::getDiffDensityMatrix() {
  SymmetricMatrixType diffP;
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
          RUN_RKS, this->m_nIteration);
      break;

    case METHOD_UKS:
      diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
          RUN_UKS_ALPHA, this->m_nIteration);
      diffP += DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
          RUN_UKS_BETA, this->m_nIteration);
      break;

    case METHOD_ROKS:
      diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
          RUN_ROKS_CLOSED, this->m_nIteration);
      diffP += DfObject::getDiffDensityMatrix<SymmetricMatrixType>(
          RUN_ROKS_OPEN, this->m_nIteration);
      break;

    default:
      this->log_.critical("program error");
      break;
  }

  return diffP;
}

template <class SymmetricMatrixType>
SymmetricMatrixType DfJMatrix::getDensityMatrix() {
  SymmetricMatrixType P;
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS,
                                                      this->m_nIteration - 1);
      break;

    case METHOD_UKS:
      P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA,
                                                      this->m_nIteration - 1);
      P += DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_BETA,
                                                       this->m_nIteration - 1);
      break;

    case METHOD_ROKS:
      P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSED,
                                                      this->m_nIteration - 1);
      P += DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN,
                                                       this->m_nIteration - 1);
      break;

    default:
      this->log_.critical("program error");
      break;
  }

  return P;
}

#endif  // DFJMATRIX_H
