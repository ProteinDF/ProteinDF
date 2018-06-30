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

#include "DfFockMatrix.h"
#include "DfEriX.h"
#include "DfOverlapX.h"
#include "Fl_Geometry.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_vector_blas.h"

DfFockMatrix::DfFockMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
  this->isUseNewEngine_ = (*pPdfParam)["new_engine"].getBoolean();
}

DfFockMatrix::~DfFockMatrix() {}

void DfFockMatrix::DfFockMatrixMain() {
  if (this->m_bDiskUtilization == false) {
    // DIRECT SCHEME
    switch (this->m_nMethodType) {
      case METHOD_RKS:
        this->mainDIRECT_RKS();
        break;
      case METHOD_UKS:
        this->mainDIRECT_UKS();
        break;
      case METHOD_ROKS:
        this->mainDIRECT_ROKS();  // ROKS
        break;
      default:
        // unsupported
        CnErr.abort("unsupported SCF type. stop.");
    }
  } else {
    //     // FILE scheme
    //     switch (this->m_nMethodType) {
    //     case METHOD_RKS:
    //       this->mainFILE(RUN_RKS);  // RKS
    //       break;
    //     case METHOD_UKS:
    //       this->mainFILE(RUN_UKS_ALPHA);  // UKS alpha spin
    //       this->mainFILE(RUN_UKS_BETA);  // UKS beta spin
    //       break;
    //     case METHOD_ROKS:
    //       this->mainFILE_ROKS(); // ROKS
    //       break;
    //     default:
    //       // unsupported
    //       std::cerr << "unsupported SCF type. stop." << std::endl;
    //       CnErr.abort();
    //     }
    CnErr.abort("do coding. sorry!");
  }
}

void DfFockMatrix::mainDIRECT_RKS() {
  this->mainDIRECT_RKS<TlDenseSymmetricMatrix_BLAS_Old>();
}

void DfFockMatrix::mainDIRECT_UKS() {
  this->mainDIRECT_UKS<TlDenseSymmetricMatrix_BLAS_Old>();
}

void DfFockMatrix::mainDIRECT_ROKS() {
  this->mainDIRECT_ROKS<TlDenseGeneralMatrix_BLAS_old, TlDenseSymmetricMatrix_BLAS_Old,
                        TlVector_BLAS, DfEriX, DfOverlapX>();
}

void DfFockMatrix::setXC_RI(const RUN_TYPE nRunType,
                            TlDenseSymmetricMatrix_BLAS_Old& F) {
  this->setXC_RI<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS, DfOverlapX>(
      nRunType, F);
}

void DfFockMatrix::setXC_DIRECT(const RUN_TYPE nRunType,
                                TlDenseSymmetricMatrix_BLAS_Old& F) {
  this->setXC_DIRECT<TlDenseSymmetricMatrix_BLAS_Old>(nRunType, F);
}

void DfFockMatrix::setCoulomb(const METHOD_TYPE nMethodType,
                              TlDenseSymmetricMatrix_BLAS_Old& F) {
  TlDenseSymmetricMatrix_BLAS_Old J(this->m_nNumOfAOs);
  if (this->J_engine_ == J_ENGINE_RI_J) {
    this->setCoulomb<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS, DfEriX>(
        nMethodType, J);
    // if (this->isUseNewEngine_ == true) {
    //     this->logger(" use new engine\n");
    //     this->setCoulomb<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
    //     DfEriX>(nMethodType,
    //     J);
    // } else {
    //     this->setCoulomb<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
    //     DfEri>(nMethodType, J);
    // }
    F += J;

    // update method
    if (this->m_nIteration > 1) {
      const TlDenseSymmetricMatrix_BLAS_Old prevJ =
          DfObject::getJMatrix<TlDenseSymmetricMatrix_BLAS_Old>(this->m_nIteration -
                                                            1);
      J += prevJ;
    }

    this->saveJMatrix(this->m_nIteration, J);
  } else {
    J = this->getJMatrix<TlDenseSymmetricMatrix_BLAS_Old>(this->m_nIteration);

    // update method
    if (this->m_nIteration > 1) {
      const TlDenseSymmetricMatrix_BLAS_Old prevJ =
          DfObject::getJMatrix<TlDenseSymmetricMatrix_BLAS_Old>(this->m_nIteration -
                                                            1);
      J -= prevJ;
    }

    F += J;
  }
}

TlDenseSymmetricMatrix_BLAS_Old DfFockMatrix::getFpqMatrix(const RUN_TYPE nRunType,
                                                       const int nIteration) {
  return (DfObject::getFpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>(nRunType,
                                                              nIteration));
}

TlVector_BLAS DfFockMatrix::getRho(const RUN_TYPE nRunType,
                                   const int nIteration) {
  return (DfObject::getRho<TlVector_BLAS>(nRunType, nIteration));
}

TlVector_BLAS DfFockMatrix::getMyu(const RUN_TYPE nRunType,
                                   const int nIteration) {
  return (DfObject::getMyu<TlVector_BLAS>(nRunType, nIteration));
}

// void DfFockMatrix::saveFpqMatrix(const RUN_TYPE nRunType, const
// TlDenseSymmetricMatrix_BLAS_Old& F)
// {
//     DfObject::saveFpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>(nRunType,
//     this->m_nIteration,
//     F);
// }
