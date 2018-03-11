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

#ifndef DFCONVCHECK_H
#define DFCONVCHECK_H

#include <string>
#include "CnError.h"
#include "DfObject.h"

// used for DfConvcheck class
struct THRESHOLD {
  double den;
  double ene;
};

/// 収束判定を行うクラス
class DfConvcheck : public DfObject {
 public:
  DfConvcheck(TlSerializeData* pPdfParam, int num_iter);
  virtual ~DfConvcheck();

  virtual bool isConverged();

 protected:
  enum TargetMatrix { CONV_TARGET_DENSITY, CONV_TARGET_FOCK };

 protected:
  virtual void check();
  virtual void showResults();

 protected:
  /// convergence check
  template <class SymmetricMatrixType>
  void check(const int iteration);

  /// check Total Energy
  void checkTotalEnergy(int iteration);

  /// interface function for "check***Matrix"
  template <class SymmetricMatrixType>
  void checkMatrix(const RUN_TYPE runType, int iteration, double* pRMS,
                   double* pMax);

  /// convergence check in DensityMatrix
  template <class SymmetricMatrixType>
  void checkDensityMatrix(const RUN_TYPE runType, int iteration, double* pRMS,
                          double* pMax);

  /// convergence check in Fock
  template <class SymmetricMatrixType>
  void checkFockMatrix(const RUN_TYPE runType, int iteration, double* pRMS,
                       double* pMax);

 protected:
  std::string getTargetMatrixStr(const RUN_TYPE runType) const;
  std::string getYN(bool yn) const;

 protected:
  bool isChecked_;
  bool isConverged_;

  TargetMatrix targetMatrix_;

  double reqRmsMatrix_;
  double rmsMatrixA_;
  double rmsMatrixB_;
  bool judgeRmsMatrixA_;
  bool judgeRmsMatrixB_;

  double reqMaxMatrix_;
  double maxMatrixA_;
  double maxMatrixB_;
  bool judgeMaxMatrixA_;
  bool judgeMaxMatrixB_;

  double deltaTotalEnergy_;
  double reqTotalEnergy_;
  bool judgeTotalEnergy_;
};

template <class SymmetricMatrixType>
void DfConvcheck::check(const int iteration) {
  this->checkTotalEnergy(iteration);

  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      double rmsMatrix = 0.0;
      double maxMatrix = 0.0;
      this->checkMatrix<SymmetricMatrixType>(RUN_RKS, iteration, &rmsMatrix,
                                             &maxMatrix);

      this->rmsMatrixA_ = rmsMatrix;
      this->maxMatrixA_ = maxMatrix;
      this->judgeRmsMatrixA_ = (this->rmsMatrixA_ < this->reqRmsMatrix_);
      this->judgeMaxMatrixA_ = (this->maxMatrixA_ < this->reqMaxMatrix_);
    }
      {
        this->judgeRmsMatrixB_ = true;
        this->judgeMaxMatrixB_ = true;
      }
      break;

    case METHOD_UKS: {
      double rmsMatrix = 0.0;
      double maxMatrix = 0.0;
      this->checkMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, iteration,
                                             &rmsMatrix, &maxMatrix);

      this->rmsMatrixA_ = rmsMatrix;
      this->maxMatrixA_ = maxMatrix;
      this->judgeRmsMatrixA_ = (this->rmsMatrixA_ < this->reqRmsMatrix_);
      this->judgeMaxMatrixA_ = (this->maxMatrixA_ < this->reqMaxMatrix_);
    }
      {
        double rmsMatrix = 0.0;
        double maxMatrix = 0.0;
        this->checkMatrix<SymmetricMatrixType>(RUN_UKS_BETA, iteration,
                                               &rmsMatrix, &maxMatrix);
        this->rmsMatrixB_ = rmsMatrix;
        this->maxMatrixB_ = maxMatrix;
        this->judgeRmsMatrixB_ = (this->rmsMatrixB_ < this->reqRmsMatrix_);
        this->judgeMaxMatrixB_ = (this->maxMatrixB_ < this->reqMaxMatrix_);
      }
      break;

    case METHOD_ROKS: {
      double rmsMatrix = 0.0;
      double maxMatrix = 0.0;
      this->checkMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSED, iteration,
                                             &rmsMatrix, &maxMatrix);
      this->rmsMatrixA_ = rmsMatrix;
      this->maxMatrixA_ = maxMatrix;
      this->judgeRmsMatrixA_ = (this->rmsMatrixA_ < this->reqRmsMatrix_);
      this->judgeMaxMatrixA_ = (this->maxMatrixA_ < this->reqMaxMatrix_);
    }
      {
        double rmsMatrix = 0.0;
        double maxMatrix = 0.0;
        this->checkMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN, iteration,
                                               &rmsMatrix, &maxMatrix);
        this->rmsMatrixB_ = rmsMatrix;
        this->maxMatrixB_ = maxMatrix;
        this->judgeRmsMatrixB_ = (this->rmsMatrixB_ < this->reqRmsMatrix_);
        this->judgeMaxMatrixB_ = (this->maxMatrixB_ < this->reqMaxMatrix_);
      }
      break;

    default:
      CnErr.abort(TlUtils::format("program error: %s.%d", __FILE__, __LINE__));
      break;
  }
}

template <class SymmetricMatrixType>
void DfConvcheck::checkMatrix(const RUN_TYPE runType, const int iteration,
                              double* pRMS, double* pMax) {
  switch (this->targetMatrix_) {
    case CONV_TARGET_DENSITY:
      this->checkDensityMatrix<SymmetricMatrixType>(runType, iteration, pRMS,
                                                    pMax);
      break;

    case CONV_TARGET_FOCK:
      this->checkFockMatrix<SymmetricMatrixType>(runType, iteration, pRMS,
                                                 pMax);
      break;

    default:
      CnErr.abort(TlUtils::format("program error: %s.%d", __FILE__, __LINE__));
      break;
  }
}

template <class SymmetricMatrixType>
void DfConvcheck::checkDensityMatrix(const RUN_TYPE runType,
                                     const int iteration, double* pRMS,
                                     double* pMax) {
  assert(iteration > 1);
  assert(pRMS != NULL);
  assert(pMax != NULL);

  SymmetricMatrixType deltaP =
      DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration);
  deltaP -= DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration - 1);

  *pRMS = deltaP.getRMS();
  index_type row, col;
  *pMax = deltaP.getMaxAbsoluteElement(&row, &col);
}

template <class SymmetricMatrixType>
void DfConvcheck::checkFockMatrix(const RUN_TYPE runType, const int iteration,
                                  double* pRMS, double* pMax) {
  assert(iteration > 1);
  assert(pRMS != NULL);
  assert(pMax != NULL);

  SymmetricMatrixType deltaF =
      DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration);
  deltaF -= DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration - 1);

  *pRMS = deltaF.getRMS();
  index_type row, col;
  *pMax = deltaF.getMaxAbsoluteElement(&row, &col);
}

#endif  // DFCONVCHECK_H
