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

#include "DfConvcheck.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "CnError.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

DfConvcheck::DfConvcheck(TlSerializeData* pPdfParam, int num_iter)
    : DfObject(pPdfParam) {
  const TlSerializeData& pdfParam = *pPdfParam;

  const std::string convergenceType =
      TlUtils::toUpper(pdfParam["convergence/type"].getStr());
  this->targetMatrix_ = CONV_TARGET_DENSITY;
  if (convergenceType == "FOCK") {
    this->targetMatrix_ = CONV_TARGET_FOCK;
  }

  this->reqRmsMatrix_ = pdfParam["convergence/threshold/rms"].getDouble();
  this->reqMaxMatrix_ = pdfParam["convergence/threshold"].getDouble();
  this->reqTotalEnergy_ = pdfParam["convergence/threshold_energy"].getDouble();

  this->isChecked_ = false;
}

DfConvcheck::~DfConvcheck() {}

void DfConvcheck::check() {
  // no judgement for the first iteration
  if (this->m_nIteration == 1) {
    this->isConverged_ = false;
    return;
  }

  this->log_.info("convgergence check (serial) using LAPACK");
  this->check<TlDenseSymmetricMatrix_BLAS_Old>(this->m_nIteration);

  // check for convergence
  this->isConverged_ = ((this->judgeRmsMatrixA_) && (this->judgeRmsMatrixB_) &&
                        (this->judgeMaxMatrixA_) && (this->judgeMaxMatrixB_) &&
                        (this->judgeTotalEnergy_));

  this->showResults();
}

bool DfConvcheck::isConverged() {
  if (this->isChecked_ == false) {
    this->check();
    this->isChecked_ = true;
  }

  return this->isConverged_;
}

void DfConvcheck::showResults() {
  // header
  {
    { this->log_.info(TlUtils::format("Convergence Information")); }
  }

  // matrix info
  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      this->log_.info(TlUtils::format(
          "  RMS(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_RKS).c_str(), this->rmsMatrixA_,
          this->reqRmsMatrix_, this->getYN(this->judgeRmsMatrixA_).c_str()));
      this->log_.info(TlUtils::format(
          "  MAX(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_RKS).c_str(), this->maxMatrixA_,
          this->reqMaxMatrix_, this->getYN(this->judgeMaxMatrixA_).c_str()));
    } break;

    case METHOD_UKS: {
      this->log_.info(TlUtils::format(
          "  RMS(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_UKS_ALPHA).c_str(), this->rmsMatrixA_,
          this->reqRmsMatrix_, this->getYN(this->judgeRmsMatrixA_).c_str()));
      this->log_.info(TlUtils::format(
          "  MAX(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_UKS_ALPHA).c_str(), this->maxMatrixA_,
          this->reqMaxMatrix_, this->getYN(this->judgeMaxMatrixA_).c_str()));
    }
      {
        this->log_.info(TlUtils::format(
            "  RMS(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
            this->getTargetMatrixStr(RUN_UKS_BETA).c_str(), this->rmsMatrixB_,
            this->reqRmsMatrix_, this->getYN(this->judgeRmsMatrixB_).c_str()));
        this->log_.info(TlUtils::format(
            "  MAX(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
            this->getTargetMatrixStr(RUN_UKS_BETA).c_str(), this->maxMatrixB_,
            this->reqMaxMatrix_, this->getYN(this->judgeMaxMatrixB_).c_str()));
      }
      break;

    case METHOD_ROKS: {
      this->log_.info(TlUtils::format(
          "  RMS(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_ROKS_CLOSED).c_str(), this->rmsMatrixA_,
          this->reqRmsMatrix_, this->getYN(this->judgeRmsMatrixA_).c_str()));
      this->log_.info(TlUtils::format(
          "  MAX(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
          this->getTargetMatrixStr(RUN_ROKS_CLOSED).c_str(), this->maxMatrixA_,
          this->reqMaxMatrix_, this->getYN(this->judgeMaxMatrixA_).c_str()));
    }
      {
        this->log_.info(TlUtils::format(
            "  RMS(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
            this->getTargetMatrixStr(RUN_ROKS_OPEN).c_str(), this->rmsMatrixB_,
            this->reqRmsMatrix_, this->getYN(this->judgeRmsMatrixB_).c_str()));
        this->log_.info(TlUtils::format(
            "  MAX(matrix: %7s) = %14.4le (req. %14.4le) [%3s]",
            this->getTargetMatrixStr(RUN_ROKS_OPEN).c_str(), this->maxMatrixB_,
            this->reqMaxMatrix_, this->getYN(this->judgeMaxMatrixB_).c_str()));
      }
      break;

    default:
      CnErr.abort();
      break;
  }

  // for Total Energy
  {
    {
      this->log_.info(TlUtils::format(
          "  T.E.                 = %14.4le (req. %14.4le) [%3s]",
          this->deltaTotalEnergy_, this->reqTotalEnergy_,
          this->getYN(this->judgeTotalEnergy_).c_str()));
    }
  }
}

std::string DfConvcheck::getTargetMatrixStr(const RUN_TYPE runType) const {
  std::string mat;
  switch (this->targetMatrix_) {
    case CONV_TARGET_DENSITY:
      mat = "DENS";
      break;

    case CONV_TARGET_FOCK:
      mat = "FOCK";
      break;

    default:
      CnErr.abort(TlUtils::format("program error: %s.%d", __FILE__, __LINE__));
      break;
  }

  std::string suffix;
  switch (runType) {
    case RUN_RKS:
      suffix = "   ";
      break;

    case RUN_UKS_ALPHA:
      suffix = "(a)";
      break;

    case RUN_UKS_BETA:
      suffix = "(b)";
      break;

    case RUN_ROKS_CLOSED:
      suffix = "(c)";
      break;

    case RUN_ROKS_OPEN:
      suffix = "(o)";
      break;

    default:
      CnErr.abort(TlUtils::format("program error: %s.%d", __FILE__, __LINE__));
      break;
  }

  return mat + suffix;
}

std::string DfConvcheck::getYN(bool yn) const { return (yn) ? "YES" : "NO "; }

// total energy convergence
void DfConvcheck::checkTotalEnergy(const int iteration) {
  assert(iteration > 1);

  const double TE = (*this->pPdfParam_)["TEs"][iteration].getDouble();
  const double prevTE = (*this->pPdfParam_)["TEs"][iteration - 1].getDouble();

  this->deltaTotalEnergy_ = std::fabs(TE - prevTE);
  this->judgeTotalEnergy_ = (this->deltaTotalEnergy_ < this->reqTotalEnergy_);
}
