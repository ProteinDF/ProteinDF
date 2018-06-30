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

#ifndef DFDENSITYFITTINGOBJECT_H
#define DFDENSITYFITTINGOBJECT_H

#include <cassert>
#include <cstdlib>
#include <string>

#include "CnError.h"
#include "DfObject.h"
#include "TlFile.h"
#include "TlTime.h"
#include "TlUtils.h"

// interface
class DfDensityFittingObject : public DfObject {
 public:
  DfDensityFittingObject(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {}

  virtual ~DfDensityFittingObject() {}

 public:
  virtual void exec() = 0;
};

template <class SymmetricMatrix, class Vector, class DfERI_Class>
class DfDensityFittingTmpl : public DfDensityFittingObject {
 public:
  DfDensityFittingTmpl(TlSerializeData* pPdfParam)
      : DfDensityFittingObject(pPdfParam) {}

  virtual ~DfDensityFittingTmpl() {}

 protected:
  virtual void calc();

  virtual Vector getNalpha();
  virtual SymmetricMatrix getSabinv();
  virtual Vector calcTAlpha_DIRECT(const SymmetricMatrix& P);

  virtual void getTalpha_ROKS(Vector* pT_alphaA, Vector* pT_alphaB);
  virtual SymmetricMatrix getDiffDensityMatrix(RUN_TYPE runType);
  virtual SymmetricMatrix getP1pq(const int nIteration);
  virtual SymmetricMatrix getP2pq(const int nIteration);

  Vector getTalpha(RUN_TYPE runType);
  virtual Vector getTalpha(RUN_TYPE runType, int iteration);

  virtual double getLamda(const Vector& SabinvN, const Vector& t_alpha,
                          const Vector& N_alpha, double dNumOfElec);

  virtual void saveRho(const Vector& rRho, RUN_TYPE runType);
};

// member functions
template <class SymmetricMatrix, class Vector, class DfERI_Class>
void DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::calc() {
  const TlSerializeData& pdfParam = *(this->pPdfParam_);
  this->log_.info("start");

  const Vector N_alpha = this->getNalpha();
  const SymmetricMatrix Sabinv = this->getSabinv();
  this->log_.info(" calc S^-1*N_alpha");
  const Vector SabinvN = Sabinv * N_alpha;

  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      // RKS
      this->log_.info(" calc t_alpha");
      const Vector t_alpha = this->getTalpha(RUN_RKS);

      // calcRamda
      this->log_.info(" calc Ramda");
      const double dNumOfElec = this->m_nNumOfElectrons;
      const double lambda =
          this->getLamda(SabinvN, t_alpha, N_alpha, dNumOfElec);
      this->log_.info(
          TlUtils::format(" the number of electrons = %5.2lf", dNumOfElec));
      this->log_.info(TlUtils::format(" lambda = % 8.3f", lambda));

      // calculate Sabinv * (tAlpha - lamda * nAlpha)
      const Vector rho = Sabinv * (t_alpha - lambda * N_alpha);
      this->saveRho(rho, RUN_RKS);

      const double Sum = rho * N_alpha;
      this->log_.info(
          " following parameter must be equal to the number of electrons.");
      this->log_.info(TlUtils::format(" Sum RouAlpha*nAlpha = %f", Sum));
    } break;

    case METHOD_UKS: {
      // UKS
      this->log_.info(" calc t_alpha");
      const Vector t_alphaA = this->getTalpha(RUN_UKS_ALPHA);
      const Vector t_alphaB = this->getTalpha(RUN_UKS_BETA);

      this->log_.info(" calc Ramda");
      const double dNumOfElecA = this->m_nNumOfAlphaElectrons;
      const double dNumOfElecB = this->m_nNumOfBetaElectrons;
      const double lamdaA =
          this->getLamda(SabinvN, t_alphaA, N_alpha, dNumOfElecA);
      const double lamdaB =
          this->getLamda(SabinvN, t_alphaB, N_alpha, dNumOfElecB);
      this->log_.info(TlUtils::format(" the number of alpha elctrons = %5.2lf",
                                      dNumOfElecA));
      this->log_.info(TlUtils::format(" the number of beta  elctrons = %5.2lf",
                                      dNumOfElecB));

      // calculate Sabinv * (tAlpha - lamda * nAlpha)
      const Vector rhoA = Sabinv * (t_alphaA - lamdaA * N_alpha);
      const Vector rhoB = Sabinv * (t_alphaB - lamdaB * N_alpha);
      this->saveRho(rhoA, RUN_UKS_ALPHA);
      this->saveRho(rhoB, RUN_UKS_BETA);

      const double SumA = rhoA * N_alpha;
      const double SumB = rhoB * N_alpha;
      this->log_.info(
          " following parameter must be equal to the number of electrons.");
      this->log_.info(TlUtils::format(" Sum RouAlphaA*nAlpha = %f", SumA));
      this->log_.info(TlUtils::format(" Sum RouAlphaB*nAlpha = %f", SumB));
    } break;

    case METHOD_ROKS: {
      // ROKS
      this->log_.info(" calc t_alpha");
      Vector t_alphaA, t_alphaB;
      this->getTalpha_ROKS(&t_alphaA, &t_alphaB);

      this->log_.info(" calc Ramda");
      const double dNumOfElecA =
          pdfParam["method/roks/electron-number-alpha"].getDouble();
      const double dNumOfElecB =
          pdfParam["method/roks/electron-number-beta"].getDouble();
      const double lamdaA =
          this->getLamda(SabinvN, t_alphaA, N_alpha, dNumOfElecA);
      const double lamdaB =
          this->getLamda(SabinvN, t_alphaB, N_alpha, dNumOfElecB);
      this->log_.info(TlUtils::format(
          " the number of alpha elctrons = %5.2lf, %5.2lf", dNumOfElecA));
      this->log_.info(TlUtils::format(
          " the number of beta  elctrons = %5.2lf, %5.2lf", dNumOfElecA));

      // calculate tAlpha - lamda * nAlpha
      const Vector rhoA = Sabinv * (t_alphaA - lamdaA * N_alpha);
      const Vector rhoB = Sabinv * (t_alphaB - lamdaB * N_alpha);
      this->saveRho(rhoA, RUN_UKS_ALPHA);
      this->saveRho(rhoA, RUN_UKS_BETA);

      const double SumA = rhoA * N_alpha;
      const double SumB = rhoB * N_alpha;
      this->log_.info(
          " following parameter must be equal to the number of electrons.");
      this->log_.info(TlUtils::format("Sum RouAlphaA*nAlpha = %f", SumA));
      this->log_.info(TlUtils::format("Sum RouAlphaB*nAlpha = %f", SumB));
    } break;

    default:
      CnErr.abort();
      break;
  }

  this->log_.info("finish");
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
Vector DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getNalpha() {
  Vector N_alpha;
  N_alpha.load(DfObject::getNalphaPath());

  return N_alpha;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
SymmetricMatrix
DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getSabinv() {
  SymmetricMatrix Sabinv;
  Sabinv.load(this->getSabInvMatrixPath());

  assert(Sabinv.getNumOfRows() == this->m_nNumOfAux);

  return Sabinv;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
Vector DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getTalpha(
    const RUN_TYPE runType) {
  Vector t_alpha;

  if (this->m_bDiskUtilization == false) {
    const SymmetricMatrix diffP = this->getDiffDensityMatrix(runType);
    t_alpha = this->calcTAlpha_DIRECT(diffP);
    if (this->m_nIteration > 1) {
      t_alpha += this->getTalpha(runType, this->m_nIteration - 1);
    }

    t_alpha.save(this->getTalphaPath(runType, this->m_nIteration));
  } else {
    // To do
    abort();
  }

  return t_alpha;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
Vector
DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::calcTAlpha_DIRECT(
    const SymmetricMatrix& P) {
  Vector t_alpha(this->m_nNumOfAux);

  DfERI_Class dfEri(this->pPdfParam_);
  dfEri.getJ(P, &t_alpha);

  return t_alpha;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
Vector DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getTalpha(
    const RUN_TYPE runType, const int iteration) {
  Vector t_alpha(this->m_nNumOfAux);
  const std::string sFileName = this->getTalphaPath(runType, iteration);

  t_alpha.load(sFileName);
  return t_alpha;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
void DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getTalpha_ROKS(
    Vector* pT_alphaA, Vector* pT_alphaB) {
  assert(pT_alphaA != NULL);
  assert(pT_alphaB != NULL);

  SymmetricMatrix PA = this->getP2pq(this->m_nIteration - 1);
  SymmetricMatrix PB = this->getP1pq(this->m_nIteration - 1);
  PA += PB;

  if (this->m_bDiskUtilization == false) {
    {
      const SymmetricMatrix prevP2 = this->getP2pq(this->m_nIteration - 2);
      PA -= prevP2;
    }
    {
      const SymmetricMatrix prevP1 = this->getP1pq(this->m_nIteration - 2);
      PA -= prevP1;
      PB -= prevP1;
    }

    *pT_alphaA = this->calcTAlpha_DIRECT(PA) +
                 this->getTalpha(RUN_UKS_ALPHA, this->m_nIteration - 1);
    *pT_alphaB = this->calcTAlpha_DIRECT(PB) +
                 this->getTalpha(RUN_UKS_BETA, this->m_nIteration - 1);
  } else {
    // TODO
    abort();
  }

  // save
  if (this->m_bDiskUtilization == false) {
    pT_alphaA->save(this->getTalphaPath(RUN_UKS_ALPHA, this->m_nIteration));
    pT_alphaB->save(this->getTalphaPath(RUN_UKS_BETA, this->m_nIteration));
  }
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
SymmetricMatrix DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::
    getDiffDensityMatrix(const RUN_TYPE runType) {
  SymmetricMatrix diffP;
  diffP = DfObject::getDiffDensityMatrix<SymmetricMatrix>(runType,
                                                          this->m_nIteration);
  diffP.resize(this->m_nNumOfAOs);

  return diffP;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
SymmetricMatrix DfDensityFittingTmpl<
    SymmetricMatrix, Vector, DfERI_Class>::getP1pq(const int nIteration) {
  SymmetricMatrix P(this->m_nNumOfAOs);
  const std::string sFileName =
      TlUtils::format("fl_Work/fl_Mtr_P1pq.matrix.roks%d", nIteration);

  if (TlFile::isExistFile(sFileName) == true) {
    P.load(sFileName);
  }

  return P;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
SymmetricMatrix DfDensityFittingTmpl<
    SymmetricMatrix, Vector, DfERI_Class>::getP2pq(const int nIteration) {
  SymmetricMatrix P(this->m_nNumOfAOs);
  const std::string sFileName =
      TlUtils::format("fl_Work/fl_Mtr_P2pq.matrix.roks%d", nIteration);

  if (TlFile::isExistFile(sFileName) == true) {
    P.load(sFileName);
  }

  return P;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
double DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::getLamda(
    const Vector& SabinvN, const Vector& t_alpha, const Vector& N_alpha,
    const double dNumOfElec) {
  const double SumnSt = SabinvN * t_alpha;
  const double SumnSn = SabinvN * N_alpha;

  const double lamda = (SumnSt - dNumOfElec) / SumnSn;

  return lamda;
}

template <class SymmetricMatrix, class Vector, class DfERI_Class>
void DfDensityFittingTmpl<SymmetricMatrix, Vector, DfERI_Class>::saveRho(
    const Vector& rho, const RUN_TYPE runType) {
  DfObject::saveRho(runType, this->m_nIteration, rho);
}

#endif  // DFDENSITYFITTING_H
