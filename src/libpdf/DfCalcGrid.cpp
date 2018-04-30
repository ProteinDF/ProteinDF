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

#include "DfCalcGrid.h"
#include <cmath>
#include "Fl_Geometry.h"

#include "CnError.h"
#include "TlMath.h"
#include "TlPosition.h"
#include "TlSymmetricMatrix.h"
#include "TlTime.h"
#include "TlUtils.h"

#define SQ2 1.414213562373095049
#define SQ1_2 0.707106781186547524
#define SQ1_3 0.577350269189625765
#define SMALL 1
#define TOOBIG 30.0
#define EPS 0.000000001 /* original #define EPS         0.000000000001 */

#define F13 0.3333333333333333
#define F23 0.6666666666666667
#define F43 1.3333333333333333
#define R3_2 1.2599210498948732

DfCalcGrid::DfCalcGrid(TlSerializeData* pPdfParam, int num_iter)
    : DfObject(pPdfParam),
      gridDataFilePath_("fl_Work/grids.dat"),
      m_tlOrbInfo((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]),
      m_tlOrbInfoAuxCD_((*pPdfParam)["coordinates"],
                        (*pPdfParam)["basis_set_j"]),
      m_tlOrbInfoXC_((*pPdfParam)["coordinates"],
                     (*pPdfParam)["basis_set_xc"]) {
  const TlSerializeData& pdfParam = *pPdfParam;

  this->alphaval = pdfParam["xc-potential/gxalpha/alpha-value"].getDouble();
  if (this->alphaval <= 0.0) {
    CnErr.abort("DfGrid", "", "constructor2",
                "inputted alpha value is illegal");
  }

  this->vectorelement = numOfAuxXC_;

  // const int lessauxnumber = numOfAuxXC_ < MaxAuxNum ? numOfAuxXC_ :
  // MaxAuxNum; this->workelement   = lessauxnumber * ( lessauxnumber + 1 ) / 2;
  this->nlsd_type = -1;

  // 1 = calc XC terms used Rho~(fitting Rho),
  // 0 = calc XC terms used Rho (true Rho)
  this->tilude = 0;

  this->m_sXCFunctional = TlUtils::toLower(this->m_sXCFunctional);
  if (m_sXCFunctional == "gxalpha~") {
    // Grid-type Xalpha
    // xc = -1, Nyu
    // xc =  0, Myu
    this->tilude = 1;
    this->xc = 0;
  } else if (m_sXCFunctional == "vbh~") {
    // von Barth & Hedin
    // xc =  2, Myu
    this->tilude = 1;
    this->xc = 2;
  } else if (m_sXCFunctional == "jmw~") {
    // Janak, Moruzzi & Williams
    // xc =  4, Myu
    this->tilude = 1;
    this->xc = 4;
  } else if (m_sXCFunctional == "gl~") {
    // Gunnarsson & Lundqvist
    // xc =  6, Myu
    this->tilude = 1;
    this->xc = 6;
  } else if ((m_sXCFunctional == "vwn~") || (m_sXCFunctional == "svwn~")) {
    // Vosko, Wilk & Nusair
    // xc =  8, Myu
    this->tilude = 1;
    this->xc = 8;
  } else if (m_sXCFunctional == "pz~") {
    // Perdew & Zunger
    // xc = 10, Myu
    this->tilude = 1;
    this->xc = 10;
  } else if (m_sXCFunctional == "b88vwn~") {
    // Vosko, Wilk & Nusair
    // NLSD method: X = Becke(B88), C = nill
    // nlsd_type = 0, xc = 8, Myu
    this->tilude = 1;
    this->nlsd_type = 0;
    this->xc = 8;
  } else if (m_sXCFunctional == "b88lyp~") {
    // Lee, Yang & Parr
    // NLSD method: X = Becke(B88), C = Lee, Yang & Parr(LYP)
    // nlsd_type = 0, xc = 12, Myu
    this->tilude = 1;
    this->nlsd_type = 0;
    this->xc = 12;
  } else if (m_sXCFunctional == "g96vwn~") {
    // Vosko, Wilk & Nusair
    // NLSD method: X = Gill(G96), C = nill
    // nlsd_type = 2, xc = 8, Myu
    this->tilude = 1;
    this->nlsd_type = 2;
    this->xc = 8;
  } else if (m_sXCFunctional == "g96lyp~") {
    // Lee, Yang & Parr
    // NLSD method: X = Gill(G96), C = Lee, Yang & Parr(LYP)
    // nlsd_type = 2, xc = 12, Myu
    this->tilude = 1;
    this->nlsd_type = 2;
    this->xc = 12;
  } else if (m_sXCFunctional == "gxalpha") {
    // Grid-type Xalpha
    // xc = -1, Nyu
    // xc =  0, Myu
    this->xc = 0;
  } else if (m_sXCFunctional == "vbh") {
    // von Barth & Hedin
    // xc =  2, Myu
    this->xc = 2;
  } else if (m_sXCFunctional == "jmw") {
    // Janak, Moruzzi & Williams
    // xc =  4, Myu
    this->xc = 4;
  } else if (m_sXCFunctional == "gl") {
    // Gunnarsson & Lundqvist
    // xc =  6, Myu
    this->xc = 6;
  } else if ((m_sXCFunctional == "vwn") || (m_sXCFunctional == "svwn")) {
    // Vosko, Wilk & Nusair
    // xc =  8, Myu
    this->xc = 8;
  } else if (m_sXCFunctional == "pz") {
    // Perdew & Zunger
    // xc = 10, Myu
    this->xc = 10;
  } else if (m_sXCFunctional == "b88vwn") {
    // Vosko, Wilk & Nusair
    // NLSD method: X = Becke(B88), C = nill
    // nlsd_type = 0, xc = 8, Myu
    this->nlsd_type = 0;
    this->xc = 8;
  } else if (m_sXCFunctional == "b88lyp") {
    // Lee, Yang & Parr
    // NLSD method: X = Becke(B88), C = Lee, Yang & Parr(LYP)
    // nlsd_type = 0, xc = 12, Myu
    this->nlsd_type = 0;
    this->xc = 12;
  } else if (m_sXCFunctional == "g96vwn~") {
    // Vosko, Wilk & Nusair
    // NLSD method: X = Gill(G96), C = nill
    // nlsd_type = 2, xc = 8, Myu
    this->nlsd_type = 2;
    this->xc = 8;
  } else if (m_sXCFunctional == "g96lyp") {
    // Lee, Yang & Parr
    // NLSD method: X = Gill(G96), C = Lee, Yang & Parr(LYP)
    // nlsd_type = 2, xc = 12, Myu
    this->nlsd_type = 2;
    this->xc = 12;
  } else if (m_sXCFunctional == "hf") {
    // HF
    this->log_.info("Selection of XC Type is HF. do nothing.");
    this->nlsd_type = 0;
    this->xc = 99;  // HF
  } else {
    this->log_.error("Selection of XC Type is Wrong.");
    this->log_.error("You type in [" + std::string(m_sXCFunctional) + "].");
    CnErr.abort(
        TlUtils::format("unknown XC functional(%s). @DfCalcGrid constructer",
                        m_sXCFunctional.c_str()));
  }
}

DfCalcGrid::~DfCalcGrid() {}

int DfCalcGrid::dfGrdMain() {
  this->log_.info("start");

  // call readGrid
  TlVector tmpVectorA, tmpVectorB, eTmpVector;
  this->calcXCInteg(tmpVectorA, tmpVectorB, eTmpVector);
  this->log_.info("calcXCInteg");

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->calcXCcoef_RKS(tmpVectorA, eTmpVector);
      break;

    case METHOD_UKS:
    case METHOD_ROKS:
      this->calcXCcoef_UKS(tmpVectorA, tmpVectorB, eTmpVector);
      break;

    default:
      CnErr.abort();
      break;
  }

  this->log_.info("end");

  return 0;
}

void DfCalcGrid::calcXCInteg(TlVector& tmpVectorA, TlVector& tmpVectorB,
                             TlVector& eTmpVector) {
  tmpVectorA = TlVector(this->numOfAuxXC_);
  tmpVectorB = TlVector(this->numOfAuxXC_);
  eTmpVector = TlVector(this->numOfAuxXC_);

  TlVector RhoAlphaA, RhoAlphaB;
  {
    switch (this->m_nMethodType) {
      case METHOD_RKS:
        RhoAlphaA.load(this->getRhoPath(RUN_RKS, this->m_nIteration));
        break;

      case METHOD_UKS:
      case METHOD_ROKS:
        RhoAlphaA.load(this->getRhoPath(RUN_UKS_ALPHA, this->m_nIteration));
        RhoAlphaB.load(this->getRhoPath(RUN_UKS_BETA, this->m_nIteration));
        break;

      default:
        CnErr.abort();
        break;
    }
  }

  TlMatrix gridMat = DfObject::getGridMatrix<TlMatrix>(this->m_nIteration - 1);
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->calcXCIntegRhoTilde_RKS(RhoAlphaA, &gridMat);
      break;

    case METHOD_UKS:
    case METHOD_ROKS:
      this->calcXCIntegRhoTilde_UKS(RhoAlphaA, RhoAlphaB, &gridMat);
      break;

    default:
      CnErr.abort();
      break;
  }
  DfObject::saveGridMatrix(this->m_nIteration, gridMat);

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->calcXCIntegMyuEpsilon_RKS(gridMat, tmpVectorA, eTmpVector);
      break;

    case METHOD_UKS:
    case METHOD_ROKS:
      this->calcXCIntegMyuEpsilon_UKS(gridMat, tmpVectorA, tmpVectorB,
                                      eTmpVector);
      break;

    default:
      CnErr.abort();
      break;
  }
}

double DfCalcGrid::getPrefactor(int nType, const TlPosition& pos) {
  double prefactor = 1.0;
  switch (nType) {
    case 0:
      // prefactor = 1.0;
      break;
    case 1:
      prefactor = pos.x();
      break;
    case 2:
      prefactor = pos.y();
      break;
    case 3:
      prefactor = pos.z();
      break;
    case 4:
      prefactor = pos.x() * pos.y();
      break;
    case 5:
      prefactor = pos.z() * pos.x();
      break;
    case 6:
      prefactor = pos.y() * pos.z();
      break;
    case 7:
      prefactor = pos.x() * pos.x() - pos.y() * pos.y();
      // prefactor = 0.5 * ((pos.x() * pos.x() - pos.y() * pos.y()));
      break;
    case 8:
      prefactor = 3.0 * pos.z() * pos.z() - pos.squareDistanceFrom();
      // prefactor = (1.0 / std::sqrt(12.0)) * (2.0*pos.z()*pos.z() -
      // pos.x()*pos.x() - pos.y()*pos.y());
      break;
    default:
      std::cout << "Basis Type is Wrong." << std::endl;
      break;
  }

  return prefactor;
}

void DfCalcGrid::getPrefactorForDerivative(int nType, double alpha,
                                           const TlPosition& pos,
                                           double* pPrefactorX,
                                           double* pPrefactorY,
                                           double* pPrefactorZ) {
  assert(pPrefactorX != NULL);
  assert(pPrefactorY != NULL);
  assert(pPrefactorZ != NULL);

  switch (nType) {
    case 0:
      *pPrefactorX = -2.0 * alpha * pos.x();
      *pPrefactorY = -2.0 * alpha * pos.y();
      *pPrefactorZ = -2.0 * alpha * pos.z();
      break;
    case 1:
      *pPrefactorX = (1.0 - 2.0 * alpha * pos.x() * pos.x());
      *pPrefactorY = -2.0 * alpha * pos.x() * pos.y();
      *pPrefactorZ = -2.0 * alpha * pos.x() * pos.z();
      break;
    case 2:
      *pPrefactorX = -2.0 * alpha * pos.y() * pos.x();
      *pPrefactorY = (1.0 - 2.0 * alpha * pos.y() * pos.y());
      *pPrefactorZ = -2.0 * alpha * pos.y() * pos.z();
      break;
    case 3:
      *pPrefactorX = -2.0 * alpha * pos.z() * pos.x();
      *pPrefactorY = -2.0 * alpha * pos.z() * pos.y();
      *pPrefactorZ = (1.0 - 2.0 * alpha * pos.z() * pos.z());
      break;
    case 4: {
      const double xyz = pos.x() * pos.y() * pos.z();
      *pPrefactorX = pos.y() * (1.0 - 2.0 * alpha * pos.x() * pos.x());
      *pPrefactorY = pos.x() * (1.0 - 2.0 * alpha * pos.y() * pos.y());
      *pPrefactorZ = -2.0 * alpha * xyz;
    } break;
    case 5: {
      const double xyz = pos.x() * pos.y() * pos.z();
      *pPrefactorX = pos.z() * (1.0 - 2.0 * alpha * pos.x() * pos.x());
      *pPrefactorY = -2.0 * alpha * xyz;
      *pPrefactorZ = pos.x() * (1.0 - 2.0 * alpha * pos.z() * pos.z());
    } break;
    case 6: {
      const double xyz = pos.x() * pos.y() * pos.z();
      *pPrefactorX = -2.0 * alpha * xyz;
      *pPrefactorY = pos.z() * (1.0 - 2.0 * alpha * pos.y() * pos.y());
      *pPrefactorZ = pos.y() * (1.0 - 2.0 * alpha * pos.z() * pos.z());
    } break;
    case 7: {
      const double x2_y2 = (pos.x() * pos.x() - pos.y() * pos.y());
      *pPrefactorX = -pos.x() * (alpha * x2_y2 - 1.0);
      *pPrefactorY = -pos.y() * (alpha * x2_y2 + 1.0);
      *pPrefactorZ = -pos.z() * alpha * x2_y2;
    } break;
    case 8: {
      const double tz2_r2 =
          2.0 * pos.z() * pos.z() - pos.x() * pos.x() - pos.y() * pos.y();
      *pPrefactorX = -pos.x() * (alpha * tz2_r2 + 1.0);
      *pPrefactorY = -pos.y() * (alpha * tz2_r2 + 1.0);
      *pPrefactorZ = -pos.z() * (alpha * tz2_r2 - 2.0);
    } break;
    default:
      std::cout << "Basis Type is Wrong." << std::endl;
      break;
  }
}

void DfCalcGrid::calcXCIntegRhoTilde_RKS(const TlVector& RhoAlphaA,
                                         TlMatrix* pGridMat) {
  const index_type numOfGrids = pGridMat->getNumOfRows();
  for (index_type i = 0; i < numOfGrids; ++i) {
    const TlPosition crdPoint(pGridMat->get(i, 0), pGridMat->get(i, 1),
                              pGridMat->get(i, 2));

    double gridRhoA = 0.0;
    for (int ialpha = 0; ialpha < this->m_nNumOfAux; ++ialpha) {
      // calc gAlpha
      const TlPosition pos =
          crdPoint - this->m_tlOrbInfoAuxCD_.getPosition(ialpha);

      const double distance2 = pos.squareDistanceFrom();
      const double shoulder =
          this->m_tlOrbInfoAuxCD_.getExponent(ialpha, 0) * distance2;

      if (shoulder <= TOOBIG) {
        double gAlpha = exp(-shoulder);

        const int basisType = this->m_tlOrbInfoAuxCD_.getBasisType(ialpha);
        const double prefactor = this->getPrefactor(basisType, pos);
        double coef = this->m_tlOrbInfoAuxCD_.getCoefficient(ialpha, 0);
        switch (basisType) {
          case 7:  // dxx-yy
            coef *= 0.5;
            break;

          case 8:  // dzz
            coef *= (SQ1_3 * 0.5);
            break;

          default:
            // nothing
            break;
        }

        gAlpha *= prefactor * coef;
        gridRhoA += RhoAlphaA[ialpha] * gAlpha;
      }
    }
    pGridMat->set(i, 5, gridRhoA);
  }
}

void DfCalcGrid::calcXCIntegRhoTilde_UKS(const TlVector& RhoAlphaA,
                                         const TlVector& RhoAlphaB,
                                         TlMatrix* pGridMat) {
  const index_type numOfGrids = pGridMat->getNumOfRows();
  for (index_type i = 0; i < numOfGrids; ++i) {
    const TlPosition crdPoint(pGridMat->get(i, 0), pGridMat->get(i, 1),
                              pGridMat->get(i, 2));

    double gridRhoA = 0.0;
    double gridRhoB = 0.0;
    for (int ialpha = 0; ialpha < this->m_nNumOfAux; ++ialpha) {
      // calc gAlpha
      const TlPosition pos =
          crdPoint - this->m_tlOrbInfoAuxCD_.getPosition(ialpha);

      const double distance2 = pos.squareDistanceFrom();
      const double shoulder =
          this->m_tlOrbInfoAuxCD_.getExponent(ialpha, 0) * distance2;

      if (shoulder <= TOOBIG) {
        double gAlpha = exp(-shoulder);
        const int basisType = this->m_tlOrbInfoAuxCD_.getBasisType(ialpha);

        const double prefactor = this->getPrefactor(basisType, pos);

        double coef = this->m_tlOrbInfoAuxCD_.getCoefficient(ialpha, 0);
        switch (basisType) {
          case 7:  // dxx-yy
            coef *= 0.5;
            break;

          case 8:  // dzz
            coef *= (SQ1_3 * 0.5);
            break;

          default:
            // nothing
            break;
        }

        gAlpha *= prefactor * coef;
        gridRhoA += RhoAlphaA[ialpha] * gAlpha;
        gridRhoB += RhoAlphaB[ialpha] * gAlpha;  // for SP and ROKS case
      }
    }

    pGridMat->set(i, 5, gridRhoA);
    pGridMat->set(i, 6, gridRhoB);
  }
}

void DfCalcGrid::calcXCIntegMyuEpsilon_RKS(const TlMatrix& gridMat,
                                           TlVector& tmpVectorA,
                                           TlVector& eTmpVector) {
  double Coef = 0.0;
  if (xc == 0) {
    // case of Xalpha Myu
    Coef = -this->alphaval * 3.0 / 2.0 * pow(2.0, F23);
  }

  // Myu & Epsilon case
  const double sCoef = pow(0.75 / M_PI, F13);
  const double XCoef = -0.75 * pow(1.5 / M_PI, F23);

  const index_type numOfGrids = gridMat.getNumOfRows();
  for (index_type i = 0; i < numOfGrids; ++i) {
    double gridRhoA = gridMat.get(i, 5);
    if (gridRhoA <= 0.0) {
      continue;
    }

    const double rs = sCoef * pow(gridRhoA, -F13);
    const double rsA = sCoef * pow(gridRhoA, F13);
    const TlPosition crdPoint(gridMat.get(i, 0), gridMat.get(i, 1),
                              gridMat.get(i, 2));
    const double weight = gridMat.get(i, 3);

    for (int idelta = 0; idelta < this->numOfAuxXC_; ++idelta) {
      // calc gDelta
      const TlPosition pos =
          crdPoint - this->m_tlOrbInfoXC_.getPosition(idelta);

      const double distance2 = pos.squareDistanceFrom();

      const double shoulder =
          this->m_tlOrbInfoXC_.getExponent(idelta, 0) * distance2;
      if (shoulder <= TOOBIG) {
        double gDelta = exp(-shoulder);

        const int basisType = this->m_tlOrbInfoXC_.getBasisType(idelta);
        double prefactor = this->getPrefactor(basisType, pos);
        switch (basisType) {
          case 7:  // dxx-yy
            prefactor *= 0.5;
            break;

          case 8:  // dzz
            prefactor *= (SQ1_3 * 0.5);
            break;

          default:
            // nothing
            break;
        }

        gDelta *= prefactor * this->m_tlOrbInfoXC_.getCoefficient(idelta, 0);
        const double gDweight = gDelta * weight;
        double EpsilonXP = 0.0;
        double EpsilonCP = 0.0;
        double RouDeltaEpsilonCP = 0.0;

        switch (xc) {
          case 0:
            // Xalpha
            tmpVectorA[idelta] += Coef * rsA * gDweight;
            break;
          case 2:
            // vBH
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0252 * HLfunc(rs / 30.0) * gDweight;
            RouDeltaEpsilonCP =
                0.0252 * rs * HLdrfunc(rs / 30.0) * gDweight / 90.0;
            break;
          case 4:
            // JMW
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0225 * HLfunc(rs / 21.0) * gDweight;
            RouDeltaEpsilonCP =
                0.0225 * rs * HLdrfunc(rs / 21.0) * gDweight / 63.0;
            break;
          case 6:
            // GL
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0333 * HLfunc(rs / 11.4) * gDweight;
            RouDeltaEpsilonCP =
                0.0333 * rs * HLdrfunc(rs / 11.4) * gDweight / 34.2;
            break;
          case 8:
            // VWN
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = 0.0310907 * VWNPfunc(rs) * gDweight;
            RouDeltaEpsilonCP = -0.0310907 * VWNPdrfunc(rs) * gDweight / 6.0;
            break;
          case 10:
            // PZ
            break;
          default:
            // do nothing
            break;
        }

        // for NSP
        eTmpVector[idelta] += EpsilonXP + EpsilonCP;
        tmpVectorA[idelta] +=
            4.0 * F13 * EpsilonXP + EpsilonCP + RouDeltaEpsilonCP;
      }
    }
  }
}

void DfCalcGrid::calcXCIntegMyuEpsilon_UKS(const TlMatrix& gridMat,
                                           TlVector& tmpVectorA,
                                           TlVector& tmpVectorB,
                                           TlVector& eTmpVector) {
  double Coef = 0.0;
  if (xc == 0) {
    // case of Xalpha Myu
    Coef = -this->alphaval * 3.0 / 2.0 * pow(2.0, F23);
    // for SP and ROKS case
    Coef *= R3_2;
  }

  // Myu & Epsilon case
  const double sCoef = pow(0.75 / M_PI, F13);
  const double XCoef = -0.75 * pow(1.5 / M_PI, F23);

  const index_type numOfGrids = gridMat.getNumOfRows();
  for (index_type i = 0; i < numOfGrids; ++i) {
    const double gridRhoA = gridMat.get(i, 5);
    const double gridRhoB = gridMat.get(i, 6);
    if ((gridRhoA <= 0.0) || (gridRhoB <= 0.0)) {
      continue;
    }

    const double rsA = sCoef * pow(gridRhoA, F13);
    const double rsB = sCoef * pow(gridRhoB, F13);
    const double TotalRou = gridRhoA + gridRhoB;
    const double rs = sCoef * pow(TotalRou, -F13);
    const double Theta = (gridRhoA - gridRhoB) / TotalRou;
    const double pTheta = 1.0 + Theta;
    const double mTheta = 1.0 - Theta;
    const double pTheta13 = pow(pTheta, F13);
    const double mTheta13 = pow(mTheta, F13);
    const double fTheta = this->polfunc(Theta);
    const double mfTheta = 1.0 - fTheta;
    const double fdTheta = poldrfunc(Theta);
    const double Csp = 1.0 + (pow(2.0, F13) - 1.0) * fTheta;

    const TlPosition crdPoint(gridMat.get(i, 0), gridMat.get(i, 1),
                              gridMat.get(i, 2));
    const double weight = gridMat.get(i, 3);

    for (int idelta = 0; idelta < this->numOfAuxXC_; ++idelta) {
      // calc gDelta
      const TlPosition pos =
          crdPoint - this->m_tlOrbInfoXC_.getPosition(idelta);

      const double distance2 = pos.squareDistanceFrom();

      const double shoulder =
          this->m_tlOrbInfoXC_.getExponent(idelta, 0) * distance2;
      if (shoulder <= TOOBIG) {
        double gDelta = exp(-shoulder);

        const int basisType = this->m_tlOrbInfoXC_.getBasisType(idelta);
        // const double prefactor =
        // this->getPrefactor(this->flVctMyuTable[idelta].basisType, pos);
        double prefactor = this->getPrefactor(basisType, pos);
        switch (basisType) {
          case 7:  // dxx-yy
            prefactor *= 0.5;
            break;

          case 8:  // dzz
            prefactor *= (SQ1_3 * 0.5);
            break;

          default:
            // nothing
            break;
        }

        // gDelta *= prefactor * this->flVctMyuTable[idelta].preExp;
        gDelta *= prefactor * this->m_tlOrbInfoXC_.getCoefficient(idelta, 0);
        const double gDweight = gDelta * weight;
        double EpsilonXP = 0.0;
        double EpsilonCP = 0.0, EpsilonCF = 0.0;
        double RouDeltaEpsilonCP = 0.0, RouDeltaEpsilonCF = 0.0;

        switch (xc) {
          case 0:
            // Xalpha
            tmpVectorA[idelta] += Coef * rsA * gDweight;
            tmpVectorB[idelta] += Coef * rsB * gDweight;
            break;
          case 2:
            // vBH
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0252 * HLfunc(rs / 30.0) * gDweight;
            RouDeltaEpsilonCP =
                0.0252 * rs * HLdrfunc(rs / 30.0) * gDweight / 90.0;
            // for SP and ROKS
            EpsilonCF = -0.0127 * HLfunc(rs / 75.0) * gDweight;
            RouDeltaEpsilonCF =
                0.0127 * rs * HLdrfunc(rs / 75.0) * gDweight / 225.0;
            break;
          case 4:
            // JMW
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0225 * HLfunc(rs / 21.0) * gDweight;
            RouDeltaEpsilonCP =
                0.0225 * rs * HLdrfunc(rs / 21.0) * gDweight / 63.0;
            // for SP and ROKS
            EpsilonCF = -0.0178 * HLfunc(rs / 26.46) * gDweight;
            RouDeltaEpsilonCF =
                0.0178 * rs * HLdrfunc(rs / 26.46) * gDweight / 79.38;
            break;
          case 6:
            // GL
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = -0.0333 * HLfunc(rs / 11.4) * gDweight;
            RouDeltaEpsilonCP =
                0.0333 * rs * HLdrfunc(rs / 11.4) * gDweight / 34.2;
            // for SP and ROKS
            EpsilonCF = -0.0203 * HLfunc(rs / 15.9) * gDweight;
            RouDeltaEpsilonCF =
                0.0203 * rs * HLdrfunc(rs / 15.9) * gDweight / 47.7;
            break;
          case 8:
            // VWN
            EpsilonXP = XCoef * gDweight / rs;
            EpsilonCP = 0.0310907 * VWNPfunc(rs) * gDweight;
            RouDeltaEpsilonCP = -0.0310907 * VWNPdrfunc(rs) * gDweight / 6.0;
            // for SP and ROKS
            EpsilonCF = 0.01554535 * VWNFfunc(rs) * gDweight;
            RouDeltaEpsilonCF = -0.01554535 * VWNFdrfunc(rs) * gDweight / 6.0;
            break;
          case 10:
            // PZ
            break;
          default:
            // do nothing
            break;
        }

        // for SP and ROKS
        eTmpVector[idelta] +=
            Csp * EpsilonXP + mfTheta * EpsilonCP + fTheta * EpsilonCF;
        tmpVectorA[idelta] += 4.0 * F13 * pTheta13 * EpsilonXP +
                              (mfTheta - fdTheta * mTheta) * EpsilonCP +
                              (fTheta + fdTheta * mTheta) * EpsilonCF +
                              mfTheta * RouDeltaEpsilonCP +
                              fTheta * RouDeltaEpsilonCF;
        tmpVectorB[idelta] += 4.0 * F13 * mTheta13 * EpsilonXP +
                              (mfTheta + fdTheta * pTheta) * EpsilonCP +
                              (fTheta - fdTheta * pTheta) * EpsilonCF +
                              mfTheta * RouDeltaEpsilonCP +
                              fTheta * RouDeltaEpsilonCF;
      }
    }
  }
}

// 1. initialize MyuGamma vector.  This is common for Epsilon and Nyu.
// 2. read [GammaDelta]-1 or [SigmaSigma']-1.
//    Support these matrixes are read at once.
//    We must support the case which cannot be read at once.
// 3. calculate product between them.

void DfCalcGrid::calcXCcoef_RKS(const TlVector& tmpVector,
                                const TlVector& eTmpVector) {
  TlSymmetricMatrix Sgdinv = DfObject::getSgdInvMatrix<TlSymmetricMatrix>();
  assert(Sgdinv.getNumOfRows() == this->numOfAuxXC_);

  {
    TlVector myuGamma = Sgdinv * tmpVector;
    myuGamma.save(this->getMyuPath(RUN_RKS, this->m_nIteration));
  }

  if (xc > 0) {
    TlVector epsGamma = Sgdinv * eTmpVector;
    if ((this->m_bMemorySave == false) && (this->m_bDiskUtilization == false)) {
      epsGamma.save("fl_Work/fl_Vct_Epsilon" +
                    TlUtils::xtos(this->m_nIteration));
    } else {
      epsGamma.save("fl_Work/fl_Vct_Epsilon");
    }
  }
}

void DfCalcGrid::calcXCcoef_UKS(const TlVector& tmpVectorA,
                                const TlVector& tmpVectorB,
                                const TlVector& eTmpVector) {
  TlSymmetricMatrix Sgdinv = DfObject::getSgdInvMatrix<TlSymmetricMatrix>();
  assert(Sgdinv.getNumOfRows() == this->numOfAuxXC_);

  {
    TlVector myuGammaA = Sgdinv * tmpVectorA;
    myuGammaA.save(this->getMyuPath(RUN_UKS_ALPHA, this->m_nIteration));
    TlVector myuGammaB = Sgdinv * tmpVectorB;
    myuGammaB.save(this->getMyuPath(RUN_UKS_BETA, this->m_nIteration));
  }

  if (xc > 0) {
    TlVector epsGamma = Sgdinv * eTmpVector;
    if ((this->m_bMemorySave == false) && (this->m_bDiskUtilization == false)) {
      epsGamma.save("fl_Work/fl_Vct_Epsilona" +
                    TlUtils::xtos(this->m_nIteration));
      epsGamma.save("fl_Work/fl_Vct_Epsilonb" +
                    TlUtils::xtos(this->m_nIteration));
    } else {
      epsGamma.save("fl_Work/fl_Vct_Epsilona");
      epsGamma.save("fl_Work/fl_Vct_Epsilonb");
    }
  }
}

// calculation polarization function
double DfCalcGrid::polfunc(double z) {
  const double p = pow(1.0 + z, F43);
  const double m = pow(1.0 - z, F43);
  const double f = 0.5 * (p + m - 2.0) / (R3_2 - 1.0);

  return f;
}

// calculation drived polarization function
double DfCalcGrid::poldrfunc(double z) {
  const double p = pow(1.0 + z, F13);
  const double m = pow(1.0 - z, F13);
  const double f = 2.0 * F13 * (p - m) / (R3_2 - 1.0);

  return f;
}

// calculation Hedin & Lundqvist function
// it is also used on JMW and GL.
double DfCalcGrid::HLfunc(double z) {
  double fz = (1.0 + z * z * z) * log(1.0 + 1.0 / z);
  fz += z / 2.0 - z * z - F13;

  return fz;
}

// calculation Hedin & Lundqvist drived function
// it is also used on JMW and GL.
double DfCalcGrid::HLdrfunc(double z) {
  double dfz = 3.0 * z * z * log(1.0 + 1.0 / z);
  dfz += 1.5 - 3.0 * z - 1.0 / z;

  return dfz;
}

// calculation VWN P function
double DfCalcGrid::VWNPfunc(double z) {
  const double aP = -0.10498;
  const double bP = 3.72744;
  const double CP = 12.9352;
  const double QP = 6.15199081976;   // QP  = sqrt( 4.0*CP - bP*bP )
  const double XPa = 12.5549141492;  // XPa = aP*aP + bP*aP + CP

  const double x = sqrt(z);
  const double XPx = x * x + bP * x + CP;
  const double atQ = atan(QP / (2.0 * x + bP));

  double fz = log((x - aP) * (x - aP) / XPx);
  fz += 2.0 * (bP + 2.0 * aP) * atQ / QP;
  fz *= -bP * aP / XPa;
  fz += log(x * x / XPx) + 2.0 * bP * atQ / QP;

  return fz;
}

// calculation VWN P derived function
double DfCalcGrid::VWNPdrfunc(double z) {
  double aP = -0.10498;
  double bP = 3.72744;
  double CP = 12.9352;
  double QP = 6.15199081976;   // QP  = sqrt( 4.0*CP - bP*bP )
  double XPa = 12.5549141492;  // XPa = aP*aP + bP*aP + CP

  const double x = sqrt(z);
  const double XPx = x * x + bP * x + CP;
  const double xbP = 2.0 * x + bP;
  const double xQQ = 4.0 * x / (QP * QP + xbP * xbP);
  const double xxX = xbP * x / XPx;

  double dfz = xQQ * (bP + 2.0 * aP);
  dfz -= 2.0 * x / (x - aP) - xxX;
  dfz *= bP * aP / XPa;
  dfz += 2.0 - xxX - xQQ * bP;

  return dfz;
}

// calculation VWN F function
double DfCalcGrid::VWNFfunc(double z) {
  double aP = -0.32500;
  double bP = 7.06042;
  double CP = 18.0578;
  double QP = 4.73092690956;  // QP  = sqrt( 4.0*CP - bP*bP )
  double XPa = 15.8687885;    // XPa = aP*aP + bP*aP + CP

  const double x = sqrt(z);
  const double XPx = x * x + bP * x + CP;
  const double atQ = atan(QP / (2.0 * x + bP)) / QP;

  double fz = log((x - aP) * (x - aP) / XPx);
  fz += 2.0 * (bP + 2.0 * aP) * atQ;
  fz *= -bP * aP / XPa;
  fz += log(x * x / XPx) + 2.0 * bP * atQ;

  return fz;
}

// calculation VWN F derived function
double DfCalcGrid::VWNFdrfunc(double z) {
  double aP = -0.32500;
  double bP = 7.06042;
  double CP = 18.0578;
  double QP = 4.73092690956;  // QP  = sqrt( 4.0*CP - bP*bP )
  double XPa = 15.8687885;    // XPa = aP*aP + bP*aP + CP

  const double x = sqrt(z);
  const double XPx = x * x + bP * x + CP;
  const double xbP = 2.0 * x + bP;
  const double xQQ = 4.0 * x / (QP * QP + xbP * xbP);
  const double xxX = xbP * x / XPx;

  double dfz = xQQ * (bP + 2.0 * aP);
  dfz -= 2.0 * x / (x - aP) - xxX;
  dfz *= bP * aP / XPa;
  dfz += 2.0 - xxX - xQQ * bP;

  return dfz;
}

// calculate b88 functional
double DfCalcGrid::B88func(double RouA, double RouB, double XA, double XB) {
  double GCex, GCexA, GCexB;
  double b = 0.0042;
  double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

  GCexA = -b * XA * XA / (1.0 + 6.0 * b * XA * asinh(XA));
  GCexB = -b * XB * XB / (1.0 + 6.0 * b * XB * asinh(XB));
  GCex = pow(RouA, F43) * (GCexA - c) + pow(RouB, F43) * (GCexB - c);

  return GCex;
}

// calculate b88 functional except const
double DfCalcGrid::DB88func(double RouA, double RouB, double XA, double XB) {
  double GCex, GCexA, GCexB;
  double b = 0.0042;

  GCexA = -b * XA * XA / (1.0 + 6.0 * b * XA * asinh(XA));
  GCexB = -b * XB * XB / (1.0 + 6.0 * b * XB * asinh(XB));
  GCex = pow(RouA, F43) * GCexA + pow(RouB, F43) * GCexB;

  return GCex;
}

double DfCalcGrid::B88dfunc(double Rou, double X, double G, double gRx,
                            double gRy, double gRz, double g, double ggx,
                            double ggy, double ggz) {
  // calculate first derivative of b88 functional
  double dGCex;
  double gX, dgX;
  double dB_dR, dB_dG;
  double XX, ish, gRgg, osbXi;
  double b = 0.0042;
  double bb = 0.00001764;         // b * b
  double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

  XX = X * X;
  ish = asinh(X);
  gRgg = gRx * ggx + gRy * ggy + gRz * ggz;
  osbXi = 1.0 + 6.0 * b * X * ish;

  gX = -b * XX / osbXi - c;
  dgX = (6.0 * bb * XX * (X / sqrt(XX + 1.0) - ish) - 2.0 * b * X) /
        (osbXi * osbXi);

  dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
  dB_dG = 0.5 * dgX / sqrt(G);

  dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

  return dGCex;
}

// calculate first derivative of b88 functional except const
double DfCalcGrid::DB88dfunc(double Rou, double X, double G, double gRx,
                             double gRy, double gRz, double g, double ggx,
                             double ggy, double ggz) {
  double dGCex;
  double gX, dgX;
  double dB_dR, dB_dG;
  double XX, ish, gRgg, osbXi;
  double b = 0.0042;
  double bb = 0.00001764;  // b * b

  XX = X * X;
  ish = asinh(X);
  gRgg = gRx * ggx + gRy * ggy + gRz * ggz;
  osbXi = 1.0 + 6.0 * b * X * ish;

  gX = -b * XX / osbXi;
  dgX = (6.0 * bb * XX * (X / sqrt(XX + 1.0) - ish) - 2.0 * b * X) /
        (osbXi * osbXi);

  dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
  dB_dG = 0.5 * dgX / sqrt(G);

  dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

  return dGCex;
}

double DfCalcGrid::G96func(double RouA, double RouB, double XA, double XB) {
  // calculate g96 functional
  double GCex, GCexA, GCexB;
  double B = 0.0072992700729927;  // 1.0/137.0;
  double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

  GCexA = -B * pow(XA, 1.5);
  GCexB = -B * pow(XB, 1.5);
  GCex = pow(RouA, F43) * (GCexA - c) + pow(RouB, F43) * (GCexB - c);

  return GCex;
}

// calculate g96 functional except const
double DfCalcGrid::DG96func(double RouA, double RouB, double XA, double XB) {
  double GCex, GCexA, GCexB;
  double B = 0.0072992700729927;  // 1.0/137.0;

  GCexA = -B * pow(XA, 1.5);
  GCexB = -B * pow(XB, 1.5);
  GCex = pow(RouA, F43) * GCexA + pow(RouB, F43) * GCexB;

  return GCex;
}

// calculate first derivative of g96 functional
double DfCalcGrid::G96dfunc(double Rou, double X, double G, double gRx,
                            double gRy, double gRz, double g, double ggx,
                            double ggy, double ggz) {
  double dGCex;
  double gX, dgX;
  double dB_dR, dB_dG;
  double gRgg;
  double B = 0.0072992700729927;  // 1.0/137.0;
  double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

  gRgg = gRx * ggx + gRy * ggy + gRz * ggz;

  gX = -B * pow(X, 1.5) - c;
  dgX = -B * 1.5 * sqrt(X);

  dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
  dB_dG = 0.5 * dgX / sqrt(G);

  dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

  return dGCex;
}

// calculate first derivative of g96 functional except
double DfCalcGrid::DG96dfunc(double Rou, double X, double G, double gRx,
                             double gRy, double gRz, double g, double ggx,
                             double ggy, double ggz) {
  double dGCex;
  double gX, dgX;
  double dB_dR, dB_dG;
  double gRgg;
  double B = 0.0072992700729927;  // 1.0/137.0;

  gRgg = gRx * ggx + gRy * ggy + gRz * ggz;

  gX = -B * pow(X, 1.5);
  dgX = -B * 1.5 * sqrt(X);

  dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
  dB_dG = 0.5 * dgX / sqrt(G);

  dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

  return dGCex;
}

// calculate lyp functional
double DfCalcGrid::LYPfunc(double TRou, double RouA, double RouB, double GAA,
                           double GAB, double GBB) {
  double a, b, c, d;
  double Cfp =
      36.4623989787647670;  // pow( 2, 11/3 ) * 3/10 * pow( 3*pai*pai, 2/3 )
  double Omega, Delta;
  double Rou13, abO, RouAB, odR;
  double F19, ttDelta, Delta11;
  double dLYP_dGAA, dLYP_dGAB, dLYP_dGBB;
  double LYP1, LYP2, LYP3;
  double LYP;

  a = 0.04918;  // Case1 of Parameter   :Default
  b = 0.132;
  c = 0.2533;
  d = 0.349;
  //    a = 0.049;                // Case2
  //    b = 0.108;
  //    c = 0.24;
  //    d = 0.342;

  Rou13 = pow(TRou, -F13);
  RouAB = RouA * RouB;
  odR = 1.0 + d * Rou13;
  Omega = exp(-c * Rou13) * pow(TRou, -11.0 / 3.0) / odR;
  Delta = c * Rou13 + d * Rou13 / odR;
  abO = a * b * Omega;

  F19 = 1.0 / 9.0;
  ttDelta = 1.0 - 3.0 * Delta;
  Delta11 = Delta - 11.0;

  dLYP_dGAA =
      -abO * (F19 * RouAB * (ttDelta - Delta11 * RouA / TRou) - RouB * RouB);
  dLYP_dGBB =
      -abO * (F19 * RouAB * (ttDelta - Delta11 * RouB / TRou) - RouA * RouA);
  dLYP_dGAB = -abO * (F19 * RouAB * (47.0 - 7.0 * Delta) - F43 * TRou * TRou);

  LYP1 = -4.0 * a * RouAB / (odR * TRou);
  LYP2 = -Cfp * abO * RouAB * (pow(RouA, 8.0 / 3.0) + pow(RouB, 8.0 / 3.0));
  LYP3 = dLYP_dGAA * GAA + dLYP_dGAB * GAB + dLYP_dGBB * GBB;
  LYP = LYP1 + LYP2 + LYP3;

  return LYP;
}

// calculate first derivative of lyp functional
double DfCalcGrid::LYPdfunc(double TRou, double RouA, double RouB, double GAA,
                            double GAB, double GBB, double gRAx, double gRAy,
                            double gRAz, double gRBx, double gRBy, double gRBz,
                            double g, double ggx, double ggy, double ggz) {
  double a, b, c, d, F19, F83;
  double Cfp =
      36.4623989787647670;  // pow( 2, 11/3 ) * 3/10 * pow( 3*pai*pai, 2/3 )
  double Omega, Delta;
  double dOmega, dDelta, abdO, dOmegap;
  double ttDelta, Delta11, fssDelta;
  double Rou13, Rou43, abO, RouAA, RouAB, RouBB, odR;
  double RouA83, RouB83, RouA_R, RouB_R;
  double gRAgg, gRBgg;
  double dLYP1, dLYP2, dLYP3;
  double dLYP_dGAA, dLYP_dGAB, dLYP_dGBB;
  double ddLYP_dRdGAA, ddLYP_dRdGAB, ddLYP_dRdGBB;
  double dLYP_dR;
  double dLYP;

  a = 0.04918;  // Case1 of Parameter   :Default
  b = 0.132;
  c = 0.2533;
  d = 0.349;
  //    a = 0.049;                // Case2
  //    b = 0.108;
  //    c = 0.24;
  //    d = 0.342;

  F19 = 1.0 / 9.0;
  F83 = 8.0 / 3.0;
  gRAgg = gRAx * ggx + gRAy * ggy + gRAz * ggz;
  gRBgg = gRBx * ggx + gRBy * ggy + gRBz * ggz;

  Rou13 = pow(TRou, -F13);
  Rou43 = pow(TRou, -F43);
  RouA83 = pow(RouA, F83);
  RouB83 = pow(RouB, F83);
  RouAA = RouA * RouA;
  RouBB = RouB * RouB;
  RouAB = RouA * RouB;
  RouA_R = RouA / TRou;
  RouB_R = RouB / TRou;
  odR = 1.0 + d * Rou13;

  Omega = exp(-c * Rou13) * pow(TRou, -11.0 / 3.0) / odR;
  Delta = c * Rou13 + d * Rou13 / odR;
  abO = a * b * Omega;
  ttDelta = 1.0 - 3.0 * Delta;
  Delta11 = Delta - 11.0;
  fssDelta = 47.0 - 7.0 * Delta;

  dOmegap = -F13 * Rou43 * (11.0 / Rou13 - c - d / odR);
  dOmega = dOmegap * Omega;
  dDelta = F13 * (d * d * pow(TRou, -5.0 / 3.0) / (odR * odR) - Delta / TRou);
  abdO = a * b * dOmega;

  dLYP_dGAA = -abO * (F19 * RouAB * (ttDelta - Delta11 * RouA_R) - RouBB);
  dLYP_dGBB = -abO * (F19 * RouAB * (ttDelta - Delta11 * RouB_R) - RouAA);
  dLYP_dGAB = -abO * (F19 * RouAB * fssDelta - F43 * TRou * TRou);

  ddLYP_dRdGAA =
      dOmegap * dLYP_dGAA -
      abO * (F19 * RouB * (ttDelta - Delta11 * RouA_R) -
             F19 * RouAB * ((3.0 + RouA_R) * dDelta + Delta11 * RouB_R / TRou));
  ddLYP_dRdGAB =
      dOmegap * dLYP_dGAB -
      abO * (F19 * RouB * fssDelta - 7.0 * F19 * RouAB * dDelta - F83 * TRou);
  ddLYP_dRdGBB =
      dOmegap * dLYP_dGBB -
      abO * (F19 * RouB * (ttDelta - Delta11 * RouB_R) -
             F19 * RouAB * ((3.0 + RouB_R) * dDelta - Delta11 * RouB_R / TRou) -
             2.0 * RouA);

  dLYP1 = -4.0 * a * RouAB * (F13 * d * Rou43 / odR + 1.0 / RouA - 1.0 / TRou) /
          (odR * TRou);
  dLYP2 = -Cfp * abdO * RouAB * (RouA83 + RouB83) -
          Cfp * abO * RouB * (11.0 / 3.0 * RouA83 + RouB83);
  dLYP3 = ddLYP_dRdGAA * GAA + ddLYP_dRdGAB * GAB + ddLYP_dRdGBB * GBB;
  dLYP_dR = dLYP1 + dLYP2 + dLYP3;

  dLYP = dLYP_dR * g + 2.0 * dLYP_dGAA * gRAgg + dLYP_dGAB * gRBgg;

  return dLYP;
}
