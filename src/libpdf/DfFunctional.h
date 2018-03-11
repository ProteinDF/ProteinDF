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

#ifndef DFFUNCTIONAL_H
#define DFFUNCTIONAL_H

#include <vector>
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

//
enum { FA_R = 0, FA_X = 1, FB_R = 2, FB_X = 3, F_DIM = 4 };

struct FunctionalSets {
  FunctionalSets(int numOfTerms, TlMatrix::index_type dim)
      : FA_termR(numOfTerms, dim),
        FA_termX(numOfTerms, dim),
        FB_termR(numOfTerms, dim),
        FB_termX(numOfTerms, dim){};

  TlMatrix FA_termR;
  TlMatrix FA_termX;
  TlMatrix FB_termR;
  TlMatrix FB_termX;
};

enum {
  RA_R = 0,  // term to depend on rho in round F / round rhoA
  RA_X = 1,
  RB_R = 2,
  RB_X = 3,
  GAA_R = 4,
  GAA_X = 5,
  GAB_R = 6,
  GAB_X = 7,
  GBB_R = 8,
  GBB_X = 9,
  D_DIM = 10  //
};

struct DerivativeFunctionalSets {
  DerivativeFunctionalSets(int numOfTerms, TlMatrix::index_type dim)
      : rFrRhoA_R(numOfTerms, dim),
        rFrRhoA_X(numOfTerms, dim),
        rFrRhoB_R(numOfTerms, dim),
        rFrRhoB_X(numOfTerms, dim),
        rFrGAA_R(numOfTerms, dim),
        rFrGAA_X(numOfTerms, dim),
        rFrGAB_R(numOfTerms, dim),
        rFrGAB_X(numOfTerms, dim),
        rFrGBB_R(numOfTerms, dim),
        rFrGBB_X(numOfTerms, dim){};

  TlMatrix rFrRhoA_R;  // roundF_roundRhoA_termR
  TlMatrix rFrRhoA_X;  // roundF_roundRhoA_termX
  TlMatrix rFrRhoB_R;  // roundF_roundRhoB_termR
  TlMatrix rFrRhoB_X;  // roundF_roundRhoB_termX
  TlMatrix rFrGAA_R;   // roundF_roundGammaAA_termR
  TlMatrix rFrGAA_X;   // roundF_roundGammaAA_termX
  TlMatrix rFrGAB_R;   // roundF_roundGammaAB_termR
  TlMatrix rFrGAB_X;   // roundF_roundGammaAB_termX
  TlMatrix rFrGBB_R;   // roundF_roundGammaBB_termR
  TlMatrix rFrGBB_X;   // roundF_roundGammaBB_termX
};

/** 局所汎関数を扱うインターフェースクラス
 *
 */
class DfFunctional_LDA {
 public:
  typedef TlMatrix::index_type index_type;

 public:
  DfFunctional_LDA();
  virtual ~DfFunctional_LDA();

 public:
  /** 交換相関関数値を返す
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dRhoB beta電子密度
   */
  virtual double getFunctional(double dRhoA, double dRhoB) = 0;

  /** 交換相関関数値を返す(RKSに特殊化)
   *
   *  @param[in] dRhoA alpha電子密度
   */
  double getFunctional(double dRhoA) {
    return this->getFunctional(dRhoA, dRhoA);
  }

  /** 交換相関関数の一次微分を返す
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dRhoB beta電子密度
   *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
   *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
   */
  virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                       double* pRoundF_roundRhoA,
                                       double* pRoundF_roundRhoB) = 0;

  /** 交換相関関数の一次微分を返す(RKSに特殊化)
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dRhoB beta電子密度
   *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
   *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
   */
  void getDerivativeFunctional(double dRhoA, double* pRoundF_roundRhoA) {
    double dDummyVariable;
    this->getDerivativeFunctional(dRhoA, dRhoA, pRoundF_roundRhoA,
                                  &dDummyVariable);
  }

  // for Grid-Free ===================================================
 public:
  FunctionalSets getFunctional_GF(const TlVector& rhoAs, const TlVector& rhoBs);

  DerivativeFunctionalSets getDerivativeFunctional_GF(const TlVector& rhoAs,
                                                      const TlVector& rhoBs);

  int getNumOfFunctionalTerms() const { return this->numOfFunctionalTerms_; }

  int getNumOfDerivativeFunctionalTerms() const {
    return this->numOfDerivativeFunctionalTerms_;
  }

 public:
  virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB) {
    return TlMatrix();
  }

  virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                               const double rhoB) {
    return TlMatrix();
  }

 protected:
  int numOfFunctionalTerms_;
  int numOfDerivativeFunctionalTerms_;
};

/** 非局所汎関数を扱うインターフェースクラス
 *
 */
class DfFunctional_GGA {
 public:
  typedef TlMatrix::index_type index_type;

 public:
  DfFunctional_GGA();
  virtual ~DfFunctional_GGA();

 public:
  /** 交換相関関数値を返す
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dRhoB beta電子密度
   *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
   *  @param[in] dGammaAB (grad Rho alpha) * (grad Rho beta)
   *  @param[in] dGammaBB (grad Rho beta) * (grad Rho beta)
   */
  virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA,
                               double dGammaAB, double dGammaBB) = 0;

  /** 交換相関関数値を返す(RKSに特殊化)
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
   */
  double getFunctional(double dRhoA, double dGammaAA) {
    return this->getFunctional(dRhoA, dRhoA, dGammaAA, dGammaAA, dGammaAA);
  }

  /** 交換相関関数の一次微分を返す
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dRhoB beta電子密度
   *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
   *  @param[in] dGammaAB (grad Rho alpha) * (grad Rho beta)
   *  @param[in] dGammaBB (grad Rho beta) * (grad Rho beta)
   *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
   *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
   *  @param[out] pRoundF_roundGammaAA (gamma alpha alpha)による一次偏微分
   *  @param[out] pRoundF_roundGammaAB (gamma alpha beta)による一次偏微分
   *  @param[out] pRoundF_roundGammaBB (gamma beta beta)による一次偏微分
   */
  virtual void getDerivativeFunctional(
      double dRhoA, double dRhoB, double dGammaAA, double dGammaAB,
      double dGammaBB, double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
      double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB,
      double* pRoundF_roundGammaBB) = 0;

  /** 交換相関関数の一次微分を返す(RKSに特殊化)
   *
   *  @param[in] dRhoA alpha電子密度
   *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
   *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
   *  @param[out] pRoundF_roundGammaAA (gamma alpha alpha)による一次偏微分
   *  @param[out] pRoundF_roundGammaAB (gamma alpha beta)による一次偏微分
   */
  void getDerivativeFunctional(double dRhoA, double dGammaAA,
                               double* pRoundF_roundRhoA,
                               double* pRoundF_roundGammaAA,
                               double* pRoundF_roundGammaAB) {
    double dRoundF_roundRhoB;
    double dRoundF_roundGammaBB;

    this->getDerivativeFunctional(dRhoA, dRhoA, dGammaAA, dGammaAA, dGammaAA,
                                  pRoundF_roundRhoA, &dRoundF_roundRhoB,
                                  pRoundF_roundGammaAA, pRoundF_roundGammaAB,
                                  &dRoundF_roundGammaBB);
  }

  // for Grid-Free ===================================================
 public:
  FunctionalSets getFunctional_GF(const TlVector& rhoAs, const TlVector& rhoBs,
                                  const TlVector& xAs, const TlVector& xBs);

  DerivativeFunctionalSets getDerivativeFunctional_GF(const TlVector& rhoAs,
                                                      const TlVector& rhoBs,
                                                      const TlVector& xAs,
                                                      const TlVector& xBs);

  int getNumOfFunctionalTerms() const { return this->numOfFunctionalTerms_; }

  int getNumOfDerivativeFunctionalTerms() const {
    return this->numOfDerivativeFunctionalTerms_;
  }

 public:
  virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB,
                                     const double xA, const double xB) {
    return TlMatrix();
  }

  virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                               const double rhoB,
                                               const double xA,
                                               const double xB) {
    return TlMatrix();
  }

 protected:
  int numOfFunctionalTerms_;
  int numOfDerivativeFunctionalTerms_;
};

#endif  // DFFUNCTIONAL_H
