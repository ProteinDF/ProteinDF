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

#ifndef DFCALCGRIDX_H
#define DFCALCGRIDX_H

#include <string>
#include <vector>

#include "DfFunctional.h"
#include "DfObject.h"
#include "DfXCFunctional.h"
#include "TlOrbitalInfo.h"

class TlMatrixObject;

class DfCalcGridX : public DfObject {
 public:
  enum {
    // COMMON
    GM_X = 0,
    GM_Y = 1,
    GM_Z = 2,
    GM_WEIGHT = 3,
    GM_ATOM_INDEX = 4,
    // for LDA
    GM_LDA_RHO_ALPHA = 5,
    GM_LDA_RHO_BETA = 6,
    // for GGA
    GM_GGA_RHO_ALPHA = 5,
    GM_GGA_GRAD_RHO_X_ALPHA = 6,
    GM_GGA_GRAD_RHO_Y_ALPHA = 7,
    GM_GGA_GRAD_RHO_Z_ALPHA = 8,
    GM_GGA_RHO_BETA = 9,
    GM_GGA_GRAD_RHO_X_BETA = 10,
    GM_GGA_GRAD_RHO_Y_BETA = 11,
    GM_GGA_GRAD_RHO_Z_BETA = 12
  };

 public:
  // the value of wave function on the grid
  struct WFGrid {
    WFGrid(std::size_t i = 0, double v = 0.0) : index(i), value(v) {}

    std::size_t index;
    double value;
  };

  struct WFGrid_sort_functional {
    bool operator()(const WFGrid& a, const WFGrid& b) {
      return (std::fabs(a.value) > std::fabs(b.value));
    }
  };

 public:
  DfCalcGridX(TlSerializeData* pPdfParam);
  virtual ~DfCalcGridX();

 public:
  // for density field data
  void gridDensity(const TlSymmetricMatrix& P, const TlPosition& gridPosition,
                   double* pRho);

 public:
  // for derivative
  virtual TlMatrix energyGradient(const TlSymmetricMatrix& P_A,
                                  DfFunctional_LDA* pFunctional);
  virtual TlMatrix energyGradient(const TlSymmetricMatrix& P_A,
                                  DfFunctional_GGA* pFunctional);

 protected:
  TlMatrix selectGridMatrixByAtom(const TlMatrix& globalGridMat,
                                  const int atomIndex);

 protected:
  virtual void defineCutOffValues(const TlSymmetricMatrix& P);

  virtual void defineCutOffValues(const TlSymmetricMatrix& PA,
                                  const TlSymmetricMatrix& PB);

  /// グリッド点における波動関数値配列(std::vector<WFGrid>)1つから
  /// 交換相関項(Fxc)を求める。
  ///
  /// @param [out] pF
  /// 交換相関項行列。この行列は実対称行列オブジェクトを指定すること。
  void buildFock(std::vector<WFGrid>::const_iterator pBegin,
                 std::vector<WFGrid>::const_iterator pEnd, double coef,
                 double cutoffValue, TlMatrixObject* pFxc);

  /// グリッド点における波動関数値配列(std::vector<WFGrid>)2つから
  /// 交換相関項(Fxc)を求める。
  ///
  /// @param [out] pF
  /// 交換相関項行列。この行列は実対称行列オブジェクトを指定すること。
  void buildFock(std::vector<WFGrid>::const_iterator pBegin,
                 std::vector<WFGrid>::const_iterator pEnd,
                 std::vector<WFGrid>::const_iterator qBegin,
                 std::vector<WFGrid>::const_iterator qEnd, const double coef,
                 const double cutoffValue, TlMatrixObject* pF);

  void getPrefactor(int nType, const TlPosition& pos, double* pPrefactor);
  void getPrefactorForDerivative(int nType, double alpha, const TlPosition& pos,
                                 double* pPrefactorX, double* pPrefactorY,
                                 double* pPrefactorZ);
  void getPrefactorForSecondDerivative(const int nType, const double alpha,
                                       const TlPosition& pos, double* pXX,
                                       double* pXY, double* pXZ, double* pYY,
                                       double* pYZ, double* pZZ);

 protected:
  void getAOs(const TlPosition& gridPosition, std::vector<double>* pAO_values);
  void getDAOs(const TlPosition& gridPosition,
               std::vector<double>* p_dAO_dx_values,
               std::vector<double>* p_dAO_dy_values,
               std::vector<double>* p_dAO_dz_values);
  void getD2AOs(const TlPosition& gridPosition,
                std::vector<double>* p_d2AO_dxdx_values,
                std::vector<double>* p_d2AO_dxdy_values,
                std::vector<double>* p_d2AO_dxdz_values,
                std::vector<double>* p_d2AO_dydy_values,
                std::vector<double>* p_d2AO_dydz_values,
                std::vector<double>* p_d2AO_dzdz_values);

 protected:
  void getAOs_core(const TlPosition& gridPosition,
                   const std::vector<index_type>& AO_indeces,
                   std::vector<double>* pAO_values);
  void getDAOs_core(const TlPosition& gridPosition,
                    const std::vector<index_type>& AO_indeces,
                    std::vector<double>* p_dAO_dx_values,
                    std::vector<double>* p_dAO_dy_values,
                    std::vector<double>* p_dAO_dz_values);
  void getD2AOs_core(const TlPosition& gridPosition,
                     const std::vector<index_type>& AO_indeces,
                     std::vector<double>* p_d2AO_dxdx_values,
                     std::vector<double>* p_d2AO_dxdy_values,
                     std::vector<double>* p_d2AO_dxdz_values,
                     std::vector<double>* p_d2AO_dydy_values,
                     std::vector<double>* p_d2AO_dydz_values,
                     std::vector<double>* p_d2AO_dzdz_values);

 protected:
  // for symmetric density matrix
  void getRhoAtGridPoint(const TlMatrixObject& PA,
                         const std::vector<double>& AO_values, double* pRhoA);
  // for asymmetric density matrix
  void getRhoAtGridPoint(const TlMatrixObject& PA,
                         const std::vector<double>& row_AO_values,
                         const std::vector<double>& col_AO_values,
                         double* pRhoA);

  // for symmetric density matrix
  void getGradRhoAtGridPoint(const TlMatrixObject& PA,
                             const std::vector<double>& AO_values,
                             const std::vector<double>& dAO_dx_values,
                             const std::vector<double>& dAO_dy_values,
                             const std::vector<double>& dAO_dz_values,
                             double* pGradRhoAX, double* pGradRhoAY,
                             double* pGradRhoAZ);

  // for asymmetric density matrix
  void getGradRhoAtGridPoint(const TlMatrixObject& PA,
                             const std::vector<double>& row_AO_values,
                             const std::vector<double>& row_dAO_dx_values,
                             const std::vector<double>& row_dAO_dy_values,
                             const std::vector<double>& row_dAO_dz_values,
                             const std::vector<double>& col_AO_values,
                             const std::vector<double>& col_dAO_dx_values,
                             const std::vector<double>& col_dAO_dy_values,
                             const std::vector<double>& col_dAO_dz_values,
                             double* pGradRhoAX, double* pGradRhoAY,
                             double* pGradRhoAZ);

 public:
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     DfFunctional_LDA* pFunctional,
                                     TlSymmetricMatrix* pF_A);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     const TlSymmetricMatrix& P_B,
                                     DfFunctional_LDA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlSymmetricMatrix* pF_B);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     DfFunctional_GGA* pFunctional,
                                     TlSymmetricMatrix* pF_A);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     const TlSymmetricMatrix& P_B,
                                     DfFunctional_GGA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlSymmetricMatrix* pF_B);

  virtual void getWholeDensity(double* pRhoA, double* pRhoB) const;

 protected:
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     DfFunctional_LDA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlMatrix* pGridMatrix);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     const TlSymmetricMatrix& P_B,
                                     DfFunctional_LDA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlSymmetricMatrix* pF_B,
                                     TlMatrix* pGridMatrix);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     DfFunctional_GGA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlMatrix* pGridMatrix);
  double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                     const TlSymmetricMatrix& P_B,
                                     DfFunctional_GGA* pFunctional,
                                     TlSymmetricMatrix* pF_A,
                                     TlSymmetricMatrix* pF_B,
                                     TlMatrix* pGridMatrix);

 protected:
  virtual void calcRho_LDA(const TlSymmetricMatrix& P_A);
  virtual void calcRho_LDA(const TlSymmetricMatrix& P_A,
                           const TlSymmetricMatrix& P_B);
  virtual void calcRho_GGA(const TlSymmetricMatrix& P_A);
  virtual void calcRho_GGA(const TlSymmetricMatrix& P_A,
                           const TlSymmetricMatrix& P_B);

  void calcRho_LDA_part(const TlSymmetricMatrix& P_A, TlMatrix* pGridMat);
  void calcRho_LDA_part(const TlSymmetricMatrix& P_A,
                        const TlSymmetricMatrix& P_B, TlMatrix* pGridMat);
  void calcRho_GGA_part(const TlSymmetricMatrix& P_A, TlMatrix* pGridMat);
  void calcRho_GGA_part(const TlSymmetricMatrix& P_A,
                        const TlSymmetricMatrix& P_B, TlMatrix* pGridMat);

  virtual double buildVxc(DfFunctional_LDA* pFunctional,
                          TlSymmetricMatrix* pF_A);
  virtual double buildVxc(DfFunctional_LDA* pFunctional,
                          TlSymmetricMatrix* pF_A, TlSymmetricMatrix* pF_B);
  virtual double buildVxc(DfFunctional_GGA* pFunctional,
                          TlSymmetricMatrix* pF_A);
  virtual double buildVxc(DfFunctional_GGA* pFunctional,
                          TlSymmetricMatrix* pF_A, TlSymmetricMatrix* pF_B);

  double buildVxc(const TlMatrix& gridMatrix, DfFunctional_LDA* pFunctional,
                  TlMatrixObject* pF_A);
  double buildVxc(const TlMatrix& gridMatrix, DfFunctional_LDA* pFunctional,
                  TlMatrixObject* pF_A, TlMatrixObject* pF_B);
  double buildVxc(const TlMatrix& gridMatrix, DfFunctional_GGA* pFunctional,
                  TlMatrixObject* pF_A);
  double buildVxc(const TlMatrix& gridMatrix, DfFunctional_GGA* pFunctional,
                  TlMatrixObject* pF_A, TlMatrixObject* pF_B);

  void build_XC_Matrix(const double roundF_roundRhoA,
                       const std::vector<double>& AO_values,
                       DfFunctional_LDA* pFunctional, const double weight,
                       TlMatrixObject* pF_A);
  void build_XC_Matrix(const double roundF_roundRhoA,
                       const double roundF_roundGammaAA,
                       const double roundF_roundGammaAB, const double gradRhoAX,
                       const double gradRhoAY, const double gradRhoAZ,
                       const std::vector<double>& AO_values,
                       const std::vector<double>& dAO_dx_values,
                       const std::vector<double>& dAO_dy_values,
                       const std::vector<double>& dAO_dz_value,
                       DfFunctional_GGA* pFunctional, const double weight,
                       TlMatrixObject* pF_A);

 protected:
  double energyGradient_part(const TlSymmetricMatrix& P_A,
                             DfFunctional_LDA* pFunctional,
                             const int startAtomIndex, const int endAtomIndex,
                             TlMatrix* pFxc_f, TlMatrix* pFxc_w);
  double energyGradient_part(const TlSymmetricMatrix& P_A,
                             DfFunctional_GGA* pFunctional,
                             const int startAtomIndex, const int endAtomIndex,
                             TlMatrix* pFxc_f, TlMatrix* pFxc_w);

 protected:
  static const double TOOBIG;
  static const double EPS;
  static const double INV_SQRT3;   // = 1.0 / sqrt(3.0);
  static const double INV_SQRT12;  // = 1.0 / sqrt(12.0);

  DfXCFunctional::FUNCTIONAL_TYPE functionalType_;

  /// 入力された電子密度・勾配のカットオフ値
  double inputtedDensityCutoffValue_;

  /// 各グリッド上の電子密度・密度勾配のカットオフ値(alpha spin)
  /// RKS計算ではこの値を共通して用いる
  double m_densityCutOffValueA;

  /// 各グリッド上の電子密度・密度勾配のカットオフ値(beta spin)
  double m_densityCutOffValueB;

  /// F_xc行列要素のカットオフ値
  double m_inputedCutoffThreshold;

  TlOrbitalInfo m_tlOrbInfo;

  index_type numOfRows_gridMatrix_;
  index_type numOfCols_gridMatrix_;

 protected:
  bool isDebugOutPhiTable_;
  bool isSaveGrad_;
};

#endif  // DFCALCGRIDX_H
