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

#ifndef DFCALCGRID_H
#define DFCALCGRID_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"
#include "TlPosition.h"
#include "tl_dense_vector_lapack.h"

/** 作成されたグリッドから交換相関ポテンシャル計算を行うクラス
 */
class DfCalcGrid : DfObject {
   public:
    DfCalcGrid(TlSerializeData* pPdfParam, int niter);
    virtual ~DfCalcGrid();

    /** グリッドの重みの計算を行う
     */
    int dfGrdMain();

   private:
    /** rhoとrhoの微分を計算し、グリッド積分の計算を行い
     */
    void calcXCInteg(TlDenseVector_Lapack& tmpVectorA,
                     TlDenseVector_Lapack& tmpVectorB,
                     TlDenseVector_Lapack& eTmpVector);

    /** calcXCInteg()の結果に<gamma delta>-1
     * をかけて交換相関ポテンシャルの展開係数を求める(RKS版)
     */
    void calcXCcoef_RKS(const TlDenseVector_Lapack& tmpVector,
                        const TlDenseVector_Lapack& eTmpVector);

    /** calcXCInteg()の結果に<gamma delta>-1
     * をかけて交換相関ポテンシャルの展開係数を求める(UKS版)
     */
    void calcXCcoef_UKS(const TlDenseVector_Lapack& tmpVectorA,
                        const TlDenseVector_Lapack& tmpVectorB,
                        const TlDenseVector_Lapack& eTmpVector);

   private:
    double getPrefactor(int nType, const TlPosition& pos);
    void getPrefactorForDerivative(int nType, double alpha,
                                   const TlPosition& pos, double* pPrefactorX,
                                   double* pPrefactorY, double* pPrefactorZ);

   private:
    void calcXCIntegRhoTilde_RKS(const TlDenseVector_Lapack& RhoAlphaA,
                                 TlDenseGeneralMatrix_Lapack* pGridMat);
    void calcXCIntegRhoTilde_UKS(const TlDenseVector_Lapack& RhoAlphaA,
                                 const TlDenseVector_Lapack& RhoAlphaB,
                                 TlDenseGeneralMatrix_Lapack* pGridMat);

    void calcXCIntegMyuEpsilon_RKS(const TlDenseGeneralMatrix_Lapack& gridMat,
                                   TlDenseVector_Lapack& tmpVectorA,
                                   TlDenseVector_Lapack& eTmpVector);
    void calcXCIntegMyuEpsilon_UKS(const TlDenseGeneralMatrix_Lapack& gridMat,
                                   TlDenseVector_Lapack& tmpVectorA,
                                   TlDenseVector_Lapack& tmpVectorB,
                                   TlDenseVector_Lapack& eTmpVector);

   private:
    double polfunc(double z);
    double poldrfunc(double z);

    double HLfunc(double z);
    double HLdrfunc(double z);

    double VWNPfunc(double z);
    double VWNPdrfunc(double z);
    double VWNFfunc(double z);
    double VWNFdrfunc(double z);

    double B88func(double RouA, double RouB, double XA, double XB);
    double DB88func(double RouA, double RouB, double XA, double XB);

    double B88dfunc(double Rou, double X, double G, double gRx, double gRy,
                    double gRz, double g, double ggx, double ggy, double ggz);
    double DB88dfunc(double Rou, double X, double G, double gRx, double gRy,
                     double gRz, double g, double ggx, double ggy, double ggz);

    double G96func(double RouA, double RouB, double XA, double XB);
    double DG96func(double RouA, double RouB, double XA, double XB);

    double G96dfunc(double Rou, double X, double G, double gRx, double gRy,
                    double gRz, double g, double ggx, double ggy, double ggz);
    double DG96dfunc(double Rou, double X, double G, double gRx, double gRy,
                     double gRz, double g, double ggx, double ggy, double ggz);

    double LYPfunc(double TRou, double RouA, double RouB, double GAA,
                   double GAB, double GBB);
    double LYPdfunc(double TRou, double RouA, double RouB, double GAA,
                    double GAB, double GBB, double gRAx, double gRAy,
                    double gRAz, double gRBx, double gRBy, double gRBz,
                    double g, double ggx, double ggy, double ggz);

   private:
    std::string gridDataFilePath_;

    double alphaval;

    double* gDelta;  // BF of Myu

    int xc;         // the through-number of XC function type
    int nlsd_type;  // flag for NLSD type
    int tilude;

    int vectorelement;

    int TotEleNum;  // total element number of workarea

    TlOrbitalInfo m_tlOrbInfo;
    TlOrbitalInfo_Density m_tlOrbInfoAuxCD_;
    TlOrbitalInfo_XC m_tlOrbInfoXC_;
};

#endif  // DFCALCGRID_H
