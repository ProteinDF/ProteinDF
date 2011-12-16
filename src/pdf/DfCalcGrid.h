#ifndef DFCALCGRID_H
#define DFCALCGRID_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"
#include "GridDataManager.h"
#include "TlVector.h"
#include "TlPosition.h"

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
    void calcXCInteg(TlVector& tmpVectorA, TlVector& tmpVectorB, TlVector& eTmpVector);

    /** calcXCInteg()の結果に<gamma delta>-1 をかけて交換相関ポテンシャルの展開係数を求める(RKS版)
     */
    void calcXCcoef_RKS(const TlVector& tmpVector, const TlVector& eTmpVector);

    /** calcXCInteg()の結果に<gamma delta>-1 をかけて交換相関ポテンシャルの展開係数を求める(UKS版)
     */
    void calcXCcoef_UKS(const TlVector& tmpVectorA, const TlVector& tmpVectorB, const TlVector& eTmpVector);

private:
    double getPrefactor(int nType, const TlPosition& pos);
    void getPrefactorForDerivative(int nType, double alpha, const TlPosition& pos,
                                   double* pPrefactorX, double* pPrefactorY, double* pPrefactorZ);


private:
    void calcXCIntegRho(const std::vector<GridDataManager::GridInfo>& grids,
                        const TlSymmetricMatrix& PA, const TlSymmetricMatrix& PB,
                        TlVector& gridRhoA,  TlVector& gridRhoB,
                        TlVector& gradRhoAx, TlVector& gradRhoAy, TlVector& gradRhoAz,
                        TlVector& gradRhoBx, TlVector& gradRhoBy, TlVector& gradRhoBz);
    void calcXCIntegRhoTilde_RKS(const std::vector<GridDataManager::GridInfo>& grids,
                                 const TlVector& RhoAlphaA,
                                 TlVector& gridRhoA);
    void calcXCIntegRhoTilde_UKS(const std::vector<GridDataManager::GridInfo>& grids,
                                 const TlVector& RhoAlphaA, const TlVector& RhoAlphaB,
                                 TlVector& gridRhoA, TlVector& gridRhoB);

    void calcXCIntegGradRhoTilde(const std::vector<GridDataManager::GridInfo>& grids,
                                 const TlVector& RhoAlphaA, const TlVector& RhoAlphaB,
                                 TlVector& gradRhoAx, TlVector& gradRhoAy, TlVector& gradRhoAz,
                                 TlVector& gradRhoBx, TlVector& gradRhoBy, TlVector& gradRhoBz);

    void calcXCIntegMyuEpsilon_RKS(const std::vector<GridDataManager::GridInfo>& grids,
                                   const TlVector& gridRhoA,
                                   TlVector& tmpVectorA, TlVector& eTmpVector);
    void calcXCIntegMyuEpsilon_UKS(const std::vector<GridDataManager::GridInfo>& grids,
                                   const TlVector& gridRhoA, const TlVector& gridRhoB,
                                   TlVector& tmpVectorA, TlVector& tmpVectorB, TlVector& eTmpVector);

    void calcXCIntegLyp(const std::vector<GridDataManager::GridInfo>& grids,
                        const TlVector& gridRhoA, const TlVector& gridRhoB,
                        const TlVector& gradRhoAx, const TlVector& gradRhoAy, const TlVector& gradRhoAz,
                        const TlVector& gradRhoBx, const TlVector& gradRhoBy, const TlVector& gradRhoBz,
                        TlVector& tmpVectorA, TlVector& tmpVectorB, TlVector& eTmpVector);

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

    double B88dfunc(double Rou, double X, double G, double gRx, double gRy, double gRz, double g, double ggx, double ggy, double ggz);
    double DB88dfunc(double Rou, double X, double G, double gRx, double gRy, double gRz, double g, double ggx, double ggy, double ggz);

    double G96func(double RouA, double RouB, double XA, double XB);
    double DG96func(double RouA, double RouB, double XA, double XB);

    double G96dfunc(double Rou, double X, double G, double gRx, double gRy, double gRz, double g, double ggx, double ggy, double ggz);
    double DG96dfunc(double Rou, double X, double G, double gRx, double gRy, double gRz, double g, double ggx, double ggy, double ggz);

    double LYPfunc(double TRou, double RouA, double RouB, double GAA, double GAB, double GBB);
    double LYPdfunc(double TRou, double RouA, double RouB, double GAA, double GAB, double GBB,
                    double gRAx, double gRAy, double gRAz, double gRBx, double gRBy, double gRBz,
                    double g, double ggx, double ggy, double ggz);

private:
    std::string gridDataFilePath_;

    double alphaval;

    double* gDelta;                     // BF of Myu

    int xc;                             // the through-number of XC function type
    int nlsd_type;                      // flag for NLSD type
    int tilude;

    int vectorelement;

    int TotEleNum;                      // total element number of workarea

    TlOrbitalInfo m_tlOrbInfo;
    TlOrbitalInfo_Density m_tlOrbInfoAuxCD_;
    TlOrbitalInfo_XC m_tlOrbInfoXC_;
};

#endif // DFCALCGRID_H
