#ifndef DFCALCGRIDX_H
#define DFCALCGRIDX_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "DfXCFunctional.h"
#include "DfFunctional.h"
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
        GM_LDA_RHO_BETA  = 6,
        // for GGA
        GM_GGA_RHO_ALPHA = 5,
        GM_GGA_GRAD_RHO_X_ALPHA =  6,
        GM_GGA_GRAD_RHO_Y_ALPHA =  7,
        GM_GGA_GRAD_RHO_Z_ALPHA =  8,
        GM_GGA_RHO_BETA = 9,
        GM_GGA_GRAD_RHO_X_BETA  = 10,
        GM_GGA_GRAD_RHO_Y_BETA  = 12,
        GM_GGA_GRAD_RHO_Z_BETA  = 13
    };
    
public:
    // the value of wave function on the grid
    struct WFGrid {
        WFGrid(std::size_t i =0, double v =0.0) : index(i), value(v) {
        }

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
    void gridDensity(const TlSymmetricMatrix& P,
                     const TlPosition& gridPosition,
                     double* pRho);

public:
    // for derivative
    void makeGammaMatrix(const TlSymmetricMatrix& P,
                         DfFunctional_LDA* pFunctional,
                         TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ);

    void makeGammaMatrix(const TlSymmetricMatrix& P,
                         DfFunctional_GGA* pFunctional,
                         TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ);

protected:
    // void calcForceFromXC(const TlSymmetricMatrix& P,
    //                      DfFunctional_GGA* pFunctional,
    //                      TlVector* pFx, TlVector* pFy, TlVector* pFz);

    TlMatrix selectGridMatrixByAtom(const TlMatrix& globalGridMat,
                                    const int atomIndex);
    
protected:
    virtual void defineCutOffValues(const TlSymmetricMatrix& P);

    virtual void defineCutOffValues(const TlSymmetricMatrix& PA,
                                    const TlSymmetricMatrix& PB);

    void buildFock(double dRhoA,
                   const std::vector<WFGrid>& aPhi,
                   DfFunctional_LDA* pFunctional, double dWeight,
                   TlMatrixObject* pF);

    void buildFock(double dRhoA, double dRhoB,
                   const std::vector<WFGrid>& aPhi,
                   DfFunctional_LDA* pFunctional, double dWeight,
                   TlMatrixObject* pFA, TlMatrixObject* FB);

    void buildFock(double dRhoA,
                   double dGradRhoAX, double dGradRhoAY, double dGradRhoAZ,
                   const std::vector<WFGrid>& aPhi,
                   const std::vector<WFGrid>& aGradPhiX,
                   const std::vector<WFGrid>& aGradPhiY,
                   const std::vector<WFGrid>& aGradPhiZ,
                   DfFunctional_GGA* pFunctional, double dWeight,
                   TlMatrixObject* pF);

    void buildFock(double dRhoA,
                   double dGradRhoAX, double dGradRhoAY, double dGradRhoAZ,
                   double dRhoB,
                   double dGradRhoBX, double dGradRhoBY, double dGradRhoBZ,
                   const std::vector<WFGrid>& aPhi,
                   const std::vector<WFGrid>& aGradPhiX,
                   const std::vector<WFGrid>& aGradPhiY,
                   const std::vector<WFGrid>& aGradPhiZ,
                   DfFunctional_GGA* pFunctional, double dWeight,
                   TlMatrixObject* pFA, TlMatrixObject* pFB);

    /// グリッド点における波動関数値配列(std::vector<WFGrid>)1つから
    /// 交換相関項(Fxc)を求める。
    ///
    /// @param [out] pF 交換相関項行列。この行列は実対称行列オブジェクトを指定すること。
    void buildFock(std::vector<WFGrid>::const_iterator pBegin,
                   std::vector<WFGrid>::const_iterator pEnd,
                   double coef, double cutoffValue,
                   TlMatrixObject* pFxc);

    /// グリッド点における波動関数値配列(std::vector<WFGrid>)2つから
    /// 交換相関項(Fxc)を求める。
    ///
    /// @param [out] pF 交換相関項行列。この行列は実対称行列オブジェクトを指定すること。
    void buildFock(std::vector<WFGrid>::const_iterator pBegin,
                   std::vector<WFGrid>::const_iterator pEnd,
                   std::vector<WFGrid>::const_iterator qBegin,
                   std::vector<WFGrid>::const_iterator qEnd,
                   const double coef, const double cutoffValue,
                   TlMatrixObject* pF);

    double getPrefactor(int nType, const TlPosition& pos);
    void getPrefactorForDerivative(int nType, double alpha, const TlPosition& pos,
                                   double* pPrefactorX, double* pPrefactorY, double* pPrefactorZ);

    /// グリッド点における波動関数の値、勾配の値を求める
    void getPhiTable(const TlPosition& pos, std::vector<WFGrid>& aPhi);

    void getPhiTable(const TlPosition& pos, std::vector<WFGrid>& aPhi,
                     std::vector<WFGrid>& aGradPhiX, std::vector<WFGrid>& aGradPhiY, std::vector<WFGrid>& aGradPhiZ);
    void getPhiTable(const TlPosition& pos,
                     int startOrbIndex, int endOrbIndex,
                     std::vector<WFGrid>& aPhi);
    void getPhiTable(const TlPosition& pos,
                     int startOrbIndex, int endOrbIndex,
                     std::vector<WFGrid>& aPhi,
                     std::vector<WFGrid>& aGradPhiX, std::vector<WFGrid>& aGradPhiY, std::vector<WFGrid>& aGradPhiZ);


    /// グリッド点における波動関数の値、勾配の値を求める
    void getPhiTable(const TlPosition& gridPosition,
                     const std::vector<int>& AO_list,
                     std::vector<WFGrid>* pPhis);
    void getPhiTable(const TlPosition& gridPosition,
                     const index_type* pAO_List,
                     const std::size_t AO_ListSize,
                     std::vector<WFGrid>* pPhis,
                     std::vector<WFGrid>* pGradPhiXs,
                     std::vector<WFGrid>* pGradPhiYs,
                     std::vector<WFGrid>* pGradPhiZs);

    
    void getRhoAtGridPoint(const TlMatrixObject& P,
                           const std::vector<WFGrid>& aPhi, double* pGridRhoA);

    void getRhoAtGridPoint(const TlMatrixObject& P, const std::vector<WFGrid>& aPhi,
                           const std::vector<WFGrid>& aGradPhiX, const std::vector<WFGrid>& aGradPhiY,
                           const std::vector<WFGrid>& aGradPhiZ,
                           double* pRhoA, double* pGradRhoAX, double* pGradRhoAY, double* pGradRhoAZ);

protected:
    void build_XC_Matrix(const double roundF_roundRhoA,
                         const std::vector<WFGrid>& phi,
                         DfFunctional_LDA* pFunctional, const double weight,
                         TlMatrixObject* pF_A);
    void build_XC_Matrix(const double roundF_roundRhoA,
                         const double roundF_roundGammaAA,
                         const double roundF_roundGammaAB,
                         const double gradRhoAX, const double gradRhoAY, const double gradRhoAZ,
                         const std::vector<WFGrid>& phi,
                         const std::vector<WFGrid>& gradPhiX,
                         const std::vector<WFGrid>& gradPhiY,
                         const std::vector<WFGrid>& gradPhiZ,
                         DfFunctional_GGA* pFunctional,
                         const double weight,
                         TlMatrixObject* pF_A);

    // experimental code -------------------------------------------------------
    void calcRhoVals_LDA(const std::vector<index_type>& P_rowIndexes,
                         const std::vector<index_type>& P_colIndexes,
                         const TlMatrix& P,
                         TlMatrix* pGridMatrix);
    void calcRhoVals_GGA(const std::vector<index_type>& P_rowIndexes,
                         const std::vector<index_type>& P_colIndexes,
                         const TlMatrix& P,
                         TlMatrix* pGridMatrix);
    void getWaveFunctionValues(const std::vector<index_type>& AO_indexes,
                               const TlPosition& gridPosition,
                               TlMatrix* pWF);
    void getWaveFunctionValues(const std::vector<index_type>& AO_indexes,
                               const TlPosition& gridPosition,
                               TlMatrix* pWF,
                               TlMatrix* pGradWF_X,
                               TlMatrix* pGradWF_Y,
                               TlMatrix* pGradWF_Z);
    double buildK(const TlMatrix& gridMatrix,
                  DfFunctional_LDA* pFunctional,
                  TlSymmetricMatrix* pF);
    double buildK(const TlMatrix& gridMatrix,
                  DfFunctional_GGA* pFunctional,
                  TlSymmetricMatrix* pF);
    double buildK_2(const TlMatrix& gridMatrix,
                    const std::vector<index_type>& rowIndexes,
                    const std::vector<index_type>& colIndexes,
                    DfFunctional_LDA* pFunctional,
                    TlSymmetricMatrix* pF);
    double buildK_2(const TlMatrix& gridMatrix,
                    const std::vector<index_type>& rowIndexes,
                    const std::vector<index_type>& colIndexes,
                    DfFunctional_GGA* pFunctional,
                    TlMatrix* pF);

    // NEW IMPLIMENT -----------------------------------------------------------
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

    void calcRho_LDA(const TlSymmetricMatrix& P_A,
                     TlMatrix* pGridMat);
    void calcRho_LDA(const TlSymmetricMatrix& P_A,
                     const TlSymmetricMatrix& P_B,
                     TlMatrix* pGridMat);
    void calcRho_GGA(const TlSymmetricMatrix& P_A,
                     TlMatrix* pGridMat);
    void calcRho_GGA(const TlSymmetricMatrix& P_A,
                     const TlSymmetricMatrix& P_B,
                     TlMatrix* pGridMat);

    virtual double buildVxc(DfFunctional_LDA* pFunctional,
                            TlSymmetricMatrix* pF_A);
    virtual double buildVxc(DfFunctional_LDA* pFunctional,
                            TlSymmetricMatrix* pF_A,
                            TlSymmetricMatrix* pF_B);
    virtual double buildVxc(DfFunctional_GGA* pFunctional,
                            TlSymmetricMatrix* pF_A);
    virtual double buildVxc(DfFunctional_GGA* pFunctional,
                            TlSymmetricMatrix* pF_A,
                            TlSymmetricMatrix* pF_B);

    double buildVxc(const TlMatrix& gridMatrix,
                    DfFunctional_LDA* pFunctional,
                    TlMatrixObject* pF_A);
    double buildVxc(const TlMatrix& gridMatrix,
                    DfFunctional_LDA* pFunctional,
                    TlMatrixObject* pF_A,
                    TlMatrixObject* pF_B);
    double buildVxc(const TlMatrix& gridMatrix,
                    DfFunctional_GGA* pFunctional,
                    TlMatrixObject* pF_A);
    double buildVxc(const TlMatrix& gridMatrix,
                    DfFunctional_GGA* pFunctional,
                    TlMatrixObject* pF_A,
                    TlMatrixObject* pF_B);
    
protected:
    static const double TOOBIG;
    static const double EPS;
    static const double INV_SQRT3; // = 1.0 / sqrt(3.0);
    static const double INV_SQRT12; // = 1.0 / sqrt(12.0);

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
};


#endif //DFCALCGRIDX_H
