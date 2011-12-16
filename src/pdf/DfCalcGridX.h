#ifndef DFCALCGRIDX_H
#define DFCALCGRIDX_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "DfFunctional.h"
#include "TlOrbitalInfo.h"
#include "GridDataManager.h"

class TlMatrixObject;

class DfCalcGridX : public DfObject {
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
    /// Fockの交換相関項と、エネルギーを同時に求める(RKS, LDA用)
    virtual double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P,
                                               DfFunctional_LDA* pFunctional,
                                               TlSymmetricMatrix* pF);

    /// Fockの交換相関項と、エネルギーを同時に求める(UKS, LDA用)
    virtual double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& PA,
                                               const TlSymmetricMatrix& PB,
                                               DfFunctional_LDA* pFunctional,
                                               TlSymmetricMatrix* pFA,
                                               TlSymmetricMatrix* pFB);

    /// Fockの交換相関項と、エネルギーを同時に求める(RKS, GGA用)
    virtual double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P,
                                               DfFunctional_GGA* pFunctional,
                                               TlSymmetricMatrix* pF);

    /// Fockの交換相関項と、エネルギーを同時に求める(UKS, GGA用)
    virtual double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& PA,
                                               const TlSymmetricMatrix& PB,
                                               DfFunctional_GGA* pFunctional,
                                               TlSymmetricMatrix* pFA,
                                               TlSymmetricMatrix* pFB);

public:
    // for density field data
    void gridDensity(const TlSymmetricMatrix& P,
                     const TlPosition& gridPosition,
                     double* pRho);

public:
    // for derivative
    void calcForceFromXC(const TlSymmetricMatrix& P,
                         DfFunctional_GGA* pFunctional,
                         TlVector* pFx, TlVector* pFy, TlVector* pFz);

    void makeGammaMatrix(const TlSymmetricMatrix& P,
                         DfFunctional_LDA* pFunctional,
                         TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ);

    void makeGammaMatrix(const TlSymmetricMatrix& P,
                         DfFunctional_GGA* pFunctional,
                         TlMatrix* pGX, TlMatrix* pGY, TlMatrix* pGZ);
    
protected:
    virtual void defineCutOffValues(const TlSymmetricMatrix& P);

    virtual void defineCutOffValues(const TlSymmetricMatrix& PA,
                                    const TlSymmetricMatrix& PB);

    double calcXCIntegForFockAndEnergy1(const TlSymmetricMatrix& P,
                                        DfFunctional_LDA* pFunctional,
                                        TlSymmetricMatrix* pF);

    double calcXCIntegForFockAndEnergy1(const TlSymmetricMatrix& P,
                                        DfFunctional_GGA* pFunctional,
                                        TlSymmetricMatrix* pF);

    /// Fockの交換相関項と、エネルギーを同時に求める(RKS, LDA用; 並列計算と共通部分)
    double calcXCIntegForFockAndEnergy(int nStartAtom, int nEndAtom,
                                       const TlSymmetricMatrix& P,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pF);

    /// Fockの交換相関項と、エネルギーを同時に求める(UKS, LDA用; 並列計算と共通部分)
    double calcXCIntegForFockAndEnergy(int nStartAtom, int nEndAtom,
                                       const TlSymmetricMatrix& PA,
                                       const TlSymmetricMatrix& PB,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pFA,
                                       TlSymmetricMatrix* pFB);

    /// Fockの交換相関項と、エネルギーを同時に求める(RKS, GGA用; 並列計算と共通部分)
    double calcXCIntegForFockAndEnergy(int nStartAtom, int nEndAtom,
                                       const TlSymmetricMatrix& P,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pF);

    /// Fockの交換相関項と、エネルギーを同時に求める(UKS, GGA用; 並列計算と共通部分)
    double calcXCIntegForFockAndEnergy(int nStartAtom, int nEndAtom,
                                       const TlSymmetricMatrix& PA,
                                       const TlSymmetricMatrix& PB,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pFA,
                                       TlSymmetricMatrix* pFB);

    virtual void backupGridData();
    virtual void flushGridData();


protected:
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


protected:
    double getPrefactor(int nType, const TlPosition& pos);
    void getPrefactorForDerivative(int nType, double alpha, const TlPosition& pos,
                                   double* pPrefactorX, double* pPrefactorY, double* pPrefactorZ);

protected:
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
    // experimental code -------------------------------------------------------
    double calcXCIntegForFockAndEnergy_usemat(const TlSymmetricMatrix& P,
                                              DfFunctional_GGA* pFunctional,
                                              TlSymmetricMatrix* pF,
                                              TlMatrix* pGridMatrix);
    double calcXCIntegForFockAndEnergy3(const TlSymmetricMatrix& P,
                                        DfFunctional_GGA* pFunctional,
                                        TlSymmetricMatrix* pF);

    double calcXCIntegForFockAndEnergy2(const TlSymmetricMatrix& P,
                                        DfFunctional_LDA* pFunctional,
                                        TlSymmetricMatrix* pF);
    double calcXCIntegForFockAndEnergy2(const TlSymmetricMatrix& P,
                                        DfFunctional_GGA* pFunctional,
                                        TlSymmetricMatrix* pF);
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
    
    
private:
    void readTable();

protected:
    static const double TOOBIG;
    static const double EPS;
    static const double INV_SQRT3; // = 1.0 / sqrt(3.0);
    static const double INV_SQRT12; // = 1.0 / sqrt(12.0);

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

    /// [データ型][原子][グリッド点index]
    std::map<GridDataManager::ChunkType,
    std::map<int, std::vector<double> > > physicalValues_;

protected:
    bool isDebugOutPhiTable_;
};




#endif //DFCALCGRIDX_H
