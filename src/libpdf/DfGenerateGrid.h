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

#ifndef DFGENERATEGRID_H
#define DFGENERATEGRID_H

#include <string>
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlLebedevGrid.h"
#include "TlPosition.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

/** グリッドを生成するクラス
 */
class DfGenerateGrid : public DfObject {
   public:
    DfGenerateGrid(TlSerializeData* pPdfParam);
    virtual ~DfGenerateGrid();

    /** グリッドを作成する
     */
    int dfGrdMain();

    void getGrids(const int atomIndex, std::vector<TlPosition>* pGrids,
                  std::vector<double>* pSingleCenterWeights,
                  std::vector<double>* pPartitioningWeights);

   protected:
    virtual void makeTable();

    void generateGrid(const TlDenseGeneralMatrix_Lapack& O, const int iatom,
                      std::vector<double>* pCoordX,
                      std::vector<double>* pCoordY,
                      std::vector<double>* pCoordZ,
                      std::vector<double>* pWeight);

    void generateGrid_SG1(const TlDenseGeneralMatrix_Lapack& O, const int iatom,
                          std::vector<double>* pCoordX,
                          std::vector<double>* pCoordY,
                          std::vector<double>* pCoordZ,
                          std::vector<double>* pWeight);

    /// rotation invariant用行列を求める
    /// ref) B. G. Johnson, P. M. W. Gill, J. A. Pople, Chem. Phys. Lett., 220,
    /// 377, (1994).
    virtual TlDenseGeneralMatrix_Lapack getOMatrix();

   private:
    /** 各原子の半径またはGauss-Legendre abscissas とGauss-Legendre weights
     * を入力する
     */
    void setCellPara();

    /**
     * ファジーセル法を用いて、各グリッドの中心座標とそのグリッドにおける重みベクトルを計算する。
     */
    virtual void generateGrid(const TlDenseGeneralMatrix_Lapack& O);

    /// 球面方向のグリッドを生成する
    void points2(const int nOgrid, const double r0, const TlPosition& core,
                 const double weight, const TlDenseGeneralMatrix_Lapack& O,
                 std::vector<TlPosition>& Ogrid, std::vector<double>& w);

   private:
    enum RADIAL_GRID_TYPE { RG_EularMaclaurin, RG_GaussChebyshev };

    enum GC_MAPPING_TYPE { GC_BECKE, GC_TA, GC_KK };

    enum PARTITIONING_METHOD { Paritioning_Becke, Partitioning_SSWeight };

   private:
    ///
    /// @param [in] R atomic size
    /// @param [in] Nr the number of radial grids
    /// @param [in] i index (1 <= i <= Nr)
    /// @param [out] p_ri pointer to grid point(r_i)
    /// @param [out] p_Weight pointer to weight
    void getRadialAbscissaAndWeight_EulerMaclaurin(const double R,
                                                   const double Nr, const int i,
                                                   double* p_ri,
                                                   double* pWeight);

    ///
    /// @param [in] R atomic size
    /// @param [in] Nr the number of radial grids
    /// @param [in] i index (1 <= i <= Nr)
    /// @param [out] p_ri pointer to grid point(r_i)
    /// @param [out] p_Weight pointer to weight
    void getRadialAbscissaAndWeight_GaussChebyshev(const double R, const int Nr,
                                                   const int i, double* p_ri,
                                                   double* pWeight);

    int getNumOfPrunedAnglarPoints_SG1(const double r, const double inv_R,
                                       const std::vector<double>& alpha);

    int getNumOfPrunedAnglarPoints(const double r, const int maxNumOfAngGrids,
                                   const int atomicNumber);

    void getSphericalGrids(const int numOfGrids, const double r,
                           const double radial_weight, const TlPosition& center,
                           const TlDenseGeneralMatrix_Lapack& O,
                           std::vector<TlPosition>* pGrids,
                           std::vector<double>* pWeights);

    void calcMultiCenterWeight_Becke(const int iAtom, const int Ogrid,
                                     const std::vector<TlPosition>& points,
                                     std::vector<double>* pWeight);
    void calcMultiCenterWeight_SS(const int iAtom, const int Ogrid,
                                  const std::vector<TlPosition>& points,
                                  std::vector<double>* pWeight);

    double getCovalentRadiiForBecke(const int atomicNumber);
    double Becke_f1(const double x);
    double Becke_f3(const double x);

    void screeningGridsByWeight0(std::vector<TlPosition>* pGrids,
                                 std::vector<double>* pWeights);
    void screeningGridsByWeight(std::vector<TlPosition>* pGrids,
                                std::vector<double>* pSingleCenterWeights,
                                std::vector<double>* pPartitioningWeights);

    double Ps_uij(const int atomIndex_m, const TlPosition& gridpoint);

   public:
    TlDenseVector_Lapack JGP_nablaB_omegaA(const int atomIndexA,
                                           const int atomIndexB,
                                           const TlPosition& gridpoint);

   private:
    TlDenseVector_Lapack JGP_nablaA_PA(const int atomIndexA,
                                       const TlPosition& gridpoint);
    TlDenseVector_Lapack JGP_nablaB_PA(const int atomIndexA,
                                       const int atonIndexB,
                                       const TlPosition& gridpoint);
    double JGP_t(const double myu);
    TlDenseVector_Lapack JGP_nablaA_myuAB(const int atomIndexA,
                                          const int atomIndexB,
                                          const TlPosition& gridpoint);

    void getGrids_sub(const int atomIndex, const int gridType,
                      const int numOfRadialGrids, const int radialGridType,
                      const int maxAngularGrids,
                      const PARTITIONING_METHOD partitioningMethod,
                      std::vector<TlPosition>* pGrids,
                      std::vector<double>* pSingleCenterWeights,
                      std::vector<double>* pPartitioningWeights);

   protected:
    enum GridType { COARSE, MEDIUM, MEDIUM_FINE, FINE, SG_1, USER };

   protected:
    std::string gridtype;
    std::string xctype;
    GridType m_gridType;

    int nrgrid;  /// grid point number or diagonal
    int nOgrid;  /// grid point number of Omega

    const Fl_Geometry flGeometry_;

    /// GridDataFileのパス
    // std::string gridDataFilePath_;

    double maxRadii_;
    std::vector<TlPosition> coord_;
    TlDenseSymmetricMatrix_Lapack distanceMatrix_;
    TlDenseSymmetricMatrix_Lapack invDistanceMatrix_;

    /// Radius list of each atom
    std::vector<double> radiusList_;

    /// Gauss-Legendre abscissas
    std::vector<double> xGL_;

    // Gauss-Legendre weights
    std::vector<double> wGL_;

    double weightCutoff_;

    /// グリッド情報行列
    /// 行方向: grid
    /// 0, 1, 2列: x, y, z
    /// 3列: weight
    TlDenseGeneralMatrix_Lapack grdMat_;

    /// グリッド情報行列の列数
    int numOfColsOfGrdMat_;

    ///
    RADIAL_GRID_TYPE radialGridType_;

    /// Gauss-Chebyshev mapping type
    GC_MAPPING_TYPE GC_mappingType_;

    /// using atomic size adjustments in Becke partitioning
    bool isAtomicSizeAdjustments_;

    bool isPruning_;

    TlLebedevGrid lebGrd_;

    /// rotation matrix
    TlDenseGeneralMatrix_Lapack O_;
};

#endif  // DFGENERATEGRID_H
