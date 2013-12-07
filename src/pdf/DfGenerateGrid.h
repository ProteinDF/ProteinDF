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
#include "TlPosition.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "Fl_Geometry.h"

/** グリッドを生成するクラス
 */
class DfGenerateGrid : public DfObject {
public:
    DfGenerateGrid(TlSerializeData* pPdfParam);
    virtual ~DfGenerateGrid();

    /** グリッドを作成する
     */
    int dfGrdMain();

protected:
    virtual void makeTable();

    void generateGrid(const TlMatrix& O,
                      const int iatom,
                      std::vector<double>* pCoordX,
                      std::vector<double>* pCoordY,
                      std::vector<double>* pCoordZ,
                      std::vector<double>* pWeight);

    void generateGrid_SG1(const TlMatrix& O,
                          const int iatom,
                          std::vector<double>* pCoordX,
                          std::vector<double>* pCoordY,
                          std::vector<double>* pCoordZ,
                          std::vector<double>* pWeight);

    /// rotation invariant用行列を求める
    /// ref) B. G. Johnson, P. M. W. Gill, J. A. Pople, Chem. Phys. Lett., 220, 377, (1994).
    virtual TlMatrix getOMatrix();

private:

    /** 各原子の半径またはGauss-Legendre abscissas とGauss-Legendre weights を入力する
     */
    void setCellPara();

    /** ファジーセル法を用いて、各グリッドの中心座標とそのグリッドにおける重みベクトルを計算する。
     */
    virtual void generateGrid(const TlMatrix& O);

    /// 球面方向のグリッドを生成する
    void points2(const int nOgrid, const double r0, const TlPosition& core,
                 const double weight, 
                 const TlMatrix& O,
                 std::vector<TlPosition>& Ogrid, std::vector<double>& w);


protected:
    enum GridType {
        COARSE,
        MEDIUM,
        MEDIUM_FINE,
        FINE,
        SG_1
    };

protected:
    std::string gridtype;
    std::string xctype;
    GridType m_gridType;

    int nrgrid; /// grid point number or diagonal
    int nOgrid; /// grid point number of Omega

    const Fl_Geometry flGeometry_;

    /// GridDataFileのパス
    //std::string gridDataFilePath_;

    double maxRadii_;
    std::vector<TlPosition> coord_;
    TlSymmetricMatrix distanceMatrix_;

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
    TlMatrix grdMat_;

    /// グリッド情報行列の列数
    int numOfColsOfGrdMat_;
};

#endif // DFGENERATEGRID_H
