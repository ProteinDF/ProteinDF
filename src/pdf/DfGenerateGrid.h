#ifndef DFGENERATEGRID_H
#define DFGENERATEGRID_H

#include <string>
#include "DfObject.h"
#include "TlPosition.h"
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

    void generateGrid(int iatom,
                      std::vector<double>* pCoordX,
                      std::vector<double>* pCoordY,
                      std::vector<double>* pCoordZ,
                      std::vector<double>* pWeight);

    void generateGrid_SG1(int iatom,
                          std::vector<double>* pCoordX,
                          std::vector<double>* pCoordY,
                          std::vector<double>* pCoordZ,
                          std::vector<double>* pWeight);

private:

    /** 各原子の半径またはGauss-Legendre abscissas とGauss-Legendre weights を入力する
     */
    void setCellPara();

    /** ファジーセル法を用いて、各グリッドの中心座標とそのグリッドにおける重みベクトルを計算する。
     */
    virtual void generateGrid();

    void points2(const int nOgrid, const double r0, const TlPosition& core,
                 const double weight, std::vector<TlPosition>& Ogrid, std::vector<double>& w);

protected:
    enum GridType {
        COARSE,
        MEDIUM,
        MEDIUM_FINE,
        FINE,
        SG_1
    };

    //int outlevel;
    //int pout;

    std::string gridtype;
//     std::string nyuswtch;

    std::string xctype;

    GridType m_gridType;

    int nrgrid; /// grid point number or diagonal
    int nOgrid; /// grid point number of Omega

protected:
    const Fl_Geometry flGeometry_;

    /// GridDataFileのパス
    std::string gridDataFilePath_;

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
};

#endif // DFGENERATEGRID_H
