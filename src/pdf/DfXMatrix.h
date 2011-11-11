#ifndef DFXMATRIX_H
#define DFXMATRIX_H

#include "DfObject.h"
#include "TlMatrix.h"
#include "TlVector.h"

/// X行列を求めるクラス
/// S行列から、固有値, 固有ベクトルを求め、X 行列および−1 X 行列の計算を行い、
/// fl_Mtr_X.matrix, fl_Mtr_InvX.matrixを出力する
class DfXMatrix : public DfObject {
public:
    // flGbi は書き換えられる
    DfXMatrix(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix();

public:
    virtual void main();

protected:
    void exec();

    virtual void saveNumOfIndependentBasis();

    void checkMatrixes();

protected:
    /// 一次従属性の判定値
    double threshold_trancation;
};


#endif // DFXMATRIX_H

