#ifndef DFXCENEFITTING_H
#define DFXCENEFITTING_H

#include <string>
#include "DfObject.h"

class TlVector;

/** 交換相関ポテンシャルとしてX-alpha法を使用時における交換相関ポテンシャルの展開係数を求めるクラス
 */
class DfXcenefitting : public DfObject {
public:
    DfXcenefitting(TlSerializeData* pPdfParam, int nIter);
    virtual ~DfXcenefitting();

    int dfXceMain();

private:
    /** X-alpha法の時の展開係数の計算を行う
     */
    int calcEpsilon();

private:
//     int outlevel;

//     std::string scftype;
//     bool bMemorySave;

    /// 現在のiteration 数
//     int number_iteration;

    /// 交換相関ポテンシャルの補助基底の総数
//     int naux;
};

#endif // DFXCENEFITTING_H

