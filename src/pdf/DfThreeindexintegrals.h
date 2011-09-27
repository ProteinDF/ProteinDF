#ifndef DFTHREEINDEXINTEGRALS_H
#define DFTHREEINDEXINTEGRALS_H

#include <string>

#include "DfObject.h"
#include "TlSymmetricMatrix.h"

/** フォック行列計算および全エネルギー計算における3中心積分部分の計算を行うクラス
 */
class DfThreeindexintegrals : public DfObject {
public:
    DfThreeindexintegrals(TlSerializeData* pPdfParam);
    ~DfThreeindexintegrals();

    void DfThreeindexintegralsMain();

private:
    void mainDIRECT_RKS(int iteration);
    void mainDIRECT_RKS2(int iteration);
    void mainDIRECT_UKS(int iteration);
    void mainDIRECT_ROKS(int iteration);

private:
    TlSymmetricMatrix getPMatrix(RUN_TYPE runType, int iteration);
};

#endif // DFTHREEINDEXINTEGRALS_H
