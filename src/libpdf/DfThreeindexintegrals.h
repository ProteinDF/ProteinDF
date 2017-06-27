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
