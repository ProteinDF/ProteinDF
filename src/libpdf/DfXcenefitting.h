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

#ifndef DFXCENEFITTING_H
#define DFXCENEFITTING_H

#include <string>
#include "DfObject.h"

class TlVector;

/**
 * 交換相関ポテンシャルとしてX-alpha法を使用時における交換相関ポテンシャルの展開係数を求めるクラス
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

#endif  // DFXCENEFITTING_H
