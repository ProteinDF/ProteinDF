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

#ifndef DFLEVELSHIFT_H
#define DFLEVELSHIFT_H

#include <string>
#include "DfObject.h"

/// Fock行列に対してレベルシフト法の処理を行うクラス
class DfLevelshift : public DfObject {
   public:
    DfLevelshift(TlSerializeData* pPdfParam, int num_iter);
    virtual ~DfLevelshift();

   public:
    void DfLshiftMain();

    /// フラグメントfragname のF 行列に対してレベルシフトの処置を施す
    void DfLshiftQclo(const std::string& fragname, int norbcut);

   private:
    void main(const RUN_TYPE runType, int iteration,
              const std::string& fragname = "", bool bPdfQcloMode = false);
};

#endif  // DFLEVELSHIFT_H
