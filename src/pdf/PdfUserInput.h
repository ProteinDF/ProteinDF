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

#ifndef PDFUSERINPUT_H
#define PDFUSERINPUT_H

#include <string>
#include "TlParameter.h"
#include "TlSerializeData.h"
#include "TlLogging.h"

/// 入力データの読み込み処理を行うクラス
class PdfUserInput {
public:
    PdfUserInput(const std::string& filePath = "fl_Userinput");
    ~PdfUserInput();

public:
    TlSerializeData getSerializeData() const;
    TlParameter getParameter() const;
    
    /// aliasの変換を行う
    void alias();
    bool check();

    void load();

public:
    TlSerializeData getBasisInfo(const std::string& basisName);
    
private:
    void load_conventional();

    void molecule_geometry_cartesian_input(const std::string& str);
    void moleculeBasisSetOrbital(const std::string& str);
    void moleculeBasisSetDensityAuxiliary(const std::string& str);
    void moleculeBasisSetExchangeAuxiliary(const std::string& str);
    void moleculeBasisSetGridFree(const std::string& str);

private:
    std::string filePath_;
    TlSerializeData data_;
    TlParameter param_;

    TlLogging& log_;
};

#endif // PDFUSERINPUT_H

