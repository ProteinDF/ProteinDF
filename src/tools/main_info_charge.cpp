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

#include <cstdlib>
#include <iostream>
#include "DfInfo.h"
#include "Fl_Geometry.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hf:");

    std::string paramPath = "pdfparam.mpac";
    if (opt["f"].empty() == false) {
        paramPath = opt["f"];
    }

    // パラメータファイルの読み込み
    TlMsgPack mpac;
    mpac.load(paramPath);
    TlSerializeData param = mpac.getSerializeData();

    Fl_Geometry flGeom(param["coordinates"]);
    DfInfo dfInfo(&param);

    const int nuclCharge = flGeom.getTotalCharge();
    const int elecCharge = -dfInfo.getNumOfElectrons();

    std::cout << "nucleus charge            = " << nuclCharge << std::endl;
    std::cout << "nucleus charge(without X) = "
              << flGeom.getTotalChargeWithoutX() << std::endl;
    std::cout << "electron charge           = " << elecCharge << std::endl;
    std::cout << "charge                    = " << nuclCharge + elecCharge
              << std::endl;
    std::cout << "formula: " << flGeom.getFormula() << std::endl;

    return EXIT_SUCCESS;
}
