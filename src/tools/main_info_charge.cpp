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

#include <iostream>
#include <cstdlib>
#include <map>
#include <string>

#include "Fl_Geometry.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hf:");

    std::string paramPath = "pdfparam.mpac";
    if (opt["f"].empty() == false) {
        paramPath = opt["f"];
    }
    
    // パラメータファイルの読み込み
    TlMsgPack mpac;
    mpac.load(paramPath);
    TlSerializeData param = mpac.getSerializeData();

    //Fl_Geometry flGeom(Fl_Geometry::getDefaultFileName());
    Fl_Geometry flGeom(param["coordinates"]);

    std::map<std::string, int> component;
    double charge = 0.0;
    double chargeWithoutX = 0.0;
    const int numOfAtoms = flGeom.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string symbol = flGeom.getAtomSymbol(i);
        ++(component[symbol]);
        const double currentCharge = flGeom.getCharge(i);

        charge += currentCharge;
        if (symbol != "X") {
            chargeWithoutX += currentCharge;
        }
    }

    std::cout << "charge            = " << charge << std::endl;
    std::cout << "charge(without X) = " << chargeWithoutX << std::endl;

    for (std::map<std::string, int>::const_iterator p = component.begin();
         p != component.end(); ++p) {
        std::cout << p->first << p->second;
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}


