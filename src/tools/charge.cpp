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
    Fl_Geometry flGeom(param["model"]["coordinates"]);

    std::map<std::string, int> component;
    double charge = 0.0;
    double chargeWithoutX = 0.0;
    const int numOfAtoms = flGeom.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string symbol = flGeom.getAtom(i);
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


