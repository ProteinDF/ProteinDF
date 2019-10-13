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
#include <string>
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlRotateLCAO.h"

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("USAGE: %s <TARGET_PDFPARAM>",
                                 progname.c_str())
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout
        << "  -i <PATH>: define input matrix path. default: rotate_lcao.input"
        << std::endl;
    std::cout
        << "  -o <PATH>: define output matrix path. default: rotate_lcao.output"
        << std::endl;
    std::cout << "  -m <PATH>: define rotate matrix." << std::endl;
    std::cout << "  -h: show help(this)" << std::endl;
    std::cout << "  -v: verbose" << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hi:o:m:v");
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    const bool verbose = (opt["v"] == "defined");

    // const std::string refPdfParamPath = opt[1];
    const std::string newPdfParamPath = opt[1];
    std::string inLcaoPath = "rotate_lcao.input";
    if (opt["i"] != "") {
        inLcaoPath = opt["i"];
    }
    std::string outLcaoPath = "rotate_lcao.output";
    if (opt["o"] != "") {
        outLcaoPath = opt["o"];
    }
    std::string rotateMatrixPath = "rotate.matrix";
    if (opt["m"] != "") {
        rotateMatrixPath = opt["m"];
    }

    //     TlMsgPack inMsgPack;
    //     if (verbose == true) {
    //         std::cerr << "reference parameter: " << refPdfParamPath <<
    //         std::endl;
    //     }
    //     inMsgPack.load(refPdfParamPath);
    //     const TlSerializeData inParam = inMsgPack.getSerializeData();

    TlMsgPack outMsgPack;
    if (verbose == true) {
        std::cerr << "target parameter: " << newPdfParamPath << std::endl;
    }
    outMsgPack.load(newPdfParamPath);
    const TlSerializeData outParam = outMsgPack.getSerializeData();

    TlMatrix lcao;
    lcao.load(inLcaoPath);

    TlMatrix rot;
    rot.load(rotateMatrixPath);

    const TlOrbitalInfo orbInfo(outParam["coordinates"], outParam["basis_set"]);
    // orbInfo.printCGTOs(std::cerr);

    TlRotateLCAO rotateLCAO(orbInfo);
    const TlMatrix newLCAO = rotateLCAO.exec(lcao, rot);

    newLCAO.save(outLcaoPath);
}
