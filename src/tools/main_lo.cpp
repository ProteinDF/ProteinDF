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
#include <string>

#include "DfLocalize.h"
#include "TlGetopt.h"
#include "TlLogging.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("%s [OPTIONS] basisset_name ...", progname.c_str()) << std::endl;
    std::cout << std::endl;
    std::cout << " localize C matrix" << std::endl;
    std::cout << " -p PATH       set ProteinDF parameter file. default = "
                 "pdfparam.mpac"
              << std::endl;
    std::cout << " -c PATH       set C matrix path" << std::endl;
    std::cout << " -r            restart lo calculation" << std::endl;
    std::cout << " -g <brd_file> specify group list" << std::endl;
    std::cout << " -m <MODE>     switch pairing algorithm (ordered, big-big(defalt), big-small, small-small)" << std::endl;
    std::cout << " -h            show help" << std::endl;
    // std::cout << " -v            verbose output" << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "c:hm:g:p:rv");
    const bool isShowHelp = (opt["h"] == "defined");
    // const bool isVerbose = (opt["v"] == "defined");

    if (isShowHelp) {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }

    TlLogging& log = TlLogging::getInstance();
    log.setFilePath("lo-output.txt");

    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }

    std::string groupFilePath = "";
    if (opt["g"].empty() != true) {
        groupFilePath = opt["g"];
    }

    std::string inputCMatrixPath = "";
    if (opt["c"].empty() != true) {
        inputCMatrixPath = opt["c"];
    }

    DfLocalize::PairingOrder pairingOrder = DfLocalize::PO_BIG_BIG;
    if (opt["m"].empty() != true) {
        const std::string sInputMode = TlUtils::toUpper(opt["m"]);
        if (sInputMode == "BIG-BIG") {
            pairingOrder = DfLocalize::PO_BIG_BIG;
        } else if (sInputMode == "BIG-SMALL") {
            pairingOrder = DfLocalize::PO_BIG_SMALL;
        } else if (sInputMode == "SMALL-SMALL") {
            pairingOrder = DfLocalize::PO_SMALL_SMALL;
        } else if (sInputMode == "ORDERED") {
            pairingOrder = DfLocalize::PO_ORDERED;
        } else {
            std::cerr << "illegal parameter: -m " << opt["m"] << std::endl;
        }
    }

    const bool isRestart = (opt["r"] == "defined");

    // if (isVerbose) {
    // std::cerr << "parameter path: " << pdfParamPath << std::endl;
    // std::cerr << "iteration: " << iteration << std::endl;
    // }

    // setup
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    std::cout << "group file: " << groupFilePath << std::endl;
    TlSerializeData groupData;
    if (groupFilePath.empty() != true) {
        log.info(TlUtils::format("load group file: %s", groupFilePath.c_str()));

        TlMsgPack mpac;
        mpac.load(groupFilePath);

        groupData = mpac.getSerializeData();
        // std::cout << groupData.str() << std::endl;
    }

    // lo
    DfLocalize lo(&param);
    if (!inputCMatrixPath.empty()) {
        lo.setCMatrixPath(inputCMatrixPath);
    }
    lo.setRestart(isRestart);
    lo.setGroup(groupData);
    lo.setPairingOrder(pairingOrder);

    lo.exec();

    {
        TlMsgPack mpac(param);
        mpac.save(pdfParamPath);
    }

    return EXIT_SUCCESS;
}
