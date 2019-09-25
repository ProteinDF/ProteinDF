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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "PdfKeyword.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

void showHelp(const std::string& progname) {
    std::cout << "output ProteinDF keyword" << std::endl;
    std::cout << TlUtils::format("USAGE: %s [OPTIONS]", progname.c_str())
              << std::endl;
    std::cout << std::endl;
    std::cout << "  -a:      output hidden parameters" << std::endl;
    std::cout << "  -c:      output CSV format" << std::endl;
    std::cout << "  -r:      output ReST format" << std::endl;
    std::cout << "  -j:      output Japanese" << std::endl;
    std::cout << "  -m FILE: output MsgPack file" << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hacjm:r");

    if (opt["h"] == "defined") {
        showHelp(opt[0]);
    }

    const bool showAll = (opt["a"] == "defined");
    const bool isCSV = (opt["c"] == "defined");
    const bool isJP = (opt["j"] == "defined");
    const bool is_reST = (opt["r"] == "defined");
    const std::string mpacFilePath = opt["m"];

    PdfKeyword kwd;

    if (isCSV == true) {
        std::string csv;
        if (isJP == true) {
            csv = kwd.getCSV_jp(showAll);
        } else {
            csv = kwd.getCSV(showAll);
        }
        std::cout << csv << std::endl;
    } else if (is_reST == true) {
        std::string reST;
        reST = kwd.get_reST_jp(showAll);
        std::cout << reST << std::endl;
    } else if (mpacFilePath.empty() != true) {
        const TlSerializeData data = kwd.getSerializeData();
        TlMsgPack mpac(data);
        mpac.save(mpacFilePath);
    } else {
        kwd.printDefault(std::cout);
    }

    return EXIT_SUCCESS;
}
