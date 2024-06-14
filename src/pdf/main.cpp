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

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>

#include "ProteinDF.h"
#include "TlGetopt.h"
#include "TlLogging.h"
#include "TlSerializeData.h"

#ifdef HAVE_VIENNACL
#include "tl_viennacl.h"
#endif  // HAVE_VIENNACL

//#ifdef __FUJITSU
//#define PDF_MAIN MAIN__
//#else
#define PDF_MAIN main
//#endif  // __FUJITSU

int PDF_MAIN(int argc, char* argv[]) {
    // setup parameters
    TlGetopt opt(argc, argv, "a:dro:");

    TlLogging& log = TlLogging::getInstance();
    std::string output = "fl_Out_Std";
    if (opt["o"].empty() != true) {
        output = opt["o"];
    }
    log.setFilePath(output);

    if (opt["d"] == "defined") {
        log.setLevel(TlLogging::TL_DEBUG);
    }

    // ViennaCL
#ifdef HAVE_VIENNACL
    {
        int deviceId = 0;
        if (!opt["a"].empty()) {
            deviceId = std::atoi(opt["a"].c_str());
        }

        TlViennaCL vcl;
        vcl.setupAllAvailableDevices();
        log.info(vcl.listDevices());

        vcl.switchDevice(deviceId);
        log.info(vcl.listCurrentDevice());
    }
#endif  // HAVE_VIENNACL

    // command option
    bool isRestart = false;
    if (opt["r"] == "defined") {
        isRestart = true;
    }

    // do ProteinDF
    ProteinDF PDF;
    if (isRestart == true) {
        std::cout << "restart calculation ..." << std::endl;
        PDF.restart("pdfparam.mpac");
    } else {
        std::cout << "normal start calculation ..." << std::endl;
        PDF.run();
    }

    return EXIT_SUCCESS;
}
