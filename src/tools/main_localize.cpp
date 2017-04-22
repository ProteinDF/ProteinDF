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
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlGetopt.h"


void showHelp(const std::string& progname)
{
    std::cout << TlUtils::format("%s [OPTIONS] basisset_name ...", progname.c_str()) << std::endl;
    std::cout << std::endl;
    std::cout << " localize C matrix" << std::endl;
    std::cout << " -p PATH       set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -c PATH       set C matrix path" << std::endl;
    std::cout << " -h            show help" << std::endl;
    std::cout << " -v            verbose output" << std::endl;
}


int main(int argc, char *argv[])
{
    TlGetopt opt(argc, argv, "c:hp:v");
    const bool isVerbose = (opt["v"] == "defined");
    const bool isShowHelp = (opt["h"] == "defined");

    if (isShowHelp) {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }

    std::string inputCMatrixPath = "";
    if (opt["c"].empty() != true) {
        inputCMatrixPath = opt["c"];
    }
    
    if (isVerbose) {
        std::cerr << "parameter path: " << pdfParamPath << std::endl;
    }

    
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }
    
    DfLocalize lo(&param);
    lo.localize(inputCMatrixPath);

    {
        TlMsgPack mpac(param);
        mpac.save(pdfParamPath);
    }

    return 0;
}




