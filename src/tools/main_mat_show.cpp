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

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ghv");

    const bool bVerbose = (opt["v"] == "defined");
    const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << sPath << std::endl;
    }

    if (TlSymmetricMatrix::isLoadable(sPath) == true) {
        TlSymmetricMatrix M;
        M.load(sPath);

        if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            M.print(std::cout);
        }
    } else if (TlMatrix::isLoadable(sPath) == true) {
        TlMatrix M;
        M.load(sPath);

        if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            M.print(std::cout);
        }
    } else {
        std::cerr << "unknown file type: " << sPath << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


