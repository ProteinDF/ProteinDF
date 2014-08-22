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

//   const bool bVerbose = (opt["v"] == "defined");
//   const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath1 = opt[1];
    std::string sPath2 = opt[2];
    std::string sPath3 = opt[3];
//   if (bVerbose){
//     std::cerr << "loading... " << sPath << std::endl;
//     std::cerr << "loading... " << sPath << std::endl;
//   }

    if (TlSymmetricMatrix::isLoadable(sPath1) != true) {
        std::cerr << "can not open file: " << sPath1 << std::endl;
        return EXIT_FAILURE;
    }
    if (TlSymmetricMatrix::isLoadable(sPath2) != true) {
        std::cerr << "can not open file: " << sPath1 << std::endl;
        return EXIT_FAILURE;
    }

    TlSymmetricMatrix M1;
    M1.load(sPath1);
    const std::size_t nRows1 = M1.getNumOfRows();
    const std::size_t nCols1 = M1.getNumOfCols();

    TlSymmetricMatrix M2;
    M2.load(sPath2);
    const std::size_t nRows2 = M2.getNumOfRows();
    const std::size_t nCols2 = M2.getNumOfCols();

    const std::size_t dim = nRows1 + nRows2;
    assert(dim == nCols1 + nCols2);
    M1.resize(dim);
    for (std::size_t r = 0; r < nRows2; ++r) {
        const std::size_t row = nRows1 + r;
        for (std::size_t c = 0; c < nCols2; ++c) {
            const std::size_t col = nCols1 + c;
            M1(row, col) = M2(r, c);
        }
    }

    M1.save(sPath3);

    return EXIT_SUCCESS;
}


