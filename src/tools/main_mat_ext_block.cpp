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

void showHelp(const std::string& name)
{
    std::cout << TlUtils::format("%s [options] base_matrix_path reference_matrix_path  output_path", name.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "t:b:l:r:hv");
    
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }

    const bool isVerbose = (opt["v"] == "defined");
    if (opt.getCount() <= 2) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }
    std::string baseMatrixPath = opt[1];
    std::string refMatrixPath = opt[2];
    std::string outputMatrixPath = opt[3];

    TlSymmetricMatrix mat1;
    TlMatrix::index_type dim1 = 0;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << baseMatrixPath << std::endl;
    }
    if (TlSymmetricMatrix::isLoadable(baseMatrixPath)) {
        mat1.load(baseMatrixPath);
        dim1 = mat1.getNumOfRows();
    } else {
        std::cerr << TlUtils::format("cannot load: %s", baseMatrixPath.c_str()) << std::endl;
        std::cerr << "create new matrix." << std::endl;
    }

    TlSymmetricMatrix mat2;
    TlMatrix::index_type dim2 = 0;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << refMatrixPath << std::endl;
    }
    if (TlSymmetricMatrix::isLoadable(refMatrixPath)) {
        mat2.load(refMatrixPath);
        dim2 = mat2.getNumOfRows();
    } else {
        std::cerr << TlUtils::format("cannot load: %s", refMatrixPath.c_str()) << std::endl;
        return EXIT_FAILURE;
    }

    const TlMatrix::index_type dim3 = dim1 + dim2;

    mat1.resize(dim3);
    for (TlMatrix::index_type r = 0; r < dim2; ++r) {
        for (TlMatrix::index_type c = 0; c <= r; ++c) {
            mat1.set(dim1 + r, dim1 + c, mat2.get(r, c));
        }
    }

    if (isVerbose == true) {
        std::cerr << "save matrix: " << outputMatrixPath << std::endl;
    }
    mat1.save(outputMatrixPath);

    return EXIT_SUCCESS;
}


