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

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_matrix_utils.h"

void help(const std::string& progname) {
    std::cout << TlUtils::format("Usage: %s [options]... FILE1 FILE2",
                                 progname.c_str())
              << std::endl;
    std::cout << "compare ProteinDF matrix files" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s FILE      save difference matrix" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hs:v");

    bool bVerbose = (opt["v"] == "defined");
    if (opt.getCount() < 2) {
        help(opt[0]);
        std::exit(1);
    }

    const std::string sPath1 = opt[1];
    const std::string sPath2 = opt[2];
    if (bVerbose) {
        std::cerr << "loading... " << sPath1 << std::endl;
        std::cerr << "loading... " << sPath2 << std::endl;
    }

    std::string savePath = "";
    if (opt["s"].empty() == false) {
        savePath = opt["s"];
    }

    int errorCode = 0;

    std::ifstream ifs1;
    ifs1.open(sPath1.c_str());
    if (ifs1.fail()) {
        std::cerr << "could not open file. " << sPath1 << std::endl;
        return 1;
    }

    std::ifstream ifs2;
    ifs2.open(sPath2.c_str());
    if (ifs2.fail()) {
        std::cerr << "could not open file. " << sPath2 << std::endl;
        return 1;
    }

    if (TlMatrixUtils::isLoadable(sPath1, TlMatrixObject::RLHD) == true) {
        if (TlMatrixUtils::isLoadable(sPath2, TlMatrixObject::RLHD) == true) {
            TlDenseSymmetricMatrix_Lapack m1, m2;
            m1.load(sPath1);
            m2.load(sPath2);

            m1 -= m2;
            if (savePath.empty() == false) {
                m1.save(savePath);
            } else {
                std::cout << m1 << std::endl;
            }
        } else {
            std::cerr << "could not open: " << sPath2 << std::endl;
            errorCode = 1;
        }
    } else if (TlMatrixUtils::isLoadable(sPath1, TlMatrixObject::CSFD) ==
               true) {
        if (TlMatrixUtils::isLoadable(sPath2, TlMatrixObject::CSFD) == true) {
            TlDenseGeneralMatrix_Lapack m1, m2;
            m1.load(sPath1);
            m2.load(sPath2);

            m1 -= m2;
            if (savePath.empty() == false) {
                m1.save(savePath);
            } else {
                std::cout << m1 << std::endl;
            }
        } else {
            std::cerr << "could not open: " << sPath2 << std::endl;
            errorCode = 1;
        }
    } else {
        std::cerr << "could not open: " << sPath1 << std::endl;
        errorCode = 1;
    }

    return errorCode;
}
