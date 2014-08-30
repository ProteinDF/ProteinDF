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

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"


struct MatrixElement {
public:
    explicit MatrixElement(int r=0, int c=0, double v=0.0) : row(r), col(c), value(v) {
    }
    
public:
    int row;
    int col;
    double value;
};

struct MatrixElement_cmp {
public:
    bool operator()(const MatrixElement& a, const MatrixElement& b) {
        return (std::fabs(a.value) > std::fabs(b.value));
    }
};


void check(const MatrixElement& me, std::vector<MatrixElement>* pRanking)
{
    std::vector<MatrixElement>::iterator it = std::lower_bound(pRanking->begin(),
                                                               pRanking->end(),
                                                               me,
                                                               MatrixElement_cmp());
    if (it != pRanking->end()) {
        pRanking->insert(it, me);
        pRanking->pop_back();
    }
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ht:v");

    int numOfRanks = 10;
    if (opt["t"].empty() == false) {
        numOfRanks = std::atoi(opt["t"].c_str());
    }
    bool bVerbose = (opt["v"] == "defined");

    const std::string path = opt[1];
    if (bVerbose) {
        std::cerr << "ranking top " << numOfRanks << std::endl;
        std::cerr << "loading... " << path << std::endl;
    }

    int errorCode = 0;

    std::ifstream ifs;
    ifs.open(path.c_str());
    if (ifs.fail()) {
        std::cerr << "could not open file. " << path << std::endl;
        return 1;
    }

    std::vector<MatrixElement> ranking(numOfRanks);
    if (TlSymmetricMatrix::isLoadable(ifs) == true) {
        TlSymmetricMatrix matrix;
        matrix.load(path);

        const int dim = matrix.getNumOfRows();
        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c <= r; ++c) {
                check(MatrixElement(r, c, matrix.get(r, c)), &ranking);
            }
        }
    } else if (TlMatrix::isLoadable(ifs) == true) {
        TlMatrix matrix;
        matrix.load(path);

        const int row = matrix.getNumOfRows();
        const int col = matrix.getNumOfCols();
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                check(MatrixElement(r, c, matrix.get(r, c)), &ranking);
            }
        }
    } else {
        std::cerr << "could not open: " << path << std::endl;
        errorCode = 1;
    }

    // display
    for (int i = 0; i < numOfRanks; ++i) {
        const MatrixElement me = ranking[i];
        std::cout << TlUtils::format("No.%2d (%6d, %6d) = % 8.3f",
                                     i, me.row, me.col, me.value)
                  << std::endl;
    }
    
    
    return errorCode;
}


