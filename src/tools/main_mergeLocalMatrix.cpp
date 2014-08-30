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
#include <vector>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlFile.h"
#include "TlMatrixObject.h"
#include "TlDistributeMatrix.h"
#include "TlFileMatrix.h"
#include "TlFileSymmetricMatrix.h"

typedef TlMatrixObject::index_type index_type;

void showHelp(const std::string& progname)
{
    std::cout << TlUtils::format("USAGE: %s file... path", progname.c_str()) << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -h: show help(this)" << std::endl;
    std::cout << "  -v: verbose" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");
    const int numOfArgs = opt.getCount();
    if ((numOfArgs < 3) || (opt["h"] == "defined"))  {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }

    std::vector<std::string> readMatrixPaths(numOfArgs -2);
    for (int i = 0; i < (numOfArgs -2); ++i) {
        readMatrixPaths[i] = opt[i +1];
    }
    const std::string outputMatrixPath = opt[numOfArgs -1];


    std::vector<std::string>::const_iterator itEnd = readMatrixPaths.end();
    for (std::vector<std::string>::const_iterator it = readMatrixPaths.begin(); it != itEnd; ++it) {
        const std::string& readFilePath = *it;
        if (isVerbose == true) {
            std::cerr << TlUtils::format("reading %s.", readFilePath.c_str())
                      << std::endl;
        }

        if (TlFile::isExist(readFilePath) == false) {
            std::cerr << TlUtils::format("not found: %s.", readFilePath.c_str()) << std::endl;
            return EXIT_FAILURE;
        }

        std::ifstream ifs;
        ifs.open(readFilePath.c_str());
        TlDistributeMatrix::LocalMatrixHeader lmh;
        lmh.load(&ifs);

        // check symmetric matrix or not
        bool isSymmetricMatrix = false;
        switch (lmh.type) {
        case 16:
            // normal matrix
            isSymmetricMatrix = false;
            break;

        case 18:
            // symmetric matrix
            isSymmetricMatrix = true;
            assert(lmh.globalRow == lmh.globalCol);
            break;

        default:
            // unknown
            std::cerr << TlUtils::format("unknown file type: %s.", readFilePath.c_str()) << std::endl;
            return EXIT_FAILURE;
            break;
        }

        // load local matrix elements
        std::size_t numOfElements = lmh.myRows * lmh.myCols; 
        double* pValues = new double[numOfElements];
        ifs.read((char*)pValues, sizeof(double) * numOfElements);
        
        // create output matrix
        TlFileMatrix* pOutputMatrix = NULL;
        if (TlFile::isExist(outputMatrixPath) == false) {
            if (isVerbose == true) {
                std::cerr << TlUtils::format("create output matrix: %s.", outputMatrixPath.c_str())
                          << std::endl;
            }

            if (isSymmetricMatrix == true) {
                pOutputMatrix = new TlFileSymmetricMatrix(outputMatrixPath,
                                                          lmh.globalRow);
            } else {
                pOutputMatrix = new TlFileMatrix(outputMatrixPath,
                                                 lmh.globalRow, lmh.globalCol);
            }
        } else {
            if (isVerbose == true) {
                std::cerr << TlUtils::format("loading output matrix: %s.", outputMatrixPath.c_str())
                          << std::endl;
            }
            
            if (isSymmetricMatrix == true) {
                pOutputMatrix = new TlFileSymmetricMatrix(outputMatrixPath);
            } else {
                pOutputMatrix = new TlFileMatrix(outputMatrixPath);
            }
        }

        // output matrix check
        if ((pOutputMatrix->getNumOfRows() != lmh.globalRow) ||
            (pOutputMatrix->getNumOfCols() != lmh.globalCol)) {
            std::cerr << "ERROR: output matrix size is not consistent with input matrix."
                      << std::endl;
            return EXIT_FAILURE;
        }
        
        // make index table(row)
        std::vector<index_type> rowIndexTable;
        {
            const int blockIndex = lmh.myProcRow * lmh.blockSize;
            const int incrementBlockIndex = lmh.procGridRow * lmh.blockSize;
            rowIndexTable.reserve(lmh.myRows);
            for (index_type r = 0; r < lmh.myRows; ++r) {
                const div_t d = std::div(r, lmh.blockSize);
                const int i = blockIndex + (incrementBlockIndex * d.quot) + d.rem;
                if (i < lmh.globalRow) {
                    rowIndexTable.push_back(i);
                } else {
                    break;
                }
            }
            std::vector<index_type>(rowIndexTable).swap(rowIndexTable);
        }

        // make index table(col)
        std::vector<index_type> colIndexTable;
        {
            const int blockIndex = lmh.myProcCol * lmh.blockSize;
            const int incrementBlockIndex = lmh.procGridCol * lmh.blockSize;
            colIndexTable.reserve(lmh.myCols);
            for (index_type c = 0; c < lmh.myCols; ++c) {
                const div_t d = std::div(c, lmh.blockSize);
                const int i = blockIndex + (incrementBlockIndex * d.quot) + d.rem;
                if (i < lmh.globalCol) {
                    colIndexTable.push_back(i);
                } else {
                    break;
                }
            }
            std::vector<index_type>(colIndexTable).swap(colIndexTable);
        }
        
        // copy elements
        const int rowIndexTableSize = rowIndexTable.size();
        const int colIndexTableSize = colIndexTable.size();
        if (isSymmetricMatrix == false) {
            std::size_t index = 0;
            for (int c = 0; c < colIndexTableSize; ++c) {
                const index_type col = colIndexTable[c];
                for (int r = 0; r < rowIndexTableSize; ++r) {
                    const index_type row = rowIndexTable[r];

                    pOutputMatrix->set(row, col, pValues[index]);
                    ++index;
                }
            }
        } else {
            std::size_t index = 0;
            for (int c = 0; c < colIndexTableSize; ++c) {
                const index_type col = colIndexTable[c];
                for (int r = 0; r < rowIndexTableSize; ++r) {
                    const index_type row = rowIndexTable[r];
                    if (row >= col) {
                        pOutputMatrix->set(row, col, pValues[index]);
                    }
                    ++index;
                }
            }
        }
        
        // clean up
        ifs.close();
        delete[] pValues;
        pValues = NULL;
        delete pOutputMatrix;
        pOutputMatrix = NULL;
    }



    return EXIT_SUCCESS;
}


