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

#include "TlCommunicate.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlFile.h"
#include "TlMatrixObject.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

typedef TlMatrixObject::index_type index_type;

void showHelp()
{
    std::cout << "divideGlobalMatrix file path" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -h: show help(this)" << std::endl;
    std::cout << "  -v: verbose" << std::endl;
}


int main(int argc, char* argv[])
{
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlGetopt opt(argc, argv, "b:hv");

    // parameters
    bool isVerbose = false;
    bool isHelpMode = false;
    int blockSize = 64;
    int numOfArgs = 1;
    std::string readMatrixPath = "";
    std::string outputMatrixPath = "";
    
    if (rComm.isMaster() == true) {
        isVerbose = (opt["v"] == "defined");
        numOfArgs = opt.getCount();
        isHelpMode = (opt["h"] == "defined");
        if (opt["b"].empty() == false) {
            blockSize = std::atoi(opt["b"].c_str());
        }
    
        readMatrixPath = opt[1];
        outputMatrixPath = opt[2];

        if (isVerbose == true) {
            std::cerr << TlUtils::format("load matrix: %s.", readMatrixPath.c_str()) << std::endl;
            std::cerr << TlUtils::format("save matrix path: %s.", outputMatrixPath.c_str()) << std::endl;
            std::cerr << TlUtils::format("block size = %d", blockSize) << std::endl;
        }
    }

    // setup
    rComm.broadcast(isVerbose);
    rComm.broadcast(isHelpMode);
    rComm.broadcast(blockSize);
    rComm.broadcast(numOfArgs);
    rComm.broadcast(readMatrixPath);
    rComm.broadcast(outputMatrixPath);
    if ((numOfArgs < 3) || (isHelpMode == true))  {
        if (rComm.isMaster() == true) {
            showHelp();
        }
        rComm.finalize();
        return EXIT_FAILURE;
    }
    TlDistributeMatrix::setSystemBlockSize(blockSize);
    
    if (TlDistributeSymmetricMatrix::isLoadable(readMatrixPath) == true) {
        if (isVerbose == true) {
            if (rComm.isMaster() == true) {
                std::cerr << "using symmetric matrix loader." << std::endl;
            }
        }
        // 対称行列
        TlDistributeSymmetricMatrix M;
        TlDistributeMatrix::setUsingPartialIO(false);
        M.load(readMatrixPath);
 
        TlDistributeMatrix::setUsingPartialIO(true);
        M.save(outputMatrixPath + TlUtils::format(".%d", rComm.getRank()));
    } else if (TlDistributeMatrix::isLoadable(readMatrixPath) == true) {
        if (isVerbose == true) {
            if (rComm.isMaster() == true) {
                std::cerr << "using normal matrix loader." << std::endl;
            }
        }
        TlDistributeMatrix M;
        TlDistributeMatrix::setUsingPartialIO(false);
        M.load(readMatrixPath);

        TlDistributeMatrix::setUsingPartialIO(true);
        M.save(outputMatrixPath + TlUtils::format(".%d", rComm.getRank()));
    } else {
        std::cerr << "unknown file type: " << readMatrixPath << std::endl;
        rComm.finalize();
        return EXIT_FAILURE;
    }
        

    rComm.finalize();
    return EXIT_SUCCESS;
}


