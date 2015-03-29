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
#include "TlMatrixObject.h"
#include "TlRowVectorMatrix.h"
#include "TlUtils.h"
#include "TlGetopt.h"

typedef TlMatrixObject::index_type index_type;

void showHelp(const std::string& progname)
{
    std::cout << TlUtils::format("%s [options] INPUT_BASENAME OUTPUT_BASENAME NUM_SUBUNITS", progname.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -m:      use memory manager" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hvm");
    
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    const bool bVerbose = (opt["v"] == "defined");
    const bool isUsingMemManager = (opt["m"] == "defined");

    if (opt.getCount() <= 3) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }
    std::string inputBaseName = opt[1];
    std::string outputBaseName = opt[2];
    const int output_numOfSubunits = std::atoi(opt[3].c_str());
    if (output_numOfSubunits < 1) {
        std::cerr << TlUtils::format("output # of subunits is too small: %d", output_numOfSubunits) << std::endl;
        return EXIT_FAILURE;
    }

    // check file
    index_type numOfVectors = 0;
    index_type sizeOfVector = 0;
    int input_numOfSubunits = 1;
    {
        if (bVerbose == true) {
            std::cerr << "check input matrices..." << std::endl;
        }

        int subunitID = 0;
        const std::string inputPath0 = TlVectorMatrixObject::getFileName(inputBaseName, subunitID);
        bool isLoadable = TlVectorMatrixObject::isLoadable(inputPath0,
                                                           &numOfVectors, &sizeOfVector,
                                                           &input_numOfSubunits, &subunitID);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath0 << std::endl;
            return EXIT_FAILURE;
        }

        if (bVerbose == true) {
            std::cerr << TlUtils::format("#vectors: %d", numOfVectors) << std::endl;
            std::cerr << TlUtils::format("#size: %d", sizeOfVector) << std::endl;
            std::cerr << TlUtils::format("#subunits: %d", input_numOfSubunits) << std::endl;
        }
    }

    // load & set
    std::vector<TlRowVectorMatrix*> outputs(output_numOfSubunits);
    for (int i = 0; i < output_numOfSubunits; ++i) {
        outputs[i] = new TlRowVectorMatrix(numOfVectors, sizeOfVector,
                                           output_numOfSubunits, i,
                                           isUsingMemManager);
    }

    for (int i = 0; i < input_numOfSubunits; ++i) {
        TlRowVectorMatrix vm;
        vm.load(inputBaseName, i);
        
        for (index_type j = 0; j < numOfVectors; ++j) {
            if (vm.getSubunitID(j) == i) {
                const TlVector v = vm.getVector(j);
                const int output_subunitID = outputs[0]->getSubunitID(j);
                outputs[output_subunitID]->setVector(j, v);
            }
        }
    }

    // save
    for (int i = 0; i < output_numOfSubunits; ++i) {
        outputs[i]->save(outputBaseName);
        
        delete outputs[i];
        outputs[i] = NULL;
    }
    outputs.clear();

    return EXIT_SUCCESS;
}

