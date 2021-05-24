#include <iostream>
#include <string>

#include "TlFile.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"

typedef TlMatrixObject::index_type index_type;

void showHelp(const std::string& prog) {
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hv");

    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    const bool verbose = (opt["v"] == "defined");

    if (opt.getCount() <= 2) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }
    std::string inputPath = opt[1];
    std::string outputPath = opt[2];
    if (verbose) {
        std::cerr << "input base path: " << inputPath << std::endl;
        std::cerr << "output: " << outputPath << std::endl;
    }

    // check
    // {
    //     index_type numOfRows = 0;
    //     index_type numOfCols = 0;
    //     int numOfSubunits = 0;
    //     int sizeOfChunk = 0;
    //     {
    //         int subunitID = 0;
    //         const std::string inputPath0 = TlDenseMatrix_arrays_Object::getFileName(inputBasePath, subunitID);
    //         index_type sizeOfVector, numOfVectors;
    //         const bool isLoadable = TlDenseMatrix_arrays_Object::isLoadable(inputPath0, &numOfVectors, &sizeOfVector,
    //                                                                         &numOfSubunits, &subunitID,
    //                                                                         &sizeOfChunk);
    //         if (isLoadable != true) {
    //             std::cerr << "can not open file: " << inputPath0 << std::endl;
    //             return EXIT_FAILURE;
    //         }
    //
    //         // row-vector matrix
    //         numOfRows = numOfVectors;
    //         numOfCols = sizeOfVector;
    //
    //         if (verbose) {
    //             std::cerr << "rows: " << numOfRows << std::endl;
    //             std::cerr << "cols: " << numOfCols << std::endl;
    //             std::cerr << "units: " << numOfSubunits << std::endl;
    //             std::cerr << "size of chunk: " << sizeOfChunk << std::endl;
    //         }
    //     }
    // }

    // load
    if (verbose) {
        std::cerr << TlUtils::format("loading %s ...", inputPath.c_str()) << std::endl;
    }
    TlDenseGeneralMatrix_arrays_RowOriented m;
    {
        index_type numOfRows = 0;
        index_type numOfCols = 0;
        int numOfSubunits = 0;
        int sizeOfChunk = 0;
        index_type sizeOfVector, numOfVectors;
        int subunitId;

        const bool isLoadable = TlDenseMatrix_arrays_Object::isLoadable(inputPath, &numOfVectors, &sizeOfVector,
                                                                        &numOfSubunits, &subunitId, &sizeOfChunk);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath << std::endl;
            return EXIT_FAILURE;
        }

        m.loadSubunitFile(inputPath);
        std::cerr << TlUtils::format("row: %d", m.getNumOfRows()) << std::endl;
        std::cerr << TlUtils::format("col: %d", m.getNumOfCols()) << std::endl;
        std::cerr << TlUtils::format("sizeOfChunk: %d", m.getSizeOfChunk()) << std::endl;
        std::cerr << TlUtils::format("%d/%d", m.getSubunitID(), m.getNumOfSubunits()) << std::endl;
    }

    // save
    if (verbose) {
        std::cerr << TlUtils::format("saving %s ...", outputPath.c_str()) << std::endl;
    }
    m.saveSubunitFileWithReservedSizeOfVector(outputPath);

    return EXIT_SUCCESS;
}
