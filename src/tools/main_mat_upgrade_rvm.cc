#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_mmap_roworiented.h"
#include "tl_matrix_utils.h"

void showHelp(const std::string& name) {
    std::cout << TlUtils::format("%s [options] INPUT_RVM_FILE OUTPUT_RVM_FILE", name.c_str())
              << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose output" << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hv");

    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }

    const bool verbose = (opt["v"] == "defined");

    std::string inputPath = opt[1];
    std::string outputPath = opt[2];

    TlMatrixObject::HeaderInfo headerInfo;
    const bool isLoadable = TlMatrixUtils::getHeaderInfo(inputPath, &headerInfo);

    if (isLoadable == false) {
        std::cerr << "cannot load file: " << inputPath << std::endl;
        return EXIT_FAILURE;
    }
    if (headerInfo.matrixType != TlMatrixObject::ABGD) {
        std::cerr << "matrix type mismatch. matrix type: " << TlMatrixObject::matrixTypeStr(headerInfo.matrixType) << std::endl;
        return EXIT_FAILURE;
    }
    if (headerInfo.version == 1) {
        std::cerr << "matrix format version: " << headerInfo.version << std::endl;
        std::cerr << "No need to upgrade format." << std::endl;
        return EXIT_SUCCESS;
    }

    if (verbose) {
        std::cerr << "load matrix: " << inputPath << std::endl;
    }
    const TlDenseGeneralMatrix_arrays_mmap_RowOriented inMat(inputPath);
    const TlMatrixObject::index_type numOfRows = inMat.getNumOfRows();
    const TlMatrixObject::index_type numOfCols = inMat.getNumOfCols();
    const int sizeOfChunk = inMat.getSizeOfChunk();
    const int numOfSubunits = inMat.getNumOfSubunits();
    const int subunitId = inMat.getSubunitID();

    if (verbose) {
        std::cerr << "output matrix: " << outputPath << std::endl;
    }
    TlDenseGeneralMatrix_arrays_mmap_RowOriented outMat(outputPath, numOfRows, numOfCols, numOfSubunits, subunitId);

    if (verbose) {
        std::cerr << "copying ..." << std::endl;
    }
    for (TlMatrixObject::index_type r = 0; r < numOfRows; ++r) {
        for (TlMatrixObject::index_type c = 0; c < numOfCols; ++c) {
            outMat.set(r, c, inMat.get(r, c));
        }
    }
    if (verbose) {
        std::cerr << "done." << std::endl;
    }

    return EXIT_SUCCESS;
}
