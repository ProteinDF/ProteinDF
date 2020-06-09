#include <iostream>

#include "TlFile.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_mmap.h"

typedef TlMatrixObject::index_type index_type;

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("%s [options] INPUT_BASENAME OUTPUT_BASENAME",
                                 progname.c_str())
              << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
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
    std::string inputBaseName = opt[1];
    std::string outputPath = opt[2];
    if (verbose) {
        std::cerr << "input base path: " << inputBaseName << std::endl;
        std::cerr << "output: " << outputPath << std::endl;
    }

    // check
    index_type numOfRows = 0;
    index_type numOfCols = 0;
    int numOfSubunits = 0;
    int sizeOfChunk = 0;
    {
        int subunitID = 0;
        const std::string inputPath0 =
            TlDenseMatrix_arrays_Object::getFileName(inputBaseName, subunitID);
        index_type sizeOfVector, numOfVectors;
        const bool isLoadable = TlDenseMatrix_arrays_Object::isLoadable(
            inputPath0, &numOfVectors, &sizeOfVector, &numOfSubunits,
            &subunitID, &sizeOfChunk);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath0 << std::endl;
            return EXIT_FAILURE;
        }

        // row-vector matrix
        numOfRows = numOfVectors;
        numOfCols = sizeOfVector;

        if (verbose) {
            std::cerr << "rows: " << numOfRows << std::endl;
            std::cerr << "cols: " << numOfCols << std::endl;
            std::cerr << "units: " << numOfSubunits << std::endl;
            std::cerr << "size of chunk: " << sizeOfChunk << std::endl;
        }
    }

    // prepare output
    if (TlFile::isExistFile(outputPath)) {
        if (verbose) {
            std::cerr << "file overwrite: " << outputPath << std::endl;
        }
        TlFile::remove(outputPath);
    }
    TlDenseGeneralMatrix_mmap fileMat(outputPath, numOfRows, numOfCols);

    // load & set
    for (int unit = 0; unit < numOfSubunits; ++unit) {
        if (verbose) {
            std::cerr << TlUtils::format("%d / %d", unit + 1, numOfSubunits)
                      << std::endl;
        }

        TlDenseGeneralMatrix_arrays_RowOriented m;
        m.load(inputBaseName, unit);

        // std::vector<double> vtr(numOfCols);
        // for (index_type r = i; r < numOfRows; r += numOfSubunits) {
        //     m.getVector(r, &(vtr[0]), numOfCols);
        //     fileMat.setRowVector(r, vtr);
        //     TlUtils::progressbar(float(r) / numOfRows);
        // }
        std::vector<double> chunkBuf(numOfCols * sizeOfChunk);
        std::vector<double> transBuf(numOfCols * sizeOfChunk);
        const int numOfLocalChunks =
            TlDenseMatrix_arrays_Object::getNumOfLocalChunks(
                numOfRows, numOfSubunits, sizeOfChunk);
        for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
            const index_type chunkStartRow =
                sizeOfChunk * (numOfSubunits * chunk + unit);
            m.getChunk(chunkStartRow, &(chunkBuf[0]), numOfCols * sizeOfChunk);

            // change memory layout
            const index_type readRowChunks =
                std::min(sizeOfChunk, numOfRows - chunkStartRow);
            TlUtils::changeMemoryLayout(&(chunkBuf[0]), readRowChunks,
                                        numOfCols, &(transBuf[0]));

            TlDenseGeneralMatrix_Eigen tmpMat(readRowChunks, numOfCols,
                                              &(transBuf[0]));
            fileMat.block(chunkStartRow, 0, tmpMat);

            TlUtils::progressbar(float(chunk) / numOfLocalChunks);
        }

        TlUtils::progressbar(1.0);
        std::cout << std::endl;
    }

    if (verbose) {
        std::cerr << "end." << std::endl;
    }
    return EXIT_SUCCESS;
}
