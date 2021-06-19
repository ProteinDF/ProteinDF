#include <iostream>

#include "TlFile.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_mmap_roworiented.h"
// #include "tl_dense_general_matrix_arrays_roworiented.h"
#include "TlLogging.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_mmap.h"
#include "tl_matrix_utils.h"

typedef TlMatrixObject::index_type index_type;

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("%s [options] INPUT_BASENAME OUTPUT_PATH", progname.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

bool checkMatrixInfo(const std::string& rvmBasePath, TlMatrixObject::index_type* pNumOfRows,
                     TlMatrixObject::index_type* pNumOfCols, int* pNumOfSubunits, int* pSizeOfChunk, int* pSubunitID,
                     bool* pGuess) {
    TlLogging& log = TlLogging::getInstance();
    // TlMatrixObject::index_type numOfRows = 0;
    // TlMatrixObject::index_type numOfCols = 0;
    // int numOfSubunits = 0;
    // int sizeOfChunk = 0;

    // TlMatrixObject::index_type& sizeOfVector = numOfCols;
    // TlMatrixObject::index_type& numOfVectors = numOfRows;

    int subunitID = 0;

    TlMatrixObject::HeaderInfo headerInfo;
    bool isLoadable = false;

    // guess start unit
    {
        int subunitID = 0;
        const std::string inputPath0 = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, subunitID);
        log.debug(TlUtils::format("inputPath0: %s", inputPath0.c_str()));
        isLoadable = TlMatrixUtils::getHeaderInfo(inputPath0, &headerInfo);

        log.debug(TlUtils::format("isLoadable: %d", isLoadable));
        *pGuess = isLoadable;
    }
    // try direct load
    if (isLoadable == false) {
        log.debug(TlUtils::format("rvmBasePath: %s", rvmBasePath.c_str()));
        isLoadable = TlMatrixUtils::getHeaderInfo(rvmBasePath, &headerInfo);
        log.debug(TlUtils::format("isLoadable: %d", isLoadable));
    }

    if (isLoadable) {
        if (isVerbose) {
            std::cerr << "check successed." << std::endl;
        }
        // numOfRows = numOfVectors;
        // numOfCols = sizeOfVector;

        *pNumOfRows = headerInfo.numOfVectors;
        *pNumOfCols = headerInfo.sizeOfVector;
        *pNumOfSubunits = headerInfo.numOfSubunits;
        *pSizeOfChunk = headerInfo.sizeOfChunk;
        *pSubunitID = headerInfo.subunitId;
    }

    return isLoadable;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "ho:v");

    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    const bool verbose = (opt["v"] == "defined");

    if (opt.getCount() <= 2) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }

    TlLogging& log = TlLogging::getInstance();
    std::string output = "rvm2csfd.log";
    if (opt["o"].empty() != true) {
        output = opt["o"];
    }
    log.setFilePath(output);

    //
    std::string inputBasePath = opt[1];
    std::string outputPath = opt[2];
    log.info(TlUtils::format("input base path: %s", inputBasePath.c_str()));
    log.info(TlUtils::format("output: %s", outputPath.c_str()));

    // transpose2CSFD(inputBaseName, outputPath, verbose, true);
    TlMatrixObject::index_type numOfRows = 0;
    TlMatrixObject::index_type numOfCols = 0;
    int numOfSubunits = 0;
    int sizeOfChunk = 0;
    int subunitID = 0;
    bool sequentialMode = false;
    const bool canGetMatrixInfo = checkMatrixInfo(inputBasePath, &numOfRows, &numOfCols, &numOfSubunits, &sizeOfChunk,
                                                  &subunitID, &sequentialMode);
    if (canGetMatrixInfo == false) {
        log.critical(TlUtils::format("can not open file: %s", inputBasePath.c_str()));

        return EXIT_FAILURE;
    }

    log.info(TlUtils::format("rows: %d", numOfRows));
    log.info(TlUtils::format("cols: %d", numOfCols));
    log.info(TlUtils::format("units: %d", numOfSubunits));
    log.info(TlUtils::format("chunk: %d", sizeOfChunk));

    // prepare output file
    if (TlFile::isExistFile(outputPath)) {
        log.info(TlUtils::format("overwrite output matrix: %s", outputPath.c_str()));
        TlDenseGeneralMatrix_mmap outMat(outputPath);
        if ((outMat.getNumOfRows() != numOfRows) || (outMat.getNumOfCols() != numOfCols)) {
            log.critical(TlUtils::format("matrix size mismatch (%d, %d) != (%d, %d).", numOfRows, numOfCols,
                                         outMat.getNumOfRows(), outMat.getNumOfCols()));
            return EXIT_FAILURE;
        }
    } else {
        log.info(TlUtils::format("create output: %s", outputPath.c_str()));
        TlDenseGeneralMatrix_mmap outMat(outputPath, numOfRows, numOfCols);
        if (verbose) {
            std::cerr << "output matrix has been prepared by mmap." << std::endl;
        }
    }

    //
    TlDenseGeneralMatrix_mmap outMat(outputPath);
    if (sequentialMode) {
        log.info("run sequential mode");
        const bool showProgress = true;
        for (int unit = 0; unit < numOfSubunits; ++unit) {
            log.info(TlUtils::format(">>>> unit %d/%d", unit, numOfSubunits - 1));
            const std::string inputPath = TlDenseMatrix_arrays_mmap_Object::getFileName(inputBasePath, unit);

            log.info(TlUtils::format("load matrix: %s", inputPath.c_str()));
            TlDenseGeneralMatrix_arrays_mmap_RowOriented inMat(inputPath);
            log.info("convert layout");
            inMat.convertMemoryLayout("", verbose, showProgress);
            log.info("set");
            inMat.set2csfd(&outMat, verbose, showProgress);
        }
    } else {
        log.info("run unit mode");
        std::cerr << TlUtils::format(">>>> unit %d/%d", subunitID, numOfSubunits - 1) << std::endl;
        TlDenseGeneralMatrix_arrays_mmap_RowOriented inMat(inputBasePath);
        inMat.convertMemoryLayout("", verbose, verbose);
        inMat.set2csfd(&outMat, verbose, verbose);
    }

    if (verbose) {
        std::cerr << "end." << std::endl;
    }
    return EXIT_SUCCESS;
}
