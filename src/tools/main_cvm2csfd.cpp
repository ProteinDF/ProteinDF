#include <iostream>
#include "TlFile.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"
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
    {
        int subunitID = 0;
        const std::string inputPath0 =
            TlDenseMatrix_arrays_Object::getFileName(inputBaseName, subunitID);
        const bool isLoadable = TlDenseMatrix_arrays_Object::isLoadable(
            inputPath0, &numOfCols, &numOfRows, &numOfSubunits, &subunitID);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath0 << std::endl;
            return EXIT_FAILURE;
        }

        if (verbose) {
            std::cerr << "rows: " << numOfRows << std::endl;
            std::cerr << "cols: " << numOfCols << std::endl;
            std::cerr << "units: " << numOfSubunits << std::endl;
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
    for (int i = 0; i < numOfSubunits; ++i) {
        if (verbose) {
            std::cerr << TlUtils::format("%d / %d", i + 1, numOfSubunits)
                      << std::endl;
        }

        TlDenseGeneralMatrix_arrays_ColOriented m;
        m.load(inputBaseName, i);

        for (index_type c = i; c < numOfCols; c += numOfSubunits) {
            const TlDenseVector_Lapack vtr = m.getVector(c);
            fileMat.setColVector(c, vtr);

            TlUtils::progressbar(float(c) / numOfCols);
        }
        TlUtils::progressbar(1.0);
        std::cout << std::endl;
    }

    if (verbose) {
        std::cerr << "end." << std::endl;
    }
    return EXIT_SUCCESS;
}
