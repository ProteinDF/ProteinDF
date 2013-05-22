#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "pdf-mat-add [options] MATRIX_FILE1 MATRIX_FILE2 OUTPUT" << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "h");

    std::string path1 = opt[1];
    std::string path2 = opt[2];
    std::string path3 = opt[3];

    bool isSymMat1 = false;
    if (TlSymmetricMatrix::isLoadable(path1) == true) {
        isSymMat1 = true;
    } else if (TlMatrix::isLoadable(path1) == true) {
        isSymMat1 = false;
    } else {
        std::cerr << "can not open file: " << path1 << std::endl;
        return EXIT_FAILURE;
    }

    bool isSymMat2 = false;
    if (TlSymmetricMatrix::isLoadable(path2) == true) {
        isSymMat2 = true;
    } else if (TlMatrix::isLoadable(path2) == true) {
        isSymMat2 = false;
    } else {
        std::cerr << "can not open file: " << path2 << std::endl;
        return EXIT_FAILURE;
    }

    if ((isSymMat1 == true) && (isSymMat2 == true)) {
        TlSymmetricMatrix M1;
        M1.load(path1);

        TlSymmetricMatrix M2;
        M2.load(path2);

        if ((M1.getNumOfRows() != M2.getNumOfRows()) ||
            (M1.getNumOfCols() != M2.getNumOfCols())) {
            std::cerr << TlUtils::format("size is not consistent: (%d, %d) != (%d, %d)",
                                         M1.getNumOfRows(), M1.getNumOfCols(),
                                         M2.getNumOfRows(), M2.getNumOfCols())
                      << std::endl;
            return EXIT_FAILURE;
        }

        M1 += M2;
        M1.save(path3);
    } else {
        TlMatrix M1;
        if (isSymMat1 == true) {
            TlSymmetricMatrix tmp;
            tmp.load(path1);
            M1 = tmp;
        } else {
            M1.load(path1);
        }

        TlMatrix M2;
        if (isSymMat2 == true) {
            TlSymmetricMatrix tmp;
            tmp.load(path2);
            M2 = tmp;
        } else {
            M2.load(path2);
        }

        if ((M1.getNumOfRows() != M2.getNumOfRows()) ||
            (M1.getNumOfCols() != M2.getNumOfCols())) {
            std::cerr << TlUtils::format("size is not consistent: (%d, %d) != (%d, %d)",
                                         M1.getNumOfRows(), M1.getNumOfCols(),
                                         M2.getNumOfRows(), M2.getNumOfCols())
                      << std::endl;
            return EXIT_FAILURE;
        }

        M1 += M2;
        M1.save(path3);
    }

    return EXIT_SUCCESS;
}

