#include <iostream>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hs:v");

    bool bVerbose = (opt["v"] == "defined");

    const std::string sPath1 = opt[1];
    const std::string sPath2 = opt[2];
    if (bVerbose) {
        std::cerr << "loading... " << sPath1 << std::endl;
        std::cerr << "loading... " << sPath2 << std::endl;
    }

    std::string savePath = "";
    if (opt["s"].empty() == false) {
        savePath = opt["s"];
    }

    int nErrorCode = 0;

    std::ifstream ifs1;
    ifs1.open(sPath1.c_str());
    if (ifs1.fail()) {
        std::cerr << "could not open file. " << sPath1 << std::endl;
        return 1;
    }

    std::ifstream ifs2;
    ifs2.open(sPath2.c_str());
    if (ifs2.fail()) {
        std::cerr << "could not open file. " << sPath2 << std::endl;
        return 1;
    }


    if (TlSymmetricMatrix::isLoadable(ifs1) == true) {
        if (TlSymmetricMatrix::isLoadable(ifs2) == true) {
            TlSymmetricMatrix m1, m2;
            m1.load(sPath1);
            m2.load(sPath2);

            m1 -= m2;
            if (savePath.empty() == false) {
                m1.save(savePath);
            } else {
                m1.print(std::cout);
            }
        } else {
            std::cerr << "could not open: " << sPath2 << std::endl;
            nErrorCode = 1;
        }
    } else if (TlMatrix::isLoadable(ifs1) == true) {
        if (TlMatrix::isLoadable(ifs2) == true) {
            TlMatrix m1, m2;
            m1.load(sPath1);
            m2.load(sPath2);

            m1 -= m2;
            if (savePath.empty() == false) {
                m1.save(savePath);
            } else {
                m1.print(std::cout);
            }
        } else {
            std::cerr << "could not open: " << sPath2 << std::endl;
            nErrorCode = 1;
        }
    } else {
        std::cerr << "could not open: " << sPath1 << std::endl;
        nErrorCode = 1;
    }

    return nErrorCode;
}


