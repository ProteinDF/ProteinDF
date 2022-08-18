#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DfOverlapX.h"
#include "Fl_Geometry.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlUtils.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#endif  // HAVE_LAPACK

void usage(const std::string& progname) {
    std::cout << TlUtils::format("usage: %s \"pdfparam1.mpac\" \"C matrix path1\" \"pdfparam2.mpac\" \"C matrix path2\"",
                                 progname.c_str())
              << std::endl;
    std::cout << " calculate CSC matrix" << std::endl;
    std::cout << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "-s FILE:    save S~ matrix" << std::endl;
    std::cout << "-v:         verbose output" << std::endl;
}

template <class GeneralMatrix>
void calcCSC(const TlSerializeData& pdfParam1, const std::string& cMatPath1,
             const TlSerializeData& pdfParam2, const std::string& cMatPath2,
             const std::string& sTildeMatPath, const std::string& cscMatPath, const bool verbose) {
    const TlOrbitalInfo orbInfo1(pdfParam1["coordinates"], pdfParam1["basis_set"]);
    const TlOrbitalInfo orbInfo2(pdfParam2["coordinates"], pdfParam2["basis_set"]);

    GeneralMatrix S_tilde;
    {
        TlSerializeData dummy;
        DfOverlapX ovp(&dummy);
        ovp.getTransMat(orbInfo1, orbInfo2, &S_tilde);
        if (sTildeMatPath.empty() == false) {
            std::cerr << "S_tilde save: " << sTildeMatPath << std::endl;
            S_tilde.save(sTildeMatPath);
        }
    }
    if (verbose) {
        std::cerr << TlUtils::format("S~ size: %d, %d", S_tilde.getNumOfRows(), S_tilde.getNumOfCols()) << std::endl;
    }

    GeneralMatrix C1;
    if (verbose) {
        std::cerr << "load C matrix1: " << cMatPath1 << std::endl;
    }
    C1.load(cMatPath1);
    if (verbose) {
        std::cerr << TlUtils::format("C matrix size: %d, %d", C1.getNumOfRows(), C1.getNumOfCols()) << std::endl;
    }

    GeneralMatrix C2;
    if (verbose) {
        std::cerr << "load C matrix2: " << cMatPath2 << std::endl;
    }
    C2.load(cMatPath2);
    if (verbose) {
        std::cerr << TlUtils::format("C matrix size: %d, %d", C2.getNumOfRows(), C2.getNumOfCols()) << std::endl;
    }

    GeneralMatrix C1t = C1.transpose();
    GeneralMatrix CSC = C1t * S_tilde * C2;
    CSC.save(cscMatPath);
}

int main(int argc, char* argv[]) {
    // setup file path
    TlGetopt opt(argc, argv, "ho:s:v");

    if ((opt["h"].empty() == false) || (opt.getCount() <= 2)) {
        usage(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string cscMatPath = "CSC.mat";
    if (opt["o"].empty() == false) {
        cscMatPath = opt["o"];
    }

    std::string sTildeMatPath = "Stilde.mat";
    if (opt["s"].empty() == false) {
        sTildeMatPath = opt["s"];
    }

    const bool verbose = (opt["v"] == "defined");

    const std::string paramPath1 = opt[1];
    const std::string cMatPath1 = opt[2];
    const std::string paramPath2 = opt[3];
    const std::string cMatPath2 = opt[4];

    if (verbose) {
        std::cerr << "ProteinDF parameter file1: " << paramPath1 << std::endl;
        std::cerr << "ProteinDF parameter file2: " << paramPath2 << std::endl;
        std::cerr << "C matrix path1: " << cMatPath1 << std::endl;
        std::cerr << "C matrix path2: " << cMatPath2 << std::endl;
    }

    TlSerializeData pdfParam1;
    {
        TlMsgPack mpac;
        mpac.load(paramPath1);
        pdfParam1 = mpac.getSerializeData();
    }
    TlSerializeData pdfParam2;
    {
        TlMsgPack mpac;
        mpac.load(paramPath2);
        pdfParam2 = mpac.getSerializeData();
    }

#ifdef HAVE_LAPACK
    calcCSC<TlDenseGeneralMatrix_Lapack>(pdfParam1, cMatPath1, pdfParam2, cMatPath2, sTildeMatPath, cscMatPath, verbose);
#endif  // HAVE_LAPACK

    return EXIT_SUCCESS;
}
