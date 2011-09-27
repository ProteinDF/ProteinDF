#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "TlGetopt.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "Fl_Geometry.h"
#include "Fl_Gto_Orbital.h"
#include "DfOverlap.h"
#include "TlMsgPack.h"

void usage()
{
    std::cout << "usage: transdim \"old matrix file\" \"new matrix file\"" << std::endl;
    std::cout << " transform matrix into the different basis." << std::endl;
    std::cout << " default old matrix file is old.mtx, and new matrix file is new.mtx." << std::endl;
    std::cout << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "-o FILE:    old ProteinDF parameter file (default: pdfparam.mpac.old)" << std::endl;
    std::cout << "-n FILE:    new ProteinDF parameter file (default: pdfparam.mpac)" << std::endl;
    std::cout << "-s FILE:    overlap matrix (default: fl_Work/fl_Mtr_Spq.matrix)" << std::endl;
}

int main(int argc, char* argv[])
{
    // setup file path
    TlGetopt opt(argc, argv, "hn:o:s:S:");

    if (opt["h"].empty() == false) {
        usage();
        return EXIT_SUCCESS;
    }
    
    std::string paramPath_old = "pdfparam.mpac.old";
    if (opt["o"].empty() == false) {
        paramPath_old = opt["o"];
    }

    std::string paramPath_new = "pdfparam.mpac";
    if (opt["n"].empty() == false) {
        paramPath_new = opt["n"];
    }
    
    std::string spqMatrixPath = "fl_Work/fl_Mtr_Spq.matrix";
    if (opt["s"].empty() == false) {
        spqMatrixPath = opt["s"];
    }
    
    std::string oldMatrixPath = "old.mtx";
    if (opt.getCount() > 1) {
        oldMatrixPath = opt[1];
    }

    std::string newMatrixPath = "new.mtx";
    if (opt.getCount() > 2) {
        newMatrixPath = opt[2];
    }

    std::string sTildePath = "";
    if (opt["S"].empty() == false) {
        sTildePath = opt["S"];
    }
    
    std::cerr << "old ProteinDF parameter file: " << paramPath_old << std::endl;
    std::cerr << "new ProteinDF parameter file: " << paramPath_new << std::endl;
    std::cerr << "overlap matrix path: " << spqMatrixPath << std::endl;
    std::cerr << "old matrix path: " << oldMatrixPath << std::endl;
    std::cerr << "new matrix path: " << newMatrixPath << std::endl;
    
    TlSerializeData pdfParamOld;
    {
        TlMsgPack mpac;
        mpac.load(paramPath_old);
        pdfParamOld = mpac.getSerializeData();
    }
    TlSerializeData pdfParamNew;
    {
        TlMsgPack mpac;
        mpac.load(paramPath_new);
        pdfParamNew = mpac.getSerializeData();
    }
    
    const TlOrbitalInfo orbInfo_old(pdfParamOld["model"]["coordinates"],
                                    pdfParamOld["model"]["basis_set"]);
    const TlOrbitalInfo orbInfo_new(pdfParamNew["model"]["coordinates"],
                                    pdfParamNew["model"]["basis_set"]);

    TlMatrix omega;
    {
        // S_tilde
        DfOverlap ovp(&pdfParamNew);
        TlMatrix S_tilde = ovp.getSpq(orbInfo_old, orbInfo_new);
        if (sTildePath.empty() == false) {
            std::cerr << "S_tilde save: " << sTildePath << std::endl;
            S_tilde.save(sTildePath);
        }
        
        TlSymmetricMatrix S_inv;
        S_inv.load(spqMatrixPath);
        S_inv.inverse();

        omega = S_tilde * S_inv;
    }
    
    TlMatrix omega_t = omega;
    omega_t.transpose();
    
    TlSymmetricMatrix D;
    D.load(oldMatrixPath);
    TlSymmetricMatrix newD = omega_t * D * omega;
    newD.save(newMatrixPath);

    return EXIT_SUCCESS;
}
