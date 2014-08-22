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

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "TlGetopt.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "Fl_Geometry.h"
#include "DfOverlapX.h"
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
    
    const TlOrbitalInfo orbInfo_old(pdfParamOld["coordinates"],
                                    pdfParamOld["basis_sets"]);
    const TlOrbitalInfo orbInfo_new(pdfParamNew["coordinates"],
                                    pdfParamNew["basis_sets"]);

    TlMatrix omega;
    {
        // S_tilde
        DfOverlapX ovp(&pdfParamNew);
        TlMatrix S_tilde;
        ovp.getTransMat(orbInfo_old, orbInfo_new, &S_tilde);
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
