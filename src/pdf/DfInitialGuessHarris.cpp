#include <iostream>
#include <set>
#include <string>

#include "CnError.h"
#include "PdfUserInput.h"
#include "DfInitialguess.h"
#include "DfInitialGuessHarris.h"

#include "DfEri.h"
#include "DfCalcGrid.h"
#include "DfOverlap.h"

#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"
#include "TlCombineDensityMatrix.h"

DfInitialGuessHarris::DfInitialGuessHarris(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam)
{
    std::string harrisDbFile = "harris.mpac";
    const char* pdfHome = std::getenv("PDF_HOME");
    if (pdfHome != NULL) {
        harrisDbFile = TlUtils::format("%s/data/harris.mpac", pdfHome);
    }

    TlMsgPack mpac;
    mpac.load(harrisDbFile);
    this->pdfParam_harrisDB_ = mpac.getSerializeData();
 
    if (this->isWorkOnDisk_ == true) {
        this->logger(" initial guess by harris functional is built on disk.\n");
        TlMatrix::useMemManager(true);
    }
}


DfInitialGuessHarris::~DfInitialGuessHarris()
{
}


void DfInitialGuessHarris::main()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    const TlOrbitalInfo orbInfo_high(pdfParam["coordinates"],
                                     pdfParam["basis_sets"]);
    const int numOfAOs_high = orbInfo_high.getNumOfOrbitals();
    
    TlSerializeData pdfParam_low; // for low-level
    
    // set low-lebel geometry
    pdfParam_low["coordinates"]["_"] = pdfParam["coordinates"]["_"];
    
    // set low-level basis set
    pdfParam_low["basis_sets"] = this->pdfParam_harrisDB_["basis_sets"];
    const Fl_Geometry flGeom(pdfParam_low["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    
    if (pdfParam["save_harris_param"].getBoolean() == true) {
        TlMsgPack mpac(pdfParam_low);
        mpac.save("pdfparam_low.mpac");
    }

    
    // create low-level density matrix
    const TlOrbitalInfo orbInfo_low(pdfParam_low["coordinates"],
                                    pdfParam_low["basis_sets"]);
    const int numOfAOs_low = orbInfo_low.getNumOfOrbitals();
    TlSymmetricMatrix P_low(numOfAOs_low);
    TlCombineDensityMatrix combineDensMat;
    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const std::string atomSymbol = flGeom.getAtom(atomIndex);
        if (atomSymbol == "X") {
            continue;
        }
        
        TlSerializeData coord;
        coord["_"].pushBack( pdfParam["coordinates"]["_"].getAt(atomIndex));
        
        TlOrbitalInfo orbInfo_harrisDB(coord,
                                       this->pdfParam_harrisDB_["basis_sets"]);
        
        const TlSymmetricMatrix P_DB(this->pdfParam_harrisDB_["density_matrix"][atomSymbol]);
        combineDensMat.make(orbInfo_harrisDB, P_DB,
                            orbInfo_low, &P_low);
    }
    //P_low.save("P_low.mtx");

    // transform low-level density matrix to high-level one
    TlSymmetricMatrix P_high(numOfAOs_high);
    {
        DfOverlap ovp(this->pPdfParam_);
        TlMatrix S_tilde = ovp.getSpq(orbInfo_low, orbInfo_high);
        S_tilde.save("Stilde.mat");
        
        TlSymmetricMatrix S_inv;
        S_inv.load(this->getSpqMatrixPath());
        S_inv.inverse();

        TlMatrix omega = S_tilde * S_inv;
        TlMatrix omega_t = omega;
        omega_t.transpose();

        P_high = omega_t * P_low * omega;
    }

    P_high.save(this->getPpqMatrixPath(RUN_RKS, 0));
    this->logger(" initial density matrix is created using Harris functional.\n");
}


