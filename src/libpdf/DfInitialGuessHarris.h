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

#ifndef DFINITIALGUESSHARRIS_H
#define DFINITIALGUESSHARRIS_H

#include <string>
#include "DfObject.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "TlMsgPack.h"
#include "TlCombineDensityMatrix.h"

/// Harrisの汎関数による初期値を作成する
class DfInitialGuessHarris : public DfObject {
private:
    enum ScfType {
        RKS,
        UKS,
        ROKS
    };

public:
    DfInitialGuessHarris(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHarris();

    virtual void main();

protected:
    void loadHarrisDB();

    template<class MatrixType, class SymmetricMatrixType, class DfOverlapType, class DfPopulationType>
    void calcInitialDensityMatrix();

protected:
    /// debug時はtrueにする
    bool debug_;

    TlSerializeData pdfParam_harrisDB_;
};


template<class MatrixType, class SymmetricMatrixType,
         class DfOverlapType, class DfPopulationType>
void DfInitialGuessHarris::calcInitialDensityMatrix()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    const TlOrbitalInfo orbInfo_high(pdfParam["coordinates"],
                                     pdfParam["basis_set"]);
    const int numOfAOs_high = orbInfo_high.getNumOfOrbitals();
    
    TlSerializeData pdfParam_low; // for low-level
    
    // set low-lebel geometry
    pdfParam_low["coordinates"]["atoms"] = pdfParam["coordinates"]["atoms"];
    
    // set low-level basis set
    pdfParam_low["basis_set"] = this->pdfParam_harrisDB_["basis_set"];
    const Fl_Geometry flGeom(pdfParam_low["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    
    if (pdfParam["save_harris_param"].getBoolean() == true) {
        TlMsgPack mpac(pdfParam_low);
        mpac.save("pdfparam_low.mpac");
    }
    
    // create low-level density matrix
    const TlOrbitalInfo orbInfo_low(pdfParam_low["coordinates"],
                                    pdfParam_low["basis_set"]);
    const int numOfAOs_low = orbInfo_low.getNumOfOrbitals();
    SymmetricMatrixType P_low(numOfAOs_low);
    TlCombineDensityMatrix combineDensMat;
    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const std::string atomSymbol = flGeom.getAtomSymbol(atomIndex);
        if (atomSymbol == "X") {
            continue;
        }
        
        TlSerializeData coord;
        coord["atoms"].pushBack(pdfParam["coordinates"]["atoms"].getAt(atomIndex));
        
        TlOrbitalInfo orbInfo_harrisDB(coord,
                                       this->pdfParam_harrisDB_["basis_set"]);
        
        const TlSymmetricMatrix P_DB(this->pdfParam_harrisDB_["density_matrix"][atomSymbol]);
        combineDensMat.make(orbInfo_harrisDB, P_DB,
                            orbInfo_low, &P_low);
    }
    if (this->debug_) {
        P_low.save("P_low.mtx");
    }

    // transform low-level density matrix to high-level one
    SymmetricMatrixType P_high(numOfAOs_high);
    {
        DfOverlapType ovp(this->pPdfParam_);
        MatrixType S_tilde;
        ovp.getTransMat(orbInfo_low, orbInfo_high, &S_tilde);
        // S_tilde.save("S_tilde.mat");
        
        SymmetricMatrixType S_inv;
        S_inv.load(this->getSpqMatrixPath());
        S_inv.inverse();

        MatrixType omega = S_tilde * S_inv;
        MatrixType omega_t = omega;
        omega_t.transpose();

        if (this->debug_) {
            omega.save("omega.mtx");
        }
        P_high = omega_t * P_low * omega;
    }

    // normalize
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            // double numOfElectrons = 0.0;
            DfPopulationType dfPop(this->pPdfParam_);
            // this->savePpqMatrix(RUN_RKS, 0, P_high); // sumOfElectrons()で必要
            // dfPop.sumOfElectrons(0, &numOfElectrons, NULL);
            const double numOfElectrons = dfPop.getSumOfElectrons(P_high);
            const double coef = this->m_nNumOfElectrons / numOfElectrons;

            this->savePpqMatrix(RUN_RKS, 0, coef * P_high);
        }
        break;

    case METHOD_UKS:
        {
            DfPopulationType dfPop(this->pPdfParam_);
            double numOfAlphaElectrons = dfPop.getSumOfElectrons(P_high);
            double numOfBetaElectrons = dfPop.getSumOfElectrons(P_high);
            const double coef_alpha = this->m_nNumOfAlphaElectrons / numOfAlphaElectrons;
            const double coef_beta  = this->m_nNumOfBetaElectrons  / numOfBetaElectrons;

            this->savePpqMatrix(RUN_UKS_ALPHA, 0, coef_alpha * P_high);
            this->savePpqMatrix(RUN_UKS_BETA,  0, coef_beta  * P_high);
        }
        break;

    case METHOD_ROKS:
        {
            DfPopulationType dfPop(this->pPdfParam_);
            double numOfCloseElectrons = dfPop.getSumOfElectrons(P_high);
            double numOfOpenElectrons = dfPop.getSumOfElectrons(P_high);
            const double coef_close = this->numOfClosedShellElectrons_ / numOfCloseElectrons;
            const double coef_open = this->numOfOpenShellElectrons_  / numOfOpenElectrons;

            this->savePpqMatrix(RUN_ROKS_CLOSED, 0, coef_close * P_high);
            this->savePpqMatrix(RUN_ROKS_OPEN,   0, coef_open  * P_high);
        }
        break;

    default:
        this->log_.critical(TlUtils::format("program error: %s %s", __FILE__, __LINE__));
        abort();
        break;
    }

    this->logger(" initial density matrix is created using Harris functional.\n");
}


#endif // DFINITIALGUESSHARRIS_H
