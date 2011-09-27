#include <cassert>
#include "DfXcenefitting.h"
#include "TlTime.h"
#include "TlVector.h"
#include "TlLogX.h"

DfXcenefitting::DfXcenefitting(TlSerializeData* pPdfParam, int nItr)
    : DfObject(pPdfParam)
{
    //this->number_iteration = nItr;

    //this->scftype = flGbi["SCF"]["method"];
    //this->naux = atoi(flGbi["SCF"]["control-nauxxc"].c_str());
    const int number_rotation = (*pPdfParam)["model"]["file-rot-number"].getInt();
    this->m_nIteration  = this->m_nIteration % number_rotation;
    if (this->m_nIteration == 0) {
        this->m_nIteration = number_rotation;
    }
    //this->bMemorySave      = (flGbi["SCF"]["scf-memory-saving"] == "yes") ? true : false;
}

DfXcenefitting::~DfXcenefitting()
{
}

int DfXcenefitting::dfXceMain()
{
    TlLogX& Log = TlLogX::getInstance();

    Log << "=== START ==========\n";
    Log << "START : " << TlTime::getNow() << "\n";

    this->calcEpsilon();
    Log << "calcEpsilon : " << TlTime::getNow() << "\n";
    Log.flush();

    Log << "END : " << TlTime::getNow() << "\n";
    Log << "=== END ============\n";
    Log.flush();

    return 0;
}

// Read MyuGamma.
// Scale MyuGamma to get EpsilonGamma.
// File out EpsilonGamma.
int DfXcenefitting::calcEpsilon()
{
//     TlLogX& Log = TlLogX::getInstance();
    
    if (this->m_nMethodType == METHOD_RKS) {
        TlVector myu;
        myu.load("fl_Work/fl_Vct_Myu" + TlUtils::xtos(this->m_nIteration));

        // temporal print ==> to fl_Globaloutput
//         if (outlevel < -5) {
//             Log << "*** PRINT OUT MyuGamma in readPpq ***\n";
//             for (int j = 0; j < myu.getSize(); ++j) {
//                 Log << "MyuGamma[" << j << "] = " << myu[j] << "\n";
//             }
//             Log << "\n";
//         }

        TlVector eps = 0.75 * myu;

        if (this->m_bMemorySave == false) {
            eps.save("fl_Work/fl_Vct_Epsilon" + TlUtils::xtos(this->m_nIteration));
        } else {
            eps.save("fl_Work/fl_Vct_Epsilon");
        }

//         if (outlevel < -3) {
//             Log << "*** PRINT OUT EpsilonGamma in calcMN ***\n";
//             for (int i = 0; i < this->naux; ++i) {
//                 Log << "EpsilonGamma[" << i << "] = " << eps[i] << "\n";
//             }
//             Log << "\n";
//         }
    } else {
        TlVector myuA;
        myuA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(this->m_nIteration));

        TlVector myuB;
        myuB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(this->m_nIteration));

        // temporal print ==> to fl_Globaloutput
//         if (outlevel < -5) {
//             Log << "*** PRINT OUT MyuGammaA in readPpq ***\n";
//             for (int k = 0; k < myuA.getSize(); ++k) {
//                 Log << "MyuGammaA[" << k << "] = " << myuA[k] << "\n";
//             }
//             Log << "\n";

//             Log << "*** PRINT OUT MyuGammaB in readPpq ***\n";
//             for (int l = 0; l < myuB.getSize(); ++l) {
//                 Log << "MyuGammaB[" << l << "] = " << myuB[l] << "\n";
//             }
//             Log << "\n";
//         }

        TlVector epsA = 0.75 * myuA;
        TlVector epsB = 0.75 * myuB;

//     if(this->bMemorySave == false){
//       epsA.save("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(niteration));
//       epsB.save("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(niteration));
//     } else {
        epsA.save("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(this->m_nIteration));
        epsB.save("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(this->m_nIteration));
//     }

//         if (outlevel < -3) {
//             Log << "*** PRINT OUT EpsilonGammaA in calcMN ***\n";
//             for (int i = 0; i < this->naux; ++i) {
//                 Log << "EpsilonGammaA[" << i << "] = " << epsA[i] << "\n";
//             }
//             Log << "\n";
//             Log << "*** PRINT OUT EpsilonGammaB in calcMN ***\n";
//             for (int k = 0; k < this->naux; ++k) {
//                 Log << "EpsilonGammaB[" << k << "] = " << epsB[k] << "\n";
//             }
//             Log << "\n";
//         }
    }

    return 0;
}

