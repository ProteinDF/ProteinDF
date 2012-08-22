#include <iostream>
#include <limits>
#include "CnError.h"
#include "DfInitialGuessHuckel.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

const double DfInitialGuessHuckel::EPS = std::numeric_limits<double>::epsilon();

DfInitialGuessHuckel::DfInitialGuessHuckel(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam)
{
    this->initialize();
}

DfInitialGuessHuckel::~DfInitialGuessHuckel()
{
}

void DfInitialGuessHuckel::initialize()
{
    this->m_Hii.clear();

    this->m_Hii["H"][0] = -13.6; // H
    this->m_Hii["C"][0] = -21.4; // C
    this->m_Hii["C"][1] = -11.4;
    this->m_Hii["N"][0] = -26.0; // N
    this->m_Hii["N"][1] = -13.4;
    this->m_Hii["O"][0] = -32.3; // O
    this->m_Hii["O"][1] = -14.8;
    this->m_Hii["S"][0] = -20.0; // S
    this->m_Hii["S"][1] = -11.0;
    this->m_Hii["S"][2] = - 8.0;
}


void DfInitialGuessHuckel::createGuess()
{
    this->createGuess<TlSymmetricMatrix, TlMatrix>();
}


double DfInitialGuessHuckel::getHii(const std::string& sAtomName, const int nOrbitalType)
{
    double dAnswer = 0.0;

    std::map<std::string, std::map<int, double> >::const_iterator p = this->m_Hii.find(sAtomName);
    if (p != this->m_Hii.end()) {
        std::map<int, double>::const_iterator q = p->second.find(nOrbitalType);
        if (q != p->second.end()) {
            dAnswer = q->second;
        }
    }

    return dAnswer;
}

// TlSymmetricMatrix DfInitialGuessHuckel::getHpqMatrix()
// {
//     TlSymmetricMatrix Hpq = DfObject::getHpqMatrix<TlSymmetricMatrix>();
//     return Hpq;
// }

// output occupation number to a files in fl_Work directory
TlVector DfInitialGuessHuckel::generateOccupation(const RUN_TYPE nRunType, const int nNumOrbcut, const int nNumOfElectrons)
{
    TlVector occ(nNumOrbcut);

    const int nElec_2 = (this->m_nMethodType == METHOD_RKS) ? (nNumOfElectrons / 2) : nNumOfElectrons;
    std::cerr << "elec = " << nElec_2 << std::endl;

    const double dElectrons = (this->m_nMethodType == METHOD_RKS) ? 2.0 : 1.0;
    for (int i = 0; i <nElec_2; i++) {
        occ[i] = dElectrons;
    }

    std::string sSuffix = "";
    switch (nRunType) {
    case RUN_RKS:
        //sSuffix = "";
        break;
    case RUN_UKS_ALPHA:
        sSuffix = "_Alpha";
        break;
    case RUN_UKS_BETA:
        sSuffix = "_Beta";
        break;
    case RUN_ROKS:
        CnErr.abort("to implement. DfInitialGuessHuckel::generateOccupation()");
        break;
    default:
        CnErr.abort("unknown run type.");
        break;
    }

    const std::string sFileName = "fl_Work/fl_Occupation" + sSuffix;
    occ.save(sFileName);

    return occ;
}

