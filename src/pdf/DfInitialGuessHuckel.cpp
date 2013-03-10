#include <iostream>
#include <limits>
#include "CnError.h"
#include "DfInitialGuessHuckel.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

const double DfInitialGuessHuckel::EPS = std::numeric_limits<double>::epsilon();

DfInitialGuessHuckel::DfInitialGuessHuckel(TlSerializeData* pPdfParam)
    : DfInitialGuess(pPdfParam)
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

