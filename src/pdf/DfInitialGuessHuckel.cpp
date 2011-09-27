#include <iostream>
#include <limits>
#include "CnError.h"
#include "DfInitialGuessHuckel.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlOrbitalInfo.h"

const double DfInitialGuessHuckel::EPS = std::numeric_limits<double>::epsilon();

DfInitialGuessHuckel::DfInitialGuessHuckel(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam)
{
    this->initialize();
    this->createGuess();
}

DfInitialGuessHuckel::~DfInitialGuessHuckel()
{
    //std::cerr << "<<<< DfInitialGuessHuckel::~DfInitialGuessHuckel()" << std::endl;
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
    TlSymmetricMatrix F;
    switch (initialGuessType_) {
    case GUESS_HUCKEL:
        F = this->getHuckelMatrix();
        break;
    case GUESS_CORE:
        F = this->getHpqMatrix();
        break;
    default:
        CnErr.abort("unknown guess type");
        break;
    }

    TlMatrix X;
    X.load("fl_Work/fl_Mtr_X.matrix");
    {
        TlMatrix tX = X;
        tX.transpose();

        F = tX * F * X;
    }

    TlMatrix C;
    {
        TlVector eigval;
        F.diagonal(&eigval, &C);
    }

    C = X * C;

    // P
    switch (this->m_nMethodType) {
    case METHOD_RKS: {
        const TlVector occ = this->generateOccupation(RUN_RKS, this->m_nNumOfMOs, this->m_nNumOfElectrons);
        TlSymmetricMatrix P, P1, P2;
        this->generatePMatrix(C, occ, P1, P2);

        //P1.save(this->getP1pqMatrixPath(RUN_RKS, 0));
        //P2.save(this->getP2pqMatrixPath(RUN_RKS, 0));
        //DfObject::savePCMatrix(P1, 0);
        //DfObject::savePOMatrix(P2, 0);
        P = 2.0 * P1 + P2;
        //P.save(this->getPpqMatrixPath(RUN_RKS, 0));
        DfObject::savePpqMatrix(RUN_RKS, 0, P);
    }
    break;

    case METHOD_UKS: {
        const TlVector occA = this->generateOccupation(RUN_UKS_ALPHA, this->m_nNumOfMOs, this->m_nNumOfAlphaElectrons);
        const TlVector occB = this->generateOccupation(RUN_UKS_BETA, this->m_nNumOfMOs, this->m_nNumOfBetaElectrons);

        TlSymmetricMatrix P, P1, P2;
        this->generatePMatrix(C, occA, P1, P2);
        //P1.save(this->getP1pqMatrixPath(RUN_UKS_ALPHA, 0));
        //P2.save(this->getP2pqMatrixPath(RUN_UKS_ALPHA, 0));
        P = 2.0 * P1 + P2;
        //P.save(this->getPpqMatrixPath(RUN_UKS_ALPHA, 0));
        DfObject::savePpqMatrix(RUN_UKS_ALPHA, 0, P);

        this->generatePMatrix(C, occB, P1, P2);
        //P1.save(this->getP1pqMatrixPath(RUN_UKS_BETA, 0));
        //P2.save(this->getP2pqMatrixPath(RUN_UKS_BETA, 0));
        P = 2.0 * P1 + P2;
        //P.save(this->getPpqMatrixPath(RUN_UKS_BETA, 0));
        DfObject::savePpqMatrix(RUN_UKS_BETA, 0, P);
    }
    break;

    case METHOD_ROKS:
        CnErr.abort("to implement");
        break;

    default:
        CnErr.abort("unknown method");
        break;
    }

}

TlSymmetricMatrix DfInitialGuessHuckel::getHuckelMatrix()
{
    TlOrbitalInfo orbInfo((*this->pPdfParam_)["model"]["coordinates"],
                          (*this->pPdfParam_)["model"]["basis_set"]);

    const int nNumOfOrbs = orbInfo.getNumOfOrbitals();

    TlSymmetricMatrix Hpq;
    Hpq.load("fl_Work/fl_Mtr_Hpq.matrix");
    assert(nNumOfOrbs == Hpq.getNumOfRows());

    TlSymmetricMatrix Spq;
    Spq.load("fl_Work/fl_Mtr_Spq.matrix");
    assert(nNumOfOrbs == Spq.getNumOfRows());

    const double K = 1.75;

    // create Fock matrix using modified extended Huckel
    TlSymmetricMatrix F(nNumOfOrbs);
    int nPrevAtomI = -1;
    int nPrevShellType = -1;
    for (int i=0; i < nNumOfOrbs; i++) {
        const int nAtomI = orbInfo.getAtomIndex(i);
        const int nShellTypeI = orbInfo.getShellType(i);
        const double Hii = this->getHii(orbInfo.getAtomName(i), nShellTypeI);

        // i != j
        for (int j=0; j < i; j++) {
            const int nAtomJ = orbInfo.getAtomIndex(j);
            if (nAtomI != nAtomJ) {
                const double Hjj = this->getHii(orbInfo.getAtomName(j), orbInfo.getShellType(j));
                F(i, j) = 0.5 * K * (Hii + Hjj) * Spq(i, j);
            }
        }

        // i == j
        if ((nPrevAtomI != nAtomI) || (nPrevShellType != nShellTypeI)) {
            F(i, i) = Hii;
            nPrevAtomI = nAtomI;
            nPrevShellType = nShellTypeI;
        }
    }

    F.print(std::cout);

    return F;
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

TlSymmetricMatrix DfInitialGuessHuckel::getHpqMatrix()
{
    TlSymmetricMatrix Hpq;
    Hpq.load("fl_Work/fl_Mtr_Hpq.matrix");

    return Hpq;
}

// TlMatrix DfInitialGuessHuckel::getCmatrix(const TlSymmetricMatrix& F){
//   TlMatrix X;
//   X.load("fl_Work/fl_Mtr_X.matrix");

//   TlMatrix Cprime;
//   {
//     TlMatrix Xt = X;
//     Xt.transpose();

//     TlSymmetricMatrix Fprime = Xt * F * X;

//     TlVector eigval;
//     Fprime.diagonal(&eigval, &Cprime);
//   }

//   TlMatrix C = X * Cprime;

//   switch (this->m_nScfType){
//   case RKS:
//     C.save("fl_Work/fl_Mtr_C.matrix.rks0");
//     break;
//   default:
//     CnErr.abort("unknown SCF type. stop.");
//   }

//   return C;
// }

// TlMatrix DfInitialGuessHuckel::getCmatrix(const TlSymmetricMatrix& F){
//   TlMatrix X;
//   X.load("fl_Work/fl_Mtr_X.matrix");

//   F *= X;
//   X.transpose();
//   X *= F;

//   TlMatrix C;
//   {
//     TlVector eigval;
//     F.diagonal(&eigval, &C);
//   }

//   C = X * C;

//   switch (this->m_nMethodType){
//   case METHOD_RKS:
//     this->saveCMatrix(RUN_RKS, 0, C);
//     break;
//   case METHOD_UKS:
//     this->saveCMatrix(RUN_UKS_ALPHA, 0, C);
//     this->saveCMatrix(RUN_UKS_BETA, 0, C);
//     break;
//   case METHOD_ROKS:
//     CnErr.abort("to implement. sorry.");
//     break;
//   default:
//     CnErr.abort("unknown SCF type. stop.");
//     break;
//   }

//   return C;
// }

// // output guess in orthonormal basis to a files in fl_Work directory
// void DfInitialGuessHuckel::generateCprime(const TlMatrix& C){
//   TlMatrix Xinv;
//   Xinv.load("fl_Work/fl_Mtr_InvX.matrix");

//   Xinv *= C;

//   switch (this->m_nScfType){
//   case RKS:
//     Xinv.save("fl_Work/fl_Mtr_Cprime.matrix.rks0");
//     break;
//   default:
//     CnErr.abort("unknown SCF type. stop.");
//   }
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

void DfInitialGuessHuckel::generatePMatrix(const TlMatrix& C, const TlVector& occ,
                                           TlSymmetricMatrix& P1, TlSymmetricMatrix& P2)
{
    // generate density matrix
    const int nNumOfRows = C.getNumOfRows();
    const int nNumOfCols = C.getNumOfCols();

    P1 = TlSymmetricMatrix(nNumOfRows);
    P2 = TlSymmetricMatrix(nNumOfRows);

    for (int p = 0; p < nNumOfRows; p++) {
        for (int q = 0; q < nNumOfCols; q++) {
            double p1 = 0.0;
            double p2 = 0.0;

            for (int i = 0; i < occ.getSize(); i++) {
                double occNum = occ[i];

                if (std::fabs(occNum - 2.0) < EPS) {
                    p1 += C(p, i) * C(q, i);
                } else if (std::fabs(occNum - 1.0) < EPS) {
                    p2 += C(p, i) * C(q, i);
                }
            }

            P1(p, q) = p1;
            P2(p, q) = p2;
        }
    }
}


