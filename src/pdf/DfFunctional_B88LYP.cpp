#include "DfFunctional_B88LYP.h"
#include "TlUtils.h"

DfFunctional_B88LYP::DfFunctional_B88LYP()
{
    this->numOfFunctionalTerms_ = 
        this->m_Becke88.getNumOfFunctionalTerms() + this->m_LYP.getNumOfFunctionalTerms();
    this->numOfDerivativeFunctionalTerms_ = 
        this->m_Becke88.getNumOfDerivativeFunctionalTerms() + this->m_LYP.getNumOfDerivativeFunctionalTerms();
}

DfFunctional_B88LYP::~DfFunctional_B88LYP()
{
}

double DfFunctional_B88LYP::getFunctional(const double dRhoA, const double dRhoB,
                                          const double dGammaAA, const double dGammaAB, const double dGammaBB)
{
    const double ex = this->m_Becke88.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);
    const double ec = this->m_LYP.getFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB);

    return ex + ec;
}

void DfFunctional_B88LYP::getDerivativeFunctional(const double dRhoA, const double dRhoB,
                                                  const double dGammaAA, const double dGammaAB, const double dGammaBB,
                                                  double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                                  double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB)
{
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);
    assert(pRoundF_roundGammaAA != NULL);
    assert(pRoundF_roundGammaAB != NULL);
    assert(pRoundF_roundGammaBB != NULL);

    double dRoundF_roundRhoA_X, dRoundF_roundRhoB_X;
    double dRoundF_roundRhoA_C, dRoundF_roundRhoB_C;
    double dRoundF_roundGammaAA_X, dRoundF_roundGammaAB_X, dRoundF_roundGammaBB_X;
    double dRoundF_roundGammaAA_C, dRoundF_roundGammaAB_C, dRoundF_roundGammaBB_C;

    this->m_Becke88.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                                            &dRoundF_roundRhoA_X, &dRoundF_roundRhoB_X,
                                            &dRoundF_roundGammaAA_X, &dRoundF_roundGammaAB_X, &dRoundF_roundGammaBB_X);
    this->m_LYP.getDerivativeFunctional(dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB,
                                        &dRoundF_roundRhoA_C, &dRoundF_roundRhoB_C,
                                        &dRoundF_roundGammaAA_C, &dRoundF_roundGammaAB_C, &dRoundF_roundGammaBB_C);

    *pRoundF_roundRhoA = dRoundF_roundRhoA_X + dRoundF_roundRhoA_C;
    *pRoundF_roundRhoB = dRoundF_roundRhoB_X + dRoundF_roundRhoB_C;
    *pRoundF_roundGammaAA = dRoundF_roundGammaAA_X + dRoundF_roundGammaAA_C;
    *pRoundF_roundGammaAB = dRoundF_roundGammaAB_X + dRoundF_roundGammaAB_C;
    *pRoundF_roundGammaBB = dRoundF_roundGammaBB_X + dRoundF_roundGammaBB_C;

//   std::cerr << TlUtils::format("RA=%5e, RB=%5e, GAA=%5e, GAB=%5e, GBB=%5e",
//                 dRhoA, dRhoB, dGammaAA, dGammaAB, dGammaBB)
//      << std::endl;
//   std::cerr << TlUtils::format("X: F/A=%5e, F/B=%5e, F/AA=%5e, F/AB=%5e, F/BB=%5e",
//                 dRoundF_roundRhoA_X, dRoundF_roundRhoB_X,
//                 dRoundF_roundGammaAA_X, dRoundF_roundGammaAB_X, dRoundF_roundGammaBB_X)
//      << std::endl;
//   std::cerr << TlUtils::format("C: F/A=%5e, F/B=%5e, F/AA=%5e, F/AB=%5e, F/BB=%5e",
//                 dRoundF_roundRhoA_C, dRoundF_roundRhoB_C,
//                 dRoundF_roundGammaAA_C, dRoundF_roundGammaAB_C, dRoundF_roundGammaBB_C)
//      << std::endl;
//   std::cerr << std::endl;
}

TlVector DfFunctional_B88LYP::getFunctionalTermCoef_GF()
{
    const int dim = this->getNumOfFunctionalTerms();
    TlVector coef(dim);
    for (int i = 0; i < dim; ++i) {
        coef[i] = 1.0;
    }

    return coef;
}

TlVector DfFunctional_B88LYP::getDerivativeFunctionalTermCoef_GF()
{
    const int dim = this->getNumOfDerivativeFunctionalTerms();
    TlVector coef(dim);
    for (int i = 0; i < dim; ++i) {
        coef[i] = 1.0;
    }

    return coef;
}


TlMatrix DfFunctional_B88LYP::getFunctionalCore(const double rhoA, 
                                                const double rhoB,
                                                const double xA,
                                                const double xB)
{
    const TlMatrix x = this->m_Becke88.getFunctionalCore(rhoA, rhoB, xA, xB);
    const TlMatrix c = this->m_LYP.getFunctionalCore(rhoA, rhoB, xA, xB);

    const index_type numOfFunctionalTerms_X = this->m_Becke88.getNumOfFunctionalTerms();
    const index_type numOfFunctionalTerms_C = this->m_LYP.getNumOfFunctionalTerms();
    TlMatrix xc(F_DIM, this->getNumOfFunctionalTerms());
    for (index_type i = 0; i < F_DIM; ++i) {
        for (index_type j = 0; j < numOfFunctionalTerms_X; ++j) {
            xc(i, j) = x(i, j);
        }
        for (index_type j = 0; j < numOfFunctionalTerms_C; ++j) {
            xc(i, numOfFunctionalTerms_X + j) = c(i, j);
        }
    }

    return xc;
}

TlMatrix DfFunctional_B88LYP::getDerivativeFunctionalCore(const double rhoA,
                                                          const double rhoB,
                                                          const double xA,
                                                          const double xB)
{
    const TlMatrix x = this->m_Becke88.getDerivativeFunctionalCore(rhoA, rhoB, xA, xB);
    const TlMatrix c = this->m_LYP.getDerivativeFunctionalCore(rhoA, rhoB, xA, xB);

    const index_type numOfDerivativeFunctionalTerms_X = this->m_Becke88.getNumOfDerivativeFunctionalTerms();
    const index_type numOfDerivativeFunctionalTerms_C = this->m_LYP.getNumOfDerivativeFunctionalTerms();
    TlMatrix xc(D_DIM, this->getNumOfDerivativeFunctionalTerms());
    for (index_type i = 0; i < D_DIM; ++i) {
        for (index_type j = 0; j < numOfDerivativeFunctionalTerms_X; ++j) {
            xc(i, j) = x(i, j);
        }
        for (index_type j = 0; j < numOfDerivativeFunctionalTerms_C; ++j) {
            xc(i, numOfDerivativeFunctionalTerms_X + j) = c(i, j);
        }
    }

    return xc;
}

