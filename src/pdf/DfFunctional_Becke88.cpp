#include <cmath>
#include <limits>
#include "DfFunctional_Becke88.h"
#include "TlUtils.h"
#include "TlMath.h"

const double DfFunctional_Becke88::EPS = 1.0E-16;//std::numeric_limits<double>::epsilon();
const double DfFunctional_Becke88::INV_3 = 1.0 / 3.0;
const double DfFunctional_Becke88::M_4_3 = 4.0 / 3.0;
const double DfFunctional_Becke88::G_PARAM = - (3.0 / 2.0) * std::pow(3.0 / (4.0 * M_PI), (1.0 / 3.0));
const double DfFunctional_Becke88::BECKE_B = 0.0042;

DfFunctional_Becke88::DfFunctional_Becke88()
{
    this->numOfFunctionalTerms_ = 1;
    this->numOfDerivativeFunctionalTerms_ = 1;
}

DfFunctional_Becke88::~DfFunctional_Becke88()
{
}

double DfFunctional_Becke88::getFunctional(const double dRhoA, const double dRhoB,
                                           const double dGammaAA, const double dGammaAB, const double dGammaBB)
{
    double termA = 0.0;
    if (dRhoA > EPS) {
        const double xA = this->x(dRhoA, dGammaAA);
        termA = std::pow(dRhoA, M_4_3) * this->g(xA);
    }
    double termB = 0.0;
    if (dRhoB > EPS) {
        const double xB = this->x(dRhoB, dGammaBB);
        termB = std::pow(dRhoB, M_4_3) * this->g(xB);
    }

    return (termA + termB);
}

// for RKS
double DfFunctional_Becke88::getFunctional(const double dRhoA, const double dGammaAA)
{
    double termA = 0.0;
    if (dRhoA > EPS) {
        const double xA = this->x(dRhoA, dGammaAA);
        termA = std::pow(dRhoA, M_4_3) * this->g(xA);
    }

    return (2.0 * termA);
}

void DfFunctional_Becke88::getDerivativeFunctional(const double dRhoA, const double dRhoB,
                                                   const double dGammaAA, const double dGammaAB, const double dGammaBB,
                                                   double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                                   double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB)
{
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundRhoB != NULL);
    assert(pRoundF_roundGammaAA != NULL);
    assert(pRoundF_roundGammaAB != NULL);
    assert(pRoundF_roundGammaBB != NULL);

    *pRoundF_roundGammaAB = 0.0;

    *pRoundF_roundRhoA = 0.0;
    *pRoundF_roundGammaAA = 0.0;
    if ((dRhoA > EPS) && (dGammaAA > EPS)) {
        const double xA = this->x(dRhoA, dGammaAA);
        const double gA = this->g(xA); //G_PARAM + this->getBecke88Coef(xA);
        const double g_primeA = this->g_prime(xA);

        // roundF_roundRho
        *pRoundF_roundRhoA = M_4_3 * std::pow(dRhoA, INV_3) * (gA - xA * g_primeA);
        // roundF_roundGamma
        *pRoundF_roundGammaAA = 0.5 * (1.0 / std::sqrt(dGammaAA)) * g_primeA;
    }

    *pRoundF_roundRhoB = 0.0;
    *pRoundF_roundGammaBB = 0.0;
    if ((dRhoB > EPS) && (dGammaBB > EPS)) {
        const double xB = (dRhoB > EPS) ? this->x(dRhoB, dGammaBB) : 1.0;
        const double gB = this->g(xB); //G_PARAM + this->getBecke88Coef(xB);
        const double g_primeB = this->g_prime(xB);
        *pRoundF_roundRhoB = M_4_3 * std::pow(dRhoB, INV_3) * (gB - xB * g_primeB);
        *pRoundF_roundGammaBB = 0.5 * (1.0 / std::sqrt(dGammaBB)) * g_primeB;
    }
}

// for RKS
void DfFunctional_Becke88::getDerivativeFunctional(const double dRhoA, const double dGammaAA,
                                                   double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB)
{
    assert(pRoundF_roundRhoA != NULL);
    assert(pRoundF_roundGammaAA != NULL);
    assert(pRoundF_roundGammaAB != NULL);

    *pRoundF_roundGammaAB = 0.0;

    *pRoundF_roundRhoA = 0.0;
    *pRoundF_roundGammaAA = 0.0;

    if ((dRhoA > EPS) && (dGammaAA > EPS)) {
        const double xA = this->x(dRhoA, dGammaAA);
        const double gA = this->g(xA);
        const double g_primeA = this->g_prime(xA);

        // roundF_roundRho
        *pRoundF_roundRhoA = M_4_3 * std::pow(dRhoA, INV_3) * (gA - xA * g_primeA);

        // roundF_roundGamma
        *pRoundF_roundGammaAA = 0.5 * (1.0 / std::sqrt(dGammaAA)) * g_primeA;
    }
}

double DfFunctional_Becke88::x(const double dRhoA, const double dGammaAA)
{
    return std::sqrt(dGammaAA) / std::pow(dRhoA, M_4_3);
}

double DfFunctional_Becke88::g(const double x)
{
    const double bx = BECKE_B * x;
    const double arcsinhx = TlMath::arcsinh(x);

    const double dAnswer = G_PARAM - ((bx * x) / (1.0 + 6.0 * bx * arcsinhx));
    return dAnswer;
}

double DfFunctional_Becke88::g_prime(const double x)
{
    const double bx = BECKE_B * x;
    const double xx = x * x;
    const double arcsinhx = TlMath::arcsinh(x);
    const double M = 1.0 + 6.0 * bx * arcsinhx;
    const double invM = 1.0 / M;

    const double arg1 = 6.0 * bx * bx * (arcsinhx + x / std::sqrt(xx +1.0)) * invM *invM;
    const double arg2 = - 2.0 * bx * invM;

    const double dAnswer = arg1 + arg2;

    return dAnswer;
}

// ----------
TlMatrix DfFunctional_Becke88::getFunctionalCore(const double rhoA, 
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB)
{
    TlMatrix answer(F_DIM, this->getNumOfFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    // termA = std::pow(rhoA, M_4_3) * this->g(xA);
    // termB = std::pow(rhoB, M_4_3) * this->g(xB);
    // E = termA + termB;

    const double FA_termR = std::pow(rhoA, M_4_3);
    const double FA_termX = this->g(xA);

    const double FB_termR = std::pow(rhoB, M_4_3);
    const double FB_termX = this->g(xB);

    answer.set(FA_R, 0, FA_termR);
    answer.set(FA_X, 0, FA_termX);
    answer.set(FB_R, 0, FB_termR);
    answer.set(FB_X, 0, FB_termX);

    return answer;
}

TlMatrix DfFunctional_Becke88::getDerivativeFunctionalCore(const double rhoA,
                                                           const double rhoB,
                                                           const double xA,
                                                           const double xB)
{
    TlMatrix answer(D_DIM, this->getNumOfDerivativeFunctionalTerms());
    assert(this->getNumOfFunctionalTerms() == 1);

    const double gA = this->g(xA);
    const double gB = this->g(xA);
    const double g_primeA = this->g_prime(xA);
    const double g_primeB = this->g_prime(xB);
    
    // roundF_roundRho =========================================================
    // roundF_roundRhoA = M_4_3 * std::pow(rhoA, INV_3) * (gA - xA * g_primeA);
    //                  = [M_4_3 * std::pow(rhoA, INV_3)] * [(gA - xA * g_primeA)]
    const double roundF_roundRhoA_termR = M_4_3 * std::pow(rhoA, INV_3);
    const double roundF_roundRhoA_termX = gA - xA * g_primeA;
    const double roundF_roundRhoB_termR = M_4_3 * std::pow(rhoB, INV_3);
    const double roundF_roundRhoB_termX = gB - xB * g_primeB;
    answer.set(RA_R, 0, roundF_roundRhoA_termR);
    answer.set(RA_X, 0, roundF_roundRhoA_termX);
    answer.set(RB_R, 0, roundF_roundRhoB_termR);
    answer.set(RB_X, 0, roundF_roundRhoB_termX);
    
    // roundF_roundGamma =======================================================
    // roundF_roundGammaAA = 0.5 * (1.0 / std::sqrt(gammaAA)) * g_primeA;
    // gammaAA = (rhoA^(4/3) * xA) * (rhoA^(4/3) * xA);
    // roundF_roundGammaAA = 0.5 * (1.0 / (rhoA^(4/3) * xA)) * g_primeA;
    //                     = [1.0 / rhoA^(4/3)] * [0.5 * (1.0 / xA) * g_primeA]
    if ((rhoA > 1.0E-16) && (xA > 1.0E-16)) {
        const double g_primeA = this->g_prime(xA);
        const double inv_xA = 1.0 / xA;

        const double roundF_roundGammaAA_termR = 1.0 / std::pow(rhoA, 4.0/3.0);
        const double roundF_roundGammaAA_termX = 0.5 * inv_xA *  g_primeA;
        answer.set(GAA_R, 0, roundF_roundGammaAA_termR);
        answer.set(GAA_X, 0, roundF_roundGammaAA_termX);
    }
    {
        const double roundF_roundGammaAB_termR = 0.0;
        const double roundF_roundGammaAB_termX = 0.0;
        answer.set(GAB_R, 0, roundF_roundGammaAB_termR);
        answer.set(GAB_X, 0, roundF_roundGammaAB_termX);
    }
    if ((rhoB > 1.0E-16) && (xB > 1.0E-16)) {
        const double g_primeB = this->g_prime(xB);
        const double inv_xB = 1.0 / xB;

        const double roundF_roundGammaBB_termR = 1.0 / std::pow(rhoB, 4.0/3.0);
        const double roundF_roundGammaBB_termX = 0.5 * inv_xB *  g_primeB;
        answer.set(GBB_R, 0, roundF_roundGammaBB_termR);
        answer.set(GBB_X, 0, roundF_roundGammaBB_termX);
    }

    return answer;
}
