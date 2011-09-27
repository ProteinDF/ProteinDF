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

