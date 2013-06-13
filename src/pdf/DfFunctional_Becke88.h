#ifndef DFFUNCTIONAL_BECKE88_H
#define DFFUNCTIONAL_BECKE88_H

#include "DfFunctional.h"
#include "TlMatrix.h"

class DfFunctional_Becke88 : public DfFunctional_GGA {
public:
    DfFunctional_Becke88();
    virtual ~DfFunctional_Becke88();

public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    virtual double getFunctional(double dRhoA, double dGammaAA);

    virtual void getDerivativeFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);


public:
    // for Grid-Free ===================================================
    virtual TlMatrix getFunctionalCore(const double rhoA, 
                                       const double rhoB,
                                       const double xA,
                                       const double xB);

    virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB);

protected:
    double x(double dRho, double dGamma);

    virtual double g(double x);
    virtual double g_prime(double x);

protected:
    static const double EPS;
    static const double INV_3;
    static const double M_4_3;
    static const double G_PARAM;
    static const double BECKE_B;

private:
    friend class DfFunctional_B88LYP;
};

#endif // DFFUNCTIONAL_BECKE88_H
