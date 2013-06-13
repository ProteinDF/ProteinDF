#ifndef DFFUNCTIONAL_B3LYP_H
#define DFFUNCTIONAL_B3LYP_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"
#include "DfFunctional_Becke88_ExceptLDA.h"
#include "DfFunctional_VWN3.h"
#include "DfFunctional_LYP.h"

class DfFunctional_B3LYP : public DfFunctional_GGA {
public:
    DfFunctional_B3LYP();
    virtual ~DfFunctional_B3LYP();

public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    virtual double getFunctional(double dRhoA, double dGammaAA);

    virtual void getDerivativeFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);

protected:
    virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB,
                                       const double xA, const double xB);
    virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB);

protected:
    DfFunctional_Slater m_LDA;
    DfFunctional_Becke88_ExceptLDA m_B88;
    DfFunctional_VWN3 m_VWN;
    DfFunctional_LYP m_LYP;

private:
    static const double LDA_COEF;
    static const double B88_COEF;
    static const double VWN_COEF;
    static const double LYP_COEF;
};

#endif // DFFUNCTIONAL_B3LYP_H
