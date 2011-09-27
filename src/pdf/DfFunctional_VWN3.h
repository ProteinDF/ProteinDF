#ifndef DFFUNCTIONAL_VWN3_H
#define DFFUNCTIONAL_VWN3_H

#include "DfFunctional_VWN.h"

// VWN3
class DfFunctional_VWN3 : public DfFunctional_VWN {
public:
    DfFunctional_VWN3();
    virtual ~DfFunctional_VWN3();

protected:
    virtual double epsilonC_PARA(double x);
    virtual double epsilonC_FERR(double x);

    virtual double epsilonCPrime_PARA(double x);
    virtual double epsilonCPrime_FERR(double x);

protected:
    static const double VWN3_A_PARA;
    static const double VWN3_B_PARA;
    static const double VWN3_C_PARA;
    static const double VWN3_X0_PARA;
    static const double VWN3_A_FERR;
    static const double VWN3_B_FERR;
    static const double VWN3_C_FERR;
    static const double VWN3_X0_FERR;
};

#endif // DFFUNCTIONAL_VWN3_H
