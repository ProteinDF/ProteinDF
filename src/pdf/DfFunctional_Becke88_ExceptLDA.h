#ifndef DFFUNCTIONAL_BECKE88_EXCEPTLDA_H
#define DFFUNCTIONAL_BECKE88_EXCEPTLDA_H

#include "DfFunctional_Becke88.h"

class DfFunctional_Becke88_ExceptLDA : public DfFunctional_Becke88 {
public:
    DfFunctional_Becke88_ExceptLDA();
    virtual ~DfFunctional_Becke88_ExceptLDA();

protected:
    virtual double g(double x);
    //virtual double g_prime(double x);
};

#endif // DFFUNCTIONAL_BECKE88_EXCEPTLDA_H
