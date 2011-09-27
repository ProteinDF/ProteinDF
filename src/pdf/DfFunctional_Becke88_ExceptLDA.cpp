#include "DfFunctional_Becke88_ExceptLDA.h"
#include "TlMath.h"

DfFunctional_Becke88_ExceptLDA::DfFunctional_Becke88_ExceptLDA()
{
}

DfFunctional_Becke88_ExceptLDA::~DfFunctional_Becke88_ExceptLDA()
{
}

double DfFunctional_Becke88_ExceptLDA::g(const double x)
{
    const double bx = BECKE_B * x;
    const double arcsinhx = TlMath::arcsinh(x);

    const double dAnswer = -((bx * x) / (1.0 + 6.0 * bx * arcsinhx));

    return dAnswer;
}

