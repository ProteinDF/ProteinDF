#include <cmath>
#include <cassert>
#include "DfFunctional_VWN3.h"

const double DfFunctional_VWN3::VWN3_A_PARA  =   0.0310907;
const double DfFunctional_VWN3::VWN3_B_PARA  =  13.0720;
const double DfFunctional_VWN3::VWN3_C_PARA  =  42.7198;
const double DfFunctional_VWN3::VWN3_X0_PARA = - 0.409286;
const double DfFunctional_VWN3::VWN3_A_FERR  =   0.01554535;
const double DfFunctional_VWN3::VWN3_B_FERR  =  20.1231;
const double DfFunctional_VWN3::VWN3_C_FERR  = 101.578;
const double DfFunctional_VWN3::VWN3_X0_FERR = - 0.743294;

DfFunctional_VWN3::DfFunctional_VWN3()
{
}

DfFunctional_VWN3::~DfFunctional_VWN3()
{
}

double DfFunctional_VWN3::epsilonC_PARA(const double x)
{
    return this->epsilonC(VWN3_A_PARA, VWN3_B_PARA, VWN3_C_PARA, VWN3_X0_PARA, x);
}

double DfFunctional_VWN3::epsilonC_FERR(const double x)
{
    return this->epsilonC(VWN3_A_FERR, VWN3_B_FERR, VWN3_C_FERR, VWN3_X0_FERR, x);
}

double DfFunctional_VWN3::epsilonCPrime_PARA(const double x)
{
    return this->epsilonCPrime(VWN3_A_PARA, VWN3_B_PARA, VWN3_C_PARA, VWN3_X0_PARA, x);
}

double DfFunctional_VWN3::epsilonCPrime_FERR(const double x)
{
    return this->epsilonCPrime(VWN3_A_FERR, VWN3_B_FERR, VWN3_C_FERR, VWN3_X0_FERR, x);
}


