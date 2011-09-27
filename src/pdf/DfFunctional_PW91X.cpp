#include <cmath>
#include "DfFunctional_PW91X.h"

#include "TlUtils.h"
#include "TlMath.h"

DfFunctional_PW91X::DfFunctional_PW91X()
{
}

DfFunctional_PW91X::~DfFunctional_PW91X()
{
}

double DfFunctional_PW91X::getFunctional(const double dRhoA, const double dGammaAA)
{
    return (2.0 * this->pw91x(dRhoA, dGammaAA));
}

double DfFunctional_PW91X::getFunctional(const double dRhoA, const double dRhoB,
                                         const double dGammaAA, const double dGammaAB, const double dGammaBB)
{
    return (this->pw91x(dRhoA, dGammaAA) + this->pw91x(dRhoB, dGammaBB));
}

void DfFunctional_PW91X::getDerivativeFunctional(const double dRhoA, const double dGammaAA,
                                                 double* pRoundF_roundRhoA,
                                                 double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB)
{
    *pRoundF_roundRhoA = this->pw91x_roundRho(dRhoA, dGammaAA);
    *pRoundF_roundGammaAA = this->pw91x_roundGamma(dRhoA, dGammaAA);
    *pRoundF_roundGammaAB = 0.0;
}

void DfFunctional_PW91X::getDerivativeFunctional(const double dRhoA, const double dRhoB,
                                                 const double dGammaAA, const double dGammaAB, const double dGammaBB,
                                                 double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                                 double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB)
{
    *pRoundF_roundRhoA = this->pw91x_roundRho(dRhoA, dGammaAA);
    *pRoundF_roundRhoB = this->pw91x_roundRho(dRhoB, dGammaBB);

    *pRoundF_roundGammaAA = this->pw91x_roundGamma(dRhoA, dGammaAA);
    *pRoundF_roundGammaAB = 0.0;
    *pRoundF_roundGammaBB = this->pw91x_roundGamma(dRhoB, dGammaBB);
}

double DfFunctional_PW91X::pw91x(const double r, const double s)
{
    static const double M_1_3 = 1.0 / 3.0;
    static const double A = pow(3.0*M_PI*M_PI, M_1_3);
    static const double B = -3.0 / (4.0 * M_PI);
    static const double C = 0.5;
    static const double D = 0.19645;
    static const double E = 7.7956;
    static const double F = 0.2743;
    static const double G = 0.1508;
    static const double H = 0.004;

    static const double M_2_1_3 = pow(2.0,M_1_3);
    static const double M_2_2_3 = M_2_1_3 * M_2_1_3;

    const double r_4_3 = std::pow(r, 4.0/3.0);
    const double r_8_3 = r_4_3 * r_4_3;
    const double r_16_3 = r_8_3 * r_8_3;
    const double s_1_2_C = sqrt(s) * C;
    const double M_213_Ar = M_2_1_3 * A * r_4_3;
    const double Dasinh = D*TlMath::arcsinh(s_1_2_C*E / M_213_Ar);

    static const double AA = A*A;
    static const double AAAA = AA*AA;
    static const double CC = C*C;
    static const double CCCC = CC*CC;
    const double sCC = s*CC;

    const double term1 = M_2_1_3*A*B * r_4_3;
    const double term2 = sCC*(F-exp(-100.0*sCC/(M_2_2_3*r_8_3*AA))*G) / (M_2_2_3*r_8_3*AA);
    const double term3 = s_1_2_C * Dasinh / M_213_Ar;
    const double term4 = s*s*CCCC*H/(2.0 * M_2_1_3 * r_16_3 * AAAA);

    double dAnswer = term1 * (term2 + term3 +1.0) / (term4 + term3 +1.0);

    return dAnswer;
}

double DfFunctional_PW91X::pw91x_roundRho(const double r, const double s)
{
    static const double M_1_3 = 1.0 / 3.0;
    static const double A = pow(3.0*M_PI*M_PI, M_1_3);
    static const double B = -3.0 / (4.0 * M_PI);
    static const double C = 0.5;
    static const double D = 0.19645;
    static const double E = 7.7956;
    static const double F = 0.2743;
    static const double G = 0.1508;
    static const double H = 0.004;

    static const double M_2_1_3 = pow(2.0, M_1_3);
    static const double M_2_2_3 = M_2_1_3 * M_2_1_3;

    const double r_1_3 = std::pow(r, M_1_3);
    const double r_4_3 = std::pow(r, 4.0/3.0);
    const double r_8_3 = r_4_3 * r_4_3;
    const double r_11_3 = r_8_3 * r;
    const double r_16_3 = r_8_3 * r_8_3;
    const double r_19_3 = r_16_3 * r;
    const double s_1_2_C = sqrt(s) * C;
    const double M_213_Ar = M_2_1_3 * A * r_4_3;
    const double Dasinh = D*TlMath::arcsinh(s_1_2_C*E / M_213_Ar);

    static const double AA = A*A;
    static const double AAAA = AA*AA;
    static const double AB = A*B;
    static const double CC = C*C;
    static const double CCCC = CC*CC;
    const double sCC = s*CC;
    const double sCCDE4 = 4.0 * sCC * D * E;

    const double term1 = 4.0*M_2_1_3*AB*r_1_3;
    const double term2_e = exp(-100.0*sCC/(M_2_2_3*r_8_3*AA))*G;
    const double term2_1 = F - term2_e;
    const double term2 = sCC * term2_1 / (M_2_2_3*r_8_3*AA);
    const double term3 = s_1_2_C * Dasinh / M_213_Ar;
    const double term4 = s*s*CCCC*H/(2.0 * M_2_1_3 * r_16_3 * AAAA);

    const double term5 = M_2_1_3*AB*r_4_3;
    const double term6 = 8.0 * sCC * term2_1 / (3.0 * M_2_2_3 * r_11_3 * AA);
    const double term7 = 400.0 * sCC * sCC * term2_e / (3.0 * M_2_1_3 * r_19_3 * AAAA);
    const double term8 = term3 * (4.0 / (3.0 * r));
    const double term9_1 = sCC*E*E / (M_2_2_3 * r_8_3 * AA) +1.0;
    const double term9 = sCCDE4 / (3.0 * M_2_2_3 * r_11_3 * AA * sqrt(term9_1));

    const double term10_1 = sCC*sCC*H / (r_16_3* AAAA);
    const double term10 = term10_1 / (2.0 * M_2_1_3);

    const double term11 = (16.0 * pow(2.0, -4.0/3.0) / 3.0) * term10_1 / r;
    const double term12 = (4.0 / (3.0 * r)) * term3;

    const double dAnswer1_1 = term2 + term3 +1.0;
    const double dAnswer1 = term1 * dAnswer1_1 /(3.0 * (term4 + term3 +1.0));
    const double dAnswer2_1 = term5 * (-term6 -term7 -term8 -term9);
    const double dAnswer2_2 = term10 + term3 +1.0;
    const double dAnswer2 = dAnswer2_1 / dAnswer2_2;
    const double dAnswer3_1 = term5 * dAnswer1_1 * (-term11 -term12 -term9);
    const double dAnswer3_2 = term4 + term3 +1.0;
    const double dAnswer3 = dAnswer3_1 / (dAnswer3_2 * dAnswer3_2);

    return (dAnswer1 + dAnswer2 - dAnswer3);
}

double DfFunctional_PW91X::pw91x_roundGamma(const double r, const double s)
{
    static const double M_1_3 = 1.0 / 3.0;
    static const double A = pow(3.0*M_PI*M_PI, M_1_3);
    static const double B = -3.0 / (4.0 * M_PI);
    static const double C = 0.5;
    static const double D = 0.19645;
    static const double E = 7.7956;
    static const double F = 0.2743;
    static const double G = 0.1508;
    static const double H = 0.004;

    static const double M_2_1_3 = pow(2.0, M_1_3);
    static const double M_2_2_3 = M_2_1_3 * M_2_1_3;

    //const double r_1_3 = std::pow(r, M_1_3);
    const double r_4_3 = std::pow(r, 4.0/3.0);
    const double r_8_3 = r_4_3 * r_4_3;
    //const double r_11_3 = r_8_3 * r;
    const double r_16_3 = r_8_3 * r_8_3;
    //const double r_19_3 = r_16_3 * r;
    const double s_1_2_C = sqrt(s) * C;
    const double M_213_Ar = M_2_1_3 * A * r_4_3;
    const double Dasinh = D*TlMath::arcsinh(s_1_2_C*E / M_213_Ar);

    static const double AA = A*A;
    static const double AAAA = AA*AA;
    static const double AB = A*B;
    static const double CC = C*C;
    //static const double CCCC = CC*CC;
    const double sCC = s*CC;
    //const double sCCDE4 = 4.0 * sCC * D * E;

    const double term1 = M_2_1_3 * r_4_3 * AB;
    const double term2_b = M_2_2_3*r_8_3*AA;
    const double term2_e = exp(-100.0*sCC/term2_b)*G;
    const double term2_1 = F - term2_e;
    const double term2 = CC * term2_1 / term2_b;
    const double term3 = 50.0*sCC*CC*term2_e/(M_2_1_3*r_16_3*AAAA);
    const double term4 = C*Dasinh/(2.0*M_2_1_3*r_4_3*sqrt(s)*A);
    const double term5 = CC*D*E/(2.0*M_2_2_3*r_8_3*AA*sqrt(sCC*E*E/(M_2_2_3*r_8_3*AA)+1.0));

    const double term6 = sCC*sCC*H/(2.0*M_2_1_3*r_16_3*AAAA);
    const double term7 = sqrt(s)*C*D*Dasinh/(M_2_1_3*r_4_3*A);

    const double ans1 = term1*(term2 + term3 + term4 + term5) / (term6 + term7 +1.0);

    const double term8 = M_2_1_3*r_4_3*AB;
    const double term9 = sCC*(F-term2_e)/term2_b;
    const double term10 = sqrt(s)*C*Dasinh/(M_2_1_3*r_4_3*A);

    const double term11 = 2.0*pow(2.0,-4.0/3.0)*sCC*CC*H/(r_16_3*AAAA);
    const double term12 = C*Dasinh/(2.0*M_2_1_3*r_4_3*sqrt(s)*A);
    const double term13 = CC*D*E/(2.0*M_2_2_3*r_8_3*AA*sqrt(sCC*E*E/(M_2_2_3*r_8_3*AA)+1.0));

    const double term14 = sCC*sCC*H/(2.0*M_2_1_3*r_16_3*AAAA);
    const double term15 = sqrt(s)*C*Dasinh/(M_2_1_3*r_4_3*A);
    const double ans2_b = term14 + term15 + 1.0;

    const double ans2 = term8*(term9 + term10 +1.0)*(term11+term12+term13)/(ans2_b*ans2_b);

    const double ans = ans1 - ans2;

    return ans;
}



