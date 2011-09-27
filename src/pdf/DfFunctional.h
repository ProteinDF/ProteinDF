#ifndef DFFUNCTIONAL_H
#define DFFUNCTIONAL_H

#include <vector>
#include "TlSymmetricMatrix.h"

/** 局所汎関数を扱うインターフェースクラス
 *
 */
class DfFunctional_LDA {
public:
    DfFunctional_LDA();
    virtual ~DfFunctional_LDA();

public:
    /** 交換相関関数値を返す
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dRhoB beta電子密度
     */
    virtual double getFunctional(double dRhoA, double dRhoB) = 0;

    /** 交換相関関数値を返す(RKSに特殊化)
     *
     *  @param[in] dRhoA alpha電子密度
     */
    double getFunctional(double dRhoA) {
        return this->getFunctional(dRhoA, dRhoA);
    }

    /** 交換相関関数の一次微分を返す
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dRhoB beta電子密度
     *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
     *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
     */
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB) =0;

    /** 交換相関関数の一次微分を返す(RKSに特殊化)
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dRhoB beta電子密度
     *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
     *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
     */
    void getDerivativeFunctional(double dRhoA, double* pRoundF_roundRhoA) {
        double dDummyVariable;
        this->getDerivativeFunctional(dRhoA, dRhoA, pRoundF_roundRhoA, &dDummyVariable);
    }
};

/** 非局所汎関数を扱うインターフェースクラス
 *
 */
class DfFunctional_GGA {
public:
    DfFunctional_GGA();
    virtual ~DfFunctional_GGA();

public:
    /** 交換相関関数値を返す
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dRhoB beta電子密度
     *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
     *  @param[in] dGammaAB (grad Rho alpha) * (grad Rho beta)
     *  @param[in] dGammaBB (grad Rho beta) * (grad Rho beta)
     */
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB) =0;

    /** 交換相関関数値を返す(RKSに特殊化)
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
     */
    double getFunctional(double dRhoA, double dGammaAA) {
        return this->getFunctional(dRhoA, dRhoA, dGammaAA, dGammaAA, dGammaAA);
    }

    /** 交換相関関数の一次微分を返す
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dRhoB beta電子密度
     *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
     *  @param[in] dGammaAB (grad Rho alpha) * (grad Rho beta)
     *  @param[in] dGammaBB (grad Rho beta) * (grad Rho beta)
     *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
     *  @param[out] pRoundF_roundRhoB beta電子密度による一次偏微分
     *  @param[out] pRoundF_roundGammaAA (gamma alpha alpha)による一次偏微分
     *  @param[out] pRoundF_roundGammaAB (gamma alpha beta)による一次偏微分
     *  @param[out] pRoundF_roundGammaBB (gamma beta beta)による一次偏微分
     */
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB) =0;

    /** 交換相関関数の一次微分を返す(RKSに特殊化)
     *
     *  @param[in] dRhoA alpha電子密度
     *  @param[in] dGammaAA (grad Rho alpha) * (grad Rho alpha)
     *  @param[out] pRoundF_roundRhoA alpha電子密度による一次偏微分
     *  @param[out] pRoundF_roundGammaAA (gamma alpha alpha)による一次偏微分
     *  @param[out] pRoundF_roundGammaAB (gamma alpha beta)による一次偏微分
     */
    void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                 double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB) {
        double dRoundF_roundRhoB;
        double dRoundF_roundGammaBB;

        this->getDerivativeFunctional(dRhoA, dRhoA, dGammaAA, dGammaAA, dGammaAA,
                                      pRoundF_roundRhoA, &dRoundF_roundRhoB,
                                      pRoundF_roundGammaAA, pRoundF_roundGammaAB, &dRoundF_roundGammaBB);
    }
};

#endif // DFFUNCTIONAL_H
