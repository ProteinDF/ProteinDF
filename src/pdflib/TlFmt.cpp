#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>

// STL
#include <list>
#include <numeric>

#include "TlFmt.h"

TlFmt* TlFmt::m_pInstance = NULL;

const double TlFmt::PI               = 3.14159265358979323846;
const double TlFmt::CK               = std::sqrt(2.0 * sqrt(TlFmt::PI)) * TlFmt::PI;
const double TlFmt::PI_DIV2          = TlFmt::PI / 2.0;
const double TlFmt::INV2             = 1.0 / 2.0;
const double TlFmt::INV6             = 1.0 / 6.0;
const int    TlFmt::FMT_M            = 14;
const int    TlFmt::FMT_INV_D        = 512;                                  // = pow(2.0, 9.0);
const int    TlFmt::FMT_MAX_M        = TlFmt::FMT_M + 3;                     // 4-term expansion
const double TlFmt::FMT_T            = 2 * TlFmt::FMT_M + 36;                // Tf = 2 * m + 36
const int    TlFmt::FMT_N_STEP       = (int)TlFmt::FMT_T * TlFmt::FMT_INV_D; // = Tf / d
const double TlFmt::FMT_STEPSIZE     = TlFmt::FMT_T / TlFmt::FMT_N_STEP;
const double TlFmt::INV_FMT_STEPSIZE = 1.0 / TlFmt::FMT_STEPSIZE;
const double TlFmt::EPSILON          = std::numeric_limits<double>::epsilon();    // double の計算機epsilon

TlFmt& TlFmt::getInstance()
{
    if (TlFmt::m_pInstance == NULL) {
        TlFmt::m_pInstance = new TlFmt();
    }

    return *TlFmt::m_pInstance;
}


TlFmt::TlFmt(int nFlag)
    : m_nFlags(nFlag), m_FmtTable(NULL), m_previousT(0.0), m_nNumOfFunctionCalled(0), m_nNumOfCacheHits(0)
{

    // for debug
    this->m_nFlags = TLFMT_ENABLE_CACHE;

    // initialize
    this->m_FmtTable = new double* [FMT_N_STEP +1];
    for (int n = 0; n < FMT_N_STEP +1; ++n) {
        this->m_FmtTable[n] = new double[FMT_MAX_M +1];
        for (int m = 0; m <= FMT_MAX_M; ++m) {
            this->m_FmtTable[n][m] = 0.0;
        }
    }

    this->m_aFmtCalcCache.clear();
    this->m_previousF.clear();

    // 表の作成
    this->makeFmtTable();
}


TlFmt::~TlFmt()
{
    if (this->m_FmtTable != NULL) {
        for (int n = 0; n < FMT_N_STEP +1; ++n) {
            if (this->m_FmtTable[n] != NULL) {
                delete[] this->m_FmtTable[n];
                this->m_FmtTable[n] = NULL;
            }
        }
        delete[] this->m_FmtTable;
        this->m_FmtTable = NULL;
    }
}


// this->m_FmtTable[][] の作成
//
// あらかじめ T = 0.0 〜 2m+36 まで、1/2^9 刻みでFm(T)を計算しておく。
// この値を基に、T = 0.0 〜 2m+36 の範囲では、Taylor展開を行う。
void TlFmt::makeFmtTable()
{
    // T=0の場合。
    // すなわち Fm(0) = 1 / (2*m +1)
    for (int m = 0; m <= FMT_MAX_M; ++m) {
        this->m_FmtTable[0][m] = 1.0 / (2*m +1);
    }

    // T>0 の場合
    int m = FMT_MAX_M; // m が小さいものは、後に帰納的に求める (*)
    const double THR_ZERO = 1.0E-17;

    for (int j = 1; j <= FMT_N_STEP; ++j) {
        const double t = FMT_STEPSIZE * j;
        const double expt = std::exp(-t);
        int nu = 2*m +1;
        const double t2 = t * 2.0;
        const double eps = (expt / t2) * THR_ZERO;
        double term = 1.0 / nu;
        // 変数 func に term を加算していくが、誤差を小さくするために
        // term が小さい順に加算する。後にsortする必要があるので、そのリスト。
        std::list<double> aFuncAdditionList(0); // 0 個のリスト
        aFuncAdditionList.push_back(term);       //func += term;
        int i = nu;

        do {
            i += 2;
            term *= (t2 / i);
            aFuncAdditionList.push_back(term);    //func += term;
        } while (term > eps);

        aFuncAdditionList.sort();
        double func = std::accumulate(aFuncAdditionList.begin(), aFuncAdditionList.end(), 0.0);

        this->m_FmtTable[j][m] = expt * func; // 0からTf(=2m+36)までをFMT_INV_D分割したうちのj番目のFm(T)の値 (T = FMT_STEPSIZE*j)

        // (*)
        // Fm(T) が求まれば、
        // F(m+1)(T) = -F'm(T) より、帰納的に求まる。
        for (int mm = m -1; mm >=0; mm--) {
            nu -= 2;
            this->m_FmtTable[j][mm] = (expt + t2 * this->m_FmtTable[j][mm +1]) / nu;
        }
    }
}

std::vector<double> TlFmt::getFmtSubstance(const int m, const double& t) const
{
    std::vector<double> f(m+1, 0.0);

    if (t <= (2*m + 36)) { // Tf = 2*m + 36
        const int ts       = static_cast<int>(0.5 + t * INV_FMT_STEPSIZE);
        const double delta = ts * FMT_STEPSIZE - t;
        for (int i=0; i<=m; i++) {
            f[i] = ((this->m_FmtTable[ts][i+3] * INV6 * delta
                     + this->m_FmtTable[ts][i+2] * INV2) * delta
                    + this->m_FmtTable[ts][i+1]) * delta + this->m_FmtTable[ts][i];
        }
    } else {
        const double t_inv = INV2 / t;
        f[0] = std::sqrt(PI_DIV2 * t_inv);
        double nu = 1.0;
        for (int i=1; i<=m; i++) {
            f[i] = t_inv * nu * f[i-1];
            nu += 2.0;
        }
    }

    // 実際に計算した回数をカウントアップ
    this->m_nNumOfFunctionCalled++;

    return f;
}

/// Fm(T)の値を返す
/// T が一定のとき、m が大きい方から呼び出した方が速いはず。
double TlFmt::getFmT0(const int m, const double& t) const
{
    assert(m >= 0);
    double answer = 0.0;

    if ((this->m_nFlags & TLFMT_ENABLE_CACHE) == 1) {
        // キャッシュを使う ================================================
        bool bHasCalced = false; // すでに計算されているかチェックするためのフラグ

        // キャッシュ内にすでに計算されているかチェック
        std::map<double, std::vector<double> >::const_iterator p = this->m_aFmtCalcCache.find(t);
        if (p != this->m_aFmtCalcCache.end()) {
            if (m < (int)(*p).second.size()) {
                answer = ((*p).second).at(m);

                // キャッシュにヒットした回数
                this->m_nNumOfCacheHits++;
                bHasCalced = true;
            }
        }

        if (bHasCalced == false) {
            // キャッシュにヒットしないので計算する
            std::vector<double> f = this->getFmtSubstance(m, t);
            answer = f[m];

            // 計算結果をキャッシュに入れておく
            (this->m_aFmtCalcCache)[t] = f;
        }
    } else {
        // キャッシュを使わない ============================================
//     if ((this->m_previousT == t) &&
//  (this->m_previousF.size() -1 > m)){
//       // 前回の値を呼ばれた
//       answer = this->m_previousF[m];
//     } else {
        // 新しいFを計算する
        std::vector<double> f = this->getFmtSubstance(m, t);
        answer = f[m];
        // 前回分を更新
        this->m_previousF = f;
        this->m_previousT = t;
//     }
    }

    // for debug
    //std::cout << "TlFmt::getFmT F[" << m << "](" << t << ") = " << answer << std::endl;

    return answer;
}


void TlFmt::getFmT(const int m, const double& t, double* pBuf) const
{
    if (t <= (2*m + 36)) { // Tf = 2*m + 36
        const int ts = static_cast<int>(0.5 + t * INV_FMT_STEPSIZE);
        const double delta = ts * FMT_STEPSIZE - t;
        for (int i = 0; i <= m; ++i) {
            pBuf[i] = this->m_FmtTable[ts][i]
                      + ((this->m_FmtTable[ts][i+3]*INV6*delta +this->m_FmtTable[ts][i+2]*INV2)*delta
                         +this->m_FmtTable[ts][i+1])*delta;
        }
    } else {
        const double t_inv = INV2 / t;
        pBuf[0] = std::sqrt(PI_DIV2 * t_inv);
        double nu = 1.0;
        for (int i = 1; i <= m; ++i) {
            pBuf[i] = t_inv * nu * pBuf[i-1];
            nu += 2.0;
        }
    }
}


void TlFmt::getGmT(const int m, const double& t, double* pBuf) const
{
    static const double COEF = std::sqrt(2.0 / TlFmt::PI);

    this->getFmT(m, t, pBuf);
    for (int i = 0; i <= m; ++i) {
        pBuf[i] *= COEF;
    }
}
