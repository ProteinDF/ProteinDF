#ifndef TLFMT_H
#define TLFMT_H

#include <cmath>
#include <limits>
#include <vector>
#include <map>

/** 補助積分クラス
 */
class TlFmt {
public:
    static TlFmt& getInstance();

public:
    enum TLFMT_FLAG {
        TLFMT_ENABLE_CACHE = 1 // 計算キャッシュを有効にする
    };

private:
    explicit TlFmt(int nFlag = 0);
    ~TlFmt();

public:
    /// Fm(T) を返す
    double getFmT0(const int m, const double& t) const;

    void getFmT(const int m, const double& t, double* pBuf) const;

    void getGmT(const int m, const double& t, double* pBuf) const;
    
private:
    /// Taylor展開用のテーブルを作成する
    void makeFmtTable();

    /// getFmt() の実計算部
    std::vector<double> getFmtSubstance(const int m, const double& t) const;

private:
    static TlFmt* m_pInstance;

    int m_nFlags;

    /// Taylor展開用テーブル
    ///
    /// サイズは[FMT_MAX_M +1][FMT_N_STEP +1]
    //std::vector<std::vector<double> > m_FmtTable;
    double** m_FmtTable;

    // すでに計算したFm(T)の値を保持しておくキャッシュ
    // キーは t の値、値はf[m]
    // const オブジェクトでもこのキャッシュは関係ないので、mutable宣言しています。
    mutable std::map<double, std::vector<double> > m_aFmtCalcCache;

    // 前回計算したFm(T)のTの値
    mutable double m_previousT;
    // 前回計算したFm(T)の値を格納した配列
    mutable std::vector<double> m_previousF;

    mutable unsigned int m_nNumOfFunctionCalled; // getFmt()が呼ばれた回数
    mutable unsigned int m_nNumOfCacheHits;      // 計算キャッシュが成功した回数

    // 定数
    static const double PI;
    static const double CK;
    static const double PI_DIV2;
    static const double INV2;
    static const double INV6;
    static const int    FMT_M;
    static const int    FMT_INV_D;
    static const int    FMT_MAX_M;
    static const double FMT_T;
    static const int    FMT_N_STEP;
    static const double FMT_STEPSIZE;
    static const double INV_FMT_STEPSIZE;
    static const double EPSILON;    // double の計算機ε

};

#endif // TLFMT_H
