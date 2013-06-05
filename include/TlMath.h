#ifndef TLMATH_H
#define TLMATH_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif // M_PI

/** 数学的処理を行うクラス
 */
class TlMath {
public:
    /** x!を計算する
     *
     *  limit is 170!=7.257416...e+306 (MACHINE DEPEND)
     *  note such a overflow error is not supported!
     */
    double fact(int x);
    
    /** x!!を計算する
     *
     *  limit is 300!!=8.154414e+307 (MACHINE DEPEND)
     *  note such a overflow error is not supported!
     *  and note (-1)!! is 1, ( 0)!! is 2.
     */
    double dbfact(int x);
    
    /** (−1)のx乗を計算する
     */
    double minuspow(int x);
    
    // arc sinh を返す
    // 標準C++にはasinhは定義されていない
    static double arcsinh(const double x) {
        return std::log(x + std::sqrt(x * x +1));
    }
    
    /** 組み合わせ nCr の計算をする
     *
     *  limit is MACHINE DEPEND
     *  combination,  nCr ; n>=r is required
     */
    double combination(int n, int r);
    
    /** 順列 nPr の計算を行う
     *
     *  limit is MACHINE DEPEND
     *  permutation, nPr ; n>=r is required
     */
    double permutation(int n, int r);
    
    /**
     *  素因数分解
     */
    static std::vector<int> factor(int x);
    
    /// 指数が整数のpow()関数
    static double pow(double base, int exponent);

    static double pow(const double x, const double y);

    // function below utilized math.h
    static double PI() {
        return M_PI;
    };

    double sqrt(double x) {
        return ::sqrt(x);
    };

    double exp(double x) {
        return ::exp(x);
    };

    /// 桁落ちのチェック
    /// cancellation of significant digits
    static std::string checkCancellationErr(const std::string& file, int line,
                                            const double a, const double b);
};


inline double TlMath::pow(const double x, const double y)
{
    assert(x > 0.0);
    return (std::exp(y * std::log(x)));
}


inline double TlMath::pow(double base, int exponent)
{
    double answer = 1.0;
    if (exponent < 0) {
        base = 1.0 / base;
        exponent = std::abs(exponent);
    }
    
    for (int i = 0; i < exponent; ++i) {
        answer *= base;
    }

    return answer;
}


#endif // TLMATH_H

