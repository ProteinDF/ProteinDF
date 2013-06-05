#include <cassert>
#include <cmath>
#include "TlMath.h"
#include "TlUtils.h"

/* limit is 170!=7.257416...e+306 (MACHINE DEPEND). */
/* note such a overflow error is not supported!     */
double TlMath::fact(int x)
{
    assert((0 <= x) && (x <= 170));

    double ans =1.0;
    while (x>1) {
        ans *= (double)x;
        x--;
    }

    return ans;
}

/* limit is 300!!=8.154414e+307 (MACHINE DEPEND). */
/* note such a overflow error is not supported!   */
/* and note (-1)!! is 1, ( 0)!! is 2.             */
double TlMath::dbfact(int x)
{
    assert((-1 <= x) && (x <= 300));

    double  ans;

    if (x % 2 == 0) {
        // when x is even
        ans = 2.0;
        while (x > 2) {
            ans *= (double)x;
            x -= 2;
        }
    } else {
        // when x is odd
        ans = 1.0;
        while (x > 1) {
            ans *= (double)x;
            x -= 2;
        }
    }

    return ans;
}

/* this function returns value, pow( -1.0, (double)x ). */
double TlMath::minuspow(int x)
{
    if (x%2==0) return 1.0;
    else return -1.0;
}

/* limit is MACHINE DEPEND.             */
/* combination,  nCr ; n>=r is required */
double TlMath::combination(int n, int r)
{
    assert(n >= r);
    assert(n >= 0);
    assert(r >= 0);

    int sn = n;
    int sr = r;
    double ans = 1.0;
    while (sn>(n-r) || sr>1) {
        if (sn > (n-r)) {
            ans *= (double)sn;
            sn--;
        }

        if (sr > 1) {
            ans /= (double)sr;
            sr--;
        }
    }

    return ans;
}

/* limit is MACHINE DEPEND.            */
/* permutation, nPr ; n>=r is required */
double TlMath::permutation(int n, int r)
{
    assert(n >= r);
    assert(n >= 0);
    assert(r >= 0);

    double ans = 1.0;
    for (int i = n; i>(n-r); i--) {
        ans *= (double)i;
    }

    return ans;
}

std::vector<int> TlMath::factor(const int nInput)
{
    std::vector<int> answer;

    if (nInput < 2) {
        answer.push_back(nInput);
        return answer;
    }

    int nMaxFactor = static_cast<int>(std::sqrt(static_cast<double>(nInput)));
    if ((nMaxFactor * nMaxFactor) < nInput) {
        nMaxFactor++;
    }

    int num = nInput;
    for (int nFactor = 2; nFactor <= nMaxFactor; ++nFactor) {
        int nRemain = 0;
        while ((nRemain == 0) && (num >= 2)) {
            std::ldiv_t sDiv = std::ldiv(num, nFactor);
            nRemain = sDiv.rem;
            if (nRemain == 0) {
                num = sDiv.quot;
                answer.push_back(nFactor);
            }
        }
    }

    if (num > nMaxFactor) {
        answer.push_back(num);
    }

    return answer;
}


std::string TlMath::checkCancellationErr(const std::string& file, int line,
                                         const double a, const double b)
{
    const double c = a - b;
    return TlUtils::format("%s:%d: %16.10e - %16.10e = %16.10e",
                           file.c_str(), line,
                           a, b, c);
}
