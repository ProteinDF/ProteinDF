#include <iostream>
#include <string>
#include <cstring>
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

void reverse(char* s)                           // make reverse string s
{
    int     c, i, j;
    for (i = 0, j = std::strlen(s)-1; i < j; i++, j--) {
        c    = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

void itoa(int n, char* s)                       // change int n to char* s
{
    int     i, sign;
    if ((sign = n) < 0)                     // record sign
        n = -n;                         // n becomes positive
    i = 0;
    do {                                    // make string in a reverse
        // order
        s[i++] = n % 10 + '0';
    } while ((n /= 10) > 0);
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

/******************************************************************************/
/*  densPpq ; Ppq values in the order of the through index                */
/*        Note that the diagonal part of Ppq is already devided by 2. */
/*  npqA    ; the total number of pqA elements                            */
/*  nA  ; the number of A index                                       */
/*  inA     ; A index                                                     */
/*  npq ; the number of pq index                                      */
/*  inp     ; p index                                                     */
/*  inq     ; q index                                                     */
/*  pqA     ; pqA values                                                  */
/*      naux    ; the total number of auxialiy basis functions                */
/*  vect    ; the through index for integral pq                           */
/*  tA      ; the result vector                                           */
/******************************************************************************/
// calculate inner-product between Ppq and pqA
// return tA
// diagonal part of Ppq is already devided by 2
// this routine is called only one time per cycle
void spmtrprd(double* densPpq, int npqA, int nA, int* inA, int* npq,
              int* inp, int* inq, double* pqA, int naux, int* vect, double* tA)
{
    // common though-index of pqAlpha is thrown into inp[i]
    // if the values are something wrong, check here (for example, inp <-> inq)
    for (int i = 0; i < npqA; ++i) {
        vect[i] = (inp[i] + 1) * inp[i] / 2 + inq[i];
    }

    for (int i = 0; i < naux; ++i) {
        tA[i] = 0.0;
    }

    int   point = 0;          // index point

    for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < npq[i]; ++j) {
            tA[ inA[i] ] += densPpq[ vect[point] ] * pqA[point];
            point        += 1;
        }
    }

    for (int i = 0; i < naux; ++i) {
        tA[i] *= 2.0;
    }

}

/******************************************************************************/
/*                                                                            */
/*  Ppq     ; Ppq values in the order of the through index                */
/*        Note that the diagonal part of Ppq is already devided by 2. */
/*  npqA    ; the total number of pqA elements                            */
/*  nA  ; the number of A index                                       */
/*  inA     ; A index                                                     */
/*  npq ; the number of pq index                                      */
/*  inp     ; p index                                                     */
/*  inq     ; q index                                                     */
/*  pqA     ; pqA values                                                  */
/*      naux    ; the total number of auxialiy basis functions                */
/*  vect    ; the through index for integral pq                           */
/*  tA      ; the result vector                                           */
/*                                                                            */
/******************************************************************************/
// calculate inner-product between Ppq and pqA
// return tA
// diagonal part of Ppq is already devided by 2
// this routine is called only one time per cycle

// void spmtrprd(Fl_Matrix& Ppq, int npqA, int nA, int* inA, int* npq,
//        int* inp, int* inq, double* pqA, int naux, double* tA ) {

// // common though-index of pqAlpha is thrown into inp[i]
// // if the values are something wrong, check here (for example, inp <-> inq)

//   for (int i = 0; i < naux; ++i ){
//     tA[i] = 0.0;
//   }

//   int    point = 0;          // index point

//   for (int i = 0; i < nA; ++i ) {
//     for (int j = 0; j < npq[i]; ++j ) {
//       int indp          = inp[point];
//       int indq          = inq[point];
//       tA[ inA[i] ] += Ppq( indp, indq ) * pqA[point];
//       point        += 1;
//     }
//   }

//   for (int i = 0; i < naux; ++i ){
//     tA[i] *= 2.0;
//   }

// }

void spmtrprd(const TlSymmetricMatrix& Ppq, int npqA, int nA, int* inA, int* npq,
              int* inp, int* inq, double* pqA, int naux, double* tA)
{

// common though-index of pqAlpha is thrown into inp[i]
// if the values are something wrong, check here (for example, inp <-> inq)

    for (int i = 0; i < naux; ++i) {
        tA[i] = 0.0;
    }

    int   point = 0;          // index point

    for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < npq[i]; ++j) {
            int indp          = inp[point];
            int indq          = inq[point];
            tA[ inA[i] ] += Ppq(indp, indq) * pqA[point];
            point        += 1;
        }
    }

    for (int i = 0; i < naux; ++i) {
        tA[i] *= 2.0;
    }

}

void spmtrprd(const TlSymmetricMatrix& Ppq, int npqA, int nA, int* inA, int* npq,
              int* inp, int* inq, double* pqA, int naux, TlVector* ptA)
{

// common though-index of pqAlpha is thrown into inp[i]
// if the values are something wrong, check here (for example, inp <-> inq)

    for (int i = 0; i < naux; ++i) {
        (*ptA)[i] = 0.0;
    }

    int   point = 0;          // index point

    for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < npq[i]; ++j) {
            int indp          = inp[point];
            int indq          = inq[point];
            (*ptA)[ inA[i] ] += Ppq(indp, indq) * pqA[point];
            point        += 1;
        }
    }

    for (int i = 0; i < naux; ++i) {
        (*ptA)[i] *= 2.0;
    }

}


void spvctprd(double* MyuGamma, int nGDS, int nS, int* indexS, int* nGD,
              int* indexG, int* indexD, double* GDS, int na, double* tS)
{
// calculate inner-product between MyuGamma * MyuDelta and GDS
// for diagonal part * 1, and for off-diagonal part * 2

    int   i, j;
    double    Gamma, Delta, Integ;

    for (i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i])
            GDS[i] /= 2.0;

        //      Out << "GDS[" << i << "] = " << GDS[i] << "\n";

    }

    // We do not need to make common though-index.
    // The number of MyuGamma * MyuDelta elements is bigger than that of GDS.

    //    for ( i = 0; i < na; ++i )
    //        tS[i] = 0.0;

    int   point = 0;          // index point

    for (i = 0; i < nS; ++i) {
        for (j = 0; j < nGD[i]; ++j) {
            Gamma = MyuGamma[ indexG[point] ];
            Delta = MyuGamma[ indexD[point] ];
            Integ = GDS[ point ];
            tS[ indexS[i] ] += Gamma * Delta * Integ;

            point++;

            //            cout << "j = " << j << ", Gamma = " << Gamma
            //                 << ", Delta = " << Delta << ", Integ = " << Integ << "\n";
            //            cout << "indexS[i] = " << indexS[i]
            //                << ", G*D*I = " << Gamma*Delta*Integ << "\n";

        }

        //      Out << "tS[" << i << "] = " << tS[i] << "\n";

    }

    for (i = 0; i < na; ++i)
        tS[i] *= 2.0;

    for (i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i])
            GDS[i] *= 2.0;
    }
}

void spvctprd(const TlVector& MyuGamma, int nGDS, int nS, int* indexS, int* nGD,
              int* indexG, int* indexD, double* GDS, int na, double* tS)
{
// calculate inner-product between MyuGamma * MyuDelta and GDS
// for diagonal part * 1, and for off-diagonal part * 2

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] /= 2.0;
        }
    }

    // We do not need to make common though-index.
    // The number of MyuGamma * MyuDelta elements is bigger than that of GDS.

    //    for ( i = 0; i < na; ++i )
    //        tS[i] = 0.0;

    int   point = 0;          // index point

    for (int i = 0; i < nS; ++i) {
        for (int j = 0; j < nGD[i]; ++j) {
            const double Gamma = MyuGamma[indexG[point]];
            const double Delta = MyuGamma[indexD[point]];
            const double Integ = GDS[point];
            tS[indexS[i]] += Gamma * Delta * Integ;

            point++;
        }
    }

    for (int i = 0; i < na; ++i) {
        tS[i] *= 2.0;
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}


double spvctprd2(double* NyuGamma, int num, int* inp, int* inq, double* SGS)
{
    // calculate inner-product between NyuGamma * NyuSigma and SGS
    // for diagonal part * 1, and for off-diagonal part * 2

    double total = 0.0;

    for (int i = 0; i < num; ++i) {
        if (inp[i] == inq[i]) {
            SGS[i] /= 2.0;
        }
    }

    // We do not need to make common though-index.
    // The number of NyuGamma * NyuSigma elements is bigger than that of SGS.

    for (int i = 0; i < num; ++i) {
        const double Gamma = NyuGamma[ inp[i] ];
        const double Sigma = NyuGamma[ inq[i] ];
        total += Gamma * Sigma * SGS[i];
    }

    for (int i = 0; i < num; ++i) {
        if (inp[i] == inq[i]) {
            SGS[i] *= 2.0;
        }
    }

    return 2.0 * total;
}

double spvctprd2(const TlVector& NyuGamma, const TlSymmetricMatrix& SGS)
{
    // calculate inner-product between NyuGamma * NyuSigma and SGS
    // for diagonal part * 1, and for off-diagonal part * 2

    double total = 0.0;

    // We do not need to make common though-index.
    // The number of NyuGamma * NyuSigma elements is bigger than that of SGS.

    for (int i = 0; i < SGS.getNumOfRows(); ++i) {
        for (int j = 0; j <= i; j++) {
            const double Gamma = NyuGamma[i];
            const double Sigma = NyuGamma[j];

            if (i == j) {
                total += Gamma * Sigma * SGS(i, j) * 0.5;
            } else {
                total += Gamma * Sigma * SGS(i, j);
            }
        }
    }

    return 2.0 * total;
}


void spvctprd3(double* MyuGamma, double* NyuGamma, int nGDS, int nS, int* indexS, int* nGD, \
               int* indexG, int* indexD, double* GDS, int na, double* tS)
{
// calculate inner-product between MyuGamma * NyuGamma and GDS
// for diagonal part * 1, and for off-diagonal part * 2

    int i, j;
    double  Myu1, Myu2, Sigma, Integ;
    int Delta1, Delta2;

    for (i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i])
            GDS[i] /= 2.0;
    }

// We do not need to make common though-index.
// The number of MyuGamma * NyuGamma elements is bigger than that of GDS.

//  for ( i = 0; i < na; ++i )
//      tS[i] = 0.0;

    int point = 0;          // index point

    for (i = 0; i < nS; ++i) {
        Sigma = NyuGamma[ indexS[i] ];
        for (j = 0; j < nGD[i]; ++j) {
            Delta1 = indexG[ point ];
            Delta2 = indexD[ point ];
            Myu1   = MyuGamma[ Delta2 ];
            Myu2   = MyuGamma[ Delta1 ];
            Integ  = GDS[ point ];
            tS[ Delta1 ] += Myu1 * Sigma * Integ;
            tS[ Delta2 ] += Myu2 * Sigma * Integ;

            point++;
        }
    }

    for (i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i])
            GDS[i] *= 2.0;
    }

}

void spvctprd3(const TlVector& MyuGamma, const TlVector& NyuGamma, int nGDS, int nS, int* indexS, int* nGD, \
               int* indexG, int* indexD, double* GDS, int na, double* tS)
{
    // calculate inner-product between MyuGamma * NyuGamma and GDS
    // for diagonal part * 1, and for off-diagonal part * 2

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] /= 2.0;
        }
    }

// We do not need to make common though-index.
// The number of MyuGamma * NyuGamma elements is bigger than that of GDS.

    int point = 0;            // index point

    for (int i = 0; i < nS; ++i) {
        const double Sigma = NyuGamma[ indexS[i] ];
        for (int j = 0; j < nGD[i]; ++j) {
            const int Delta1 = indexG[ point ];
            const int Delta2 = indexD[ point ];
            const double Myu1   = MyuGamma[ Delta2 ];
            const double Myu2   = MyuGamma[ Delta1 ];
            const double Integ  = GDS[ point ];
            tS[ Delta1 ] += Myu1 * Sigma * Integ;
            tS[ Delta2 ] += Myu2 * Sigma * Integ;

            point++;
        }
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}


void spvctprd4(double* NyuGamma, int num, int* inp, int* inq, double* SGD, int naux, double* tS)
{
    // calculate inner-product between NyuGamma and SGD
    // for diagonal part * 1, and for off-diagonal part * 2

    for (int i = 0; i < num; ++i) {
        if (inp[i] == inq[i]) {
            SGD[i] /= 2.0;
        }
    }

    // We do not need to make common though-index.
    // The number of NyuGamma elements is bigger than that of SGS.

    //    for ( i = 0; i < naux; ++i )
    //        tS[i] = 0.0;

    for (int i = 0; i < num; ++i) {
        const int Sigma1 = inp[i];
        const int Sigma2 = inq[i];
        const double Nyu1 = NyuGamma[Sigma2];
        const double Nyu2 = NyuGamma[Sigma1];
        tS[Sigma1] += Nyu1 * SGD[i];
        tS[Sigma2] += Nyu2 * SGD[i];
    }

    for (int i = 0; i < num; ++i) {
        if (inp[i] == inq[i]) {
            SGD[i] *= 2.0;
        }
    }
}

void spvctprd4(const TlVector& NyuGamma, const TlSymmetricMatrix& SGD, int naux, double* tS)
{
    // calculate inner-product between NyuGamma and SGD
    // for diagonal part * 1, and for off-diagonal part * 2

    // We do not need to make common though-index.
    // The number of NyuGamma elements is bigger than that of SGS.

    //    for ( i = 0; i < naux; ++i )
    //        tS[i] = 0.0;

//   for (int i = 0; i < num; ++i){
//     const int Sigma1 = inp[i];
//     const int Sigma2 = inq[i];
//     const double Nyu1 = NyuGamma[Sigma2];
//     const double Nyu2 = NyuGamma[Sigma1];

//     if (Sigma1 == Sigma2){
//       tS[Sigma1] += Nyu1 * SGD[i] * 0.5;
//       tS[Sigma2] += Nyu2 * SGD[i] * 0.5;
//     } else {
//       tS[Sigma1] += Nyu1 * SGD[i];
//       tS[Sigma2] += Nyu2 * SGD[i];
//     }
//   }

    for (int i = 0; i < SGD.getNumOfRows(); i++) {
        for (int j = 0; i <= j; j++) {
            const int Sigma1 = i;
            const int Sigma2 = j;
            const double Nyu1 = NyuGamma[Sigma2];
            const double Nyu2 = NyuGamma[Sigma1];

            if (Sigma1 == Sigma2) {
                tS[Sigma1] += Nyu1 * SGD(i, j) * 0.5;
                tS[Sigma2] += Nyu2 * SGD(i, j) * 0.5;
            } else {
                tS[Sigma1] += Nyu1 * SGD(i, j);
                tS[Sigma2] += Nyu2 * SGD(i, j);
            }
        }
    }

}
