#ifndef DFXCPOTFITTING_H
#define DFXCPOTFITTING_H

#include "TlSymmetricMatrix.h"
#include "Fl_Int_Pqg.h"
#include "DfObject.h"

/*****************************************************************************/
/*                                                                           */
/*   Note:                                                                   */
/*  We have not coded the member functions about Numerical Local Spin    */
/*  Density Methods yet.                                                 */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/*   SPECIFICATIONS                                                          */
/*                                                                           */
/*   Name:                                                                   */
/*      Class DfXcpotfitting                                                 */
/*                                                                           */
/*   Dates:                                                                  */
/*      Version 1.0; F. Sato, 1994, 4, 22                                    */
/*      Version 1.1; F. Sato, 1994, 5, 13                                    */
/*                                                                           */
/*   Work:                                                                   */
/*      Calclation of Exchange-Correlation Potential Coefficients            */
/*      by Analytical X-Alpha or Numerical Local Spin Density Methods        */
/*                                                                           */
/*      Then, this class has 2 main member functions of procedure-type       */
/*         1; Analytical X-Alpha Method (See Ref.1)                          */
/*         2; Local Spin Density Method by Mesh Integration (See Ref.2-4)    */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/*   THEORY                                                                  */
/*                                                                           */
/*  Analytical X-Alpha Method                                            */
/*                                                                           */
/*            1/3     na                                                     */
/*  1. Rou( r ) = Sum  MyuGamma * gGamma( r )                            */
/*                    Gamma                                                  */
/*                                                                           */
/*            2/3     na               y                                     */
/*  2. Rou( r ) = Sum  NyuGamma * g Gamma( r )                           */
/*                    Gamma                                                  */
/*                                                                           */
/*                    n,na                                                   */
/*  3. E1       = Sum  MyuGamma * Ppq [ pqGamma ]                        */
/*                    pqGamma                                                */
/*                                                                           */
/*                    na                                                     */
/*  4. E2       = Sum  MyuGamma * MyuDelta * NyuSigma [ GammaDeltaSigma ]*/
/*                GammaDeltaSigma                                            */
/*                                                                           */
/*                    na                                                     */
/*  5. E3       = Sum  NyuGamma * NyuSigma [ GammaSigma ]                */
/*                  GammaSigma                                               */
/*                                                                           */
/*                                  2  1/3                                   */
/*      6. MyuGamma = ( E1 * E2 / E3  )  MyuGamma                            */
/*                                                                           */
/*                        2            1/3                                   */
/*      7. NyuGamma = ( E1  / E2 * E3 )  NyuGamma                            */
/*                                                                           */
/*  where na and n are the total numbers of auxiliary basis and orbital  */
/*  functions, respectively.                                             */
/*  gGamma and gyGamma is auxiliary basis functions of Myu and Nyu,      */
/*  respectively.                                                        */
/*                                                                           */
/*  Notes:                                                               */
/*     1. Because Myu = 4 Alpha /3 Rou^1/3, MyuGamma must be scaled      */
/*        before and after this calculation.  See ==> each comment       */
/*     2. If this calculation is costed, we must consider the            */
/*        possibility of extrapolation of E1, E2 and E3.                 */
/*    ___________________________________________________________________    */
/*                                                                           */
/*  Numerical Local Spin Density Method                                  */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/*   FOR PROGRAMMERS                                                         */
/*                                                                           */
/*   Calling Sequences:                                                      */
/*                                                                           */
/*   I/O:                                                                    */
/*      I/O files:                                                           */
/*         There are 2 types of I/O file groups corresponded to the calc-    */
/*         lation methods.                                                   */
/*            1; Analytical X-alpha Method                                   */
/*               Read files-                                                 */
/*                  fl_Vct_Myu(int-1)                                        */
/*                  fl_Vct_Nyu(int-1)                                        */
/*                  fl_Mtr_Ppq                                               */
/*                  fl_Mtr_Sgd                                               */
/*                  fl_Int_Pqg            by Disk Methods                    */
/*                  fl_Int_Gds            by Disk Methods                    */
/*               Write files-                                                */
/*                  fl_Vct_Myu(int)                                          */
/*                  fl_Vct_Nyu(int)                                          */
/*               Note: fl_Vct_Myu(1) = fl_Vct_Nyu(1) = fl_Vct_Rou(1)         */
/*             or to be given at DfInitialguess                      */
/*                                                                           */
/*            2; Numerical Local Spin Density Method                         */
/*               Read files-                                                 */
/*                  fl_Vct_Rou(int)                                          */
/*                  fl_Gto_Density                                           */
/*                  fl_Gto_Xcpot                                             */
/*                  fl_Mtr_Mesh                                              */
/*                  fl_Mtr_Sgdinv                                            */
/*               Write files-                                                */
/*                  fl_Vct_Myu(int)                                          */
/*                                                                           */
/*      I/O methods:                                                         */
/*                                                                           */
/*   Interaction:                                                            */
/*      In case of Direct Methods, we have interaction with 3-index          */
/*      integrations.                                                        */
/*                                                                           */
/*      1; Analytical X-alpha Methods                                        */
/*         [ p     q     gamma ]                                             */
/*         [ gamma delta sigma ]                                             */
/*                                                                           */
/*      2; Numerical Local Spin Density Method                               */
/*         No-interaction                                                    */
/*                                                                           */
/*   Flags:                                                                  */
/*  From Fl_Globalinput, it gets Run_Level, Method_Type and Alpha_Value  */
/*  flags.  The details are not determined yet.                          */
/*                                                                           */
/*   Warnings & Errors:                                                      */
/*  1. Constructor and Destructor check                                  */
/*  2. Memory check                                                      */
/*  3. File check on reading and writing                                 */
/*  4. Element number check of Files                                     */
/*                                                                           */
/*   Debugs:                                                                 */
/*  1. Vectors and Matrixes can be printed out on display, if the        */
/*     comment-outed cout routines in member functions are taken out     */
/*     and recompiled.                                                   */
/*  2. The minus Run_Level is debug output.                              */
/*  3. This class has temporal error output member function.             */
/*                                                                           */
/*   Notes:                                                                  */
/*      We have not coded the member functions about Numerical Local Spin    */
/*      Density Methods yet.                                                 */
/*                                                                           */
/*   Reference:                                                              */
/*      1. B. I. Dunlap & N. Ro(um)sch: "On The Gaussian-Type Orbitals       */
/*         Approach to Local Density Functional Theory", Journal de chimie   */
/*         physique, 86 (1989) 671-688.                                      */
/*                                                                           */
/*   Further Reference for Extension:                                        */
/*      1. B. I. Dunlap, J. Andzelm & J. W. Mintmire: "Local-Density-        */
/*         Functional Total Energy Gradients in The Linear Combination of    */
/*         Gaussian-Type Orbitals Method", Physical Review, A42 (1990) 6354  */
/*         -6359                                                             */
/*                                                                           */
/*****************************************************************************/
class DfXcpotfitting : public DfObject {
    // for Mesh
public:
    DfXcpotfitting(TlSerializeData* pPdfParam,int nItr);
    ~DfXcpotfitting();

    int dfXcpMain();

private:
//   int dfXcpEstmemval();
    int getMemory();

    // for X-Alpha
    int initAtbl();       // Initialize Alpha Table
    int readPpq();        // read Ppq, MyuGamma and

    // NyuGamma
    int   readSgd();          // read Sgd
    int   bfscalMyu();        // scaling MyuGamma before

    // calculation
    int   newton(int);        // Newton-Raphson
    int   calcE1();       // calculation of E1
    int   calcE2();       // calculation of E2
    int   calcE3();       // calculation of E3
    int   calcMN();       // calculation of MyuGamma

    // and NyuGamma
    int   afscalMyu();        // scaling MyuGamma after
    // calculation
    void Error(const char*);      // temp error output
    void Print(int, double*); // temp vector output
    void Print(int, int*, int*, double*);

private:
    // Common functions
    double* tmpvectorA;           // temporal vector
    double* tmpvectorB;
    double* tmpE2A;
    double* tmpE2B;
    double* tAlpha;
    double* tBeta;
    // for X-Alpha
    // following vectors are
    // max 80 kB/1 vector

//   double* MyuGammaA;         // MyuGamma
//   double* MyuGammaB;
//   double* NyuGammaA;         // NyuGamma
//   double* NyuGammaB;
    TlVector MyuGammaA;
    TlVector MyuGammaB;
    TlVector NyuGammaA;
    TlVector NyuGammaB;

    // These are the work area for reading the pqGamma, GDSigma and Sgd.
    // 3-index integrals are [ pqGamma ] and [ GammaDeltaSigma ] and we call
    // p and Gamma are first index, q and Delta are second ones, and so on.
    // index0 and n12 are max 80 kB, respectively.
    int* index0;              // index for the first elements
    int* n12;             // the number of the second-
    // third elements
    int* index1;              // index for the second elements
    // or row for Sgd
    int* index2;              // index for the third elements
    // or col for Sgd
    double* Intval;               // pqGamma or GDSigma or Sgd

    // Ppq = row + col + Ppqval; Ppq is max 3 MB.
//   int* rowA;
//   int* colA;
//   double* PpqA;
    TlSymmetricMatrix PpqA;
//   int* rowB;
//   int* colB;
//   double* PpqB;
    TlSymmetricMatrix PpqB;

    // Sgd = row + col + Sgdval; Sgd is max 15 MB.
//   int* rowS;
//   int* colS;
//   double* Sgd;
    TlSymmetricMatrix Sgd;

    // Newton-Raphson
    double** Jacobi;
    double** InvJacobi;
    double** EUI;         // for check
    double* colinv;
    double* tmpvector2;
    int* indx;
    double*       X;
    double*       deltaX;
    double*       gSigmaA;
    double*       fDeltaA;
    double*       gSigmaB;
    double*       fDeltaB;
    double        NREPS;
    int       NRiter;
    double        DFactor;

    // Scaling
    double        E1A;                // E1
    double        E1B;
    double        E2A;                // E2
    double        E2B;
    double        E3A;                // E3
    double        E3B;
    double        SFactorMA;          // Scaling Factor
    double        SFactorMB;
    double        SFactorNA;
    double        SFactorNB;

    double*       densPpqA;           // densed Ppq
    double*       densPpqB;           // densed Ppq
    int*      workvector1;            // work vector
    int*      workvector2;

    int*      atomnumtbl;         // Atom Number Table
    double*       Alpha;              // scale factor
    // usually Alpha = 0.7

    int outlevel;                 // Run Level
    double alphaval;          // uniform alpha value
    int            naux;              // total number of auxiliary BFs
    int            norb;              // total number of orbital BFs
//   int            nPpqA;              // total number of Ppq elements
//   int            nPpqB;
//   int            nSgd;               // total number of Sgd elements
    //int vectorelement;          // max element number of any vectors
    int TotalPpq;         // total number of non-contract Ppq

    std::string scftype;
    int number_iteration;         // File iteration suffix (forward)
    int niteration1;          // File iteration suffix (current)
    int number_rotation;

};

#endif

