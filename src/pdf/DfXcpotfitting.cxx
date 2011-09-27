#include <cmath>
#include <cassert>

#include "DfXcpotfitting.h"
#include "Common.h"
#include "TlUtils.h"
#include "TlTime.h"
#include "TlLogX.h"
#include "Fl_Int_Gds.h"

#define F13             0.3333333333333333
// #define F23             0.6666666666666667
// #define F43             1.3333333333333333
#define R3_2            1.2599210498948732

extern void   spmtrprd(double*, int,  int,    int*, int*, int*,   int*, double*, int,  int*, double*);
extern void   spvctprd(const TlVector&, int, int, int*, int*, int*, int*, double*, int, double*);
extern double spvctprd2(const TlVector&, const TlSymmetricMatrix&);
extern void   spvctprd3(const TlVector&, const TlVector&, int, int, int*, int*, int*, int*, double*, int, double*);
extern void   spvctprd4(const TlVector&, const TlSymmetricMatrix&, int, double*);


// make matrix elements of A2 & A3 part
// if A3 = Bij, A2 = -2Bji
// diagonal parts are doubly counted
void putA23(double* MyuGamma, int nGDS, int nS, int* indexS, int* nGD, int* indexG, int* indexD, double* GDS, int na, double** Jacobi)
{
    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 0.5; // same as /= 2.0
        }
    }

    int   point = 0;          // index point

    for (int i = 0; i < nS; ++i) {
        const int Sigma = indexS[i];
        for (int j = 0; j < nGD[i]; ++j) {
            const int IGamma = indexG[ point ];
            const int IDelta = indexD[ point ];
            const double Gamma  = MyuGamma[ IGamma ];
            const double Delta  = MyuGamma[ IDelta ];
            const double IntegG = Gamma * GDS[ point ];
            const double IntegD = Delta * GDS[ point ];

            // A3
            Jacobi[Sigma][IDelta + na] += -2.0 * IntegG;
            Jacobi[Sigma][IGamma + na] += -2.0 * IntegD;

            // A2
            Jacobi[IDelta + na][Sigma] +=        IntegG;
            Jacobi[IGamma + na][Sigma] +=        IntegD;

            point++;
        }
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}

void putA23(const TlVector& MyuGamma, int nGDS, int nS, int* indexS, int* nGD, int* indexG, int* indexD, double* GDS, int na, double** Jacobi)
{
    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 0.5; // same as /= 2.0
        }
    }

    int   point = 0;          // index point

    for (int i = 0; i < nS; ++i) {
        const int Sigma = indexS[i];
        for (int j = 0; j < nGD[i]; ++j) {
            const int IGamma = indexG[ point ];
            const int IDelta = indexD[ point ];
            const double Gamma  = MyuGamma[ IGamma ];
            const double Delta  = MyuGamma[ IDelta ];
            const double IntegG = Gamma * GDS[ point ];
            const double IntegD = Delta * GDS[ point ];

            // A3
            Jacobi[Sigma][IDelta + na] += -2.0 * IntegG;
            Jacobi[Sigma][IGamma + na] += -2.0 * IntegD;

            // A2
            Jacobi[IDelta + na][Sigma] +=        IntegG;
            Jacobi[IGamma + na][Sigma] +=        IntegD;

            point++;
        }
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}


// make matrix elements of A4 part
// diagonal parts are doubly counted
void putA4(double* NyuGamma, int nGDS, int nS, int* indexS, int* nGD, int* indexG, int* indexD, double* GDS, int na, double** Jacobi)
{
    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 0.5; // same as /= 2.0
        }
    }

    int   point = 0;          // index point

    for (int i = 0; i < nS; ++i) {
        const double Sigma = NyuGamma[ indexS[i] ];
        for (int j = 0; j < nGD[i]; ++j) {
            const int Gamma = indexG[point] + na;
            const int Delta = indexD[point] + na;
            // + na is to put into A4
            const double Integ = GDS[ point ];
            Jacobi[Gamma][Delta] += Sigma * Integ;
            Jacobi[Delta][Gamma] += Sigma * Integ;

            point++;
        }
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}

void putA4(const TlVector& NyuGamma, int nGDS, int nS, int* indexS, int* nGD, int* indexG, int* indexD, double* GDS, int na, double** Jacobi)
{
    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 0.5; // same as /= 2.0
        }
    }

    int   point = 0;          // index point

    for (int i = 0; i < nS; ++i) {
        const double Sigma = NyuGamma[ indexS[i] ];
        for (int j = 0; j < nGD[i]; ++j) {
            const int Gamma = indexG[point] + na;
            const int Delta = indexD[point] + na;
            // + na is to put into A4
            const double Integ = GDS[ point ];
            Jacobi[Gamma][Delta] += Sigma * Integ;
            Jacobi[Delta][Gamma] += Sigma * Integ;

            point++;
        }
    }

    for (int i = 0; i < nGDS; ++i) {
        if (indexG[i] == indexD[i]) {
            GDS[i] *= 2.0;
        }
    }
}

void LUdcmp(double** a, int n, int* indx, double* d, double* vv)
{
    const double TINY = 1.0e-30;

    *d = 1.0;
    for (int i = 0; i < n; ++i) {
        double big = 0.0;
        for (int j = 0; j < n; ++j) {
            const double temp = fabs(a[i][j]);
            if (temp > big) {
                big = temp;
            }
        }
        if (big == 0.0) {
            std::cout << "*** Singular matrix in routine LUdcmp ***" << std::endl;
            std::cout << std::endl;
        }
        vv[i] = 1.0 / big;
    }

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) {
            double sum = a[i][j];
            for (int k = 0; k < i; ++k) {
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;
        }
        double big = 0.0;
        int imax = 0;
        for (int i = j; i < n; ++i) {
            double sum = a[i][j];
            for (int k = 0; k < j; ++k) {
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;
            const double dum = vv[i] * fabs(sum);
            if (dum >= big) {
                big  = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (int k = 0; k < n; ++k) {
                const double dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k]    = dum;
            }
            *d       = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) {
            a[j][j] = TINY;
        }

        if (j != n-1) {
            const double dum = 1.0/a[j][j];
            for (int i = j+1; i < n; ++i) {
                a[i][j] *= dum;
            }
        }
    }
}

void LUbksb(double** a, int n, int* indx, double* b)
{
    int ii = -1;
    for (int i = 0; i < n; ++i) {
        const int ip = indx[i];
        double sum = b[ip];
        b[ip] = b[i];
        if (ii != -1) {
            for (int j = ii; j <= i-1; ++j) {
                sum -= a[i][j] * b[j];
            }
        } else if (sum) {
            ii = i;
        }
        b[i] = sum;
    }

    for (int i = n-1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i+1; j < n; ++j) {
            sum -= a[i][j]*b[j];
        }
        b[i] = sum / a[i][i];
    }
}

DfXcpotfitting::DfXcpotfitting(TlSerializeData* pPdfParam, int nItr)
    : DfObject(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
//   TlLogX& Log = TlLogX::getInstance();

//     this->number_iteration = nItr;
//     this->scftype = flGbi["SCF"]["method"];
//     this->number_rotation = atoi(flGbi["SCF"]["file-rot-number"].c_str());
//     this->naux = atoi(flGbi["SCF"]["control-nauxxc"].c_str());
//     this->norb = atoi(flGbi["SCF"]["control-norb"].c_str());
    this->TotalPpq = (this->norb + 1) * this->norb / 2;
//   const int memory = this->dfXcpEstmemval();

//     this->outlevel = atoi(flGbi["SCF"]["print-level"].c_str());
    this->alphaval = pdfParam["model"]["xc-potential/xalpha/alpha-value"].getDouble();
    if (fabs(this->alphaval) <= 1.0E-16) {
        this->alphaval = 0.7;
    }

    this->NREPS    = pdfParam["model"]["xc-potential/xalpha/converge-threshold"].getDouble();
    this->NRiter   = pdfParam["model"]["xc-potential/xalpha/iteration-number"].getInt();
    this->DFactor  = pdfParam["model"]["xc-potential/xalpha/damping-factor"].getDouble();

//   if (memory > MaxMemory){
//     Error( "Cannot get memory." );
//     Log << "Used Memory is " << 8 * memory << " [Byte]" << "\n\n";
//   } else {
//     Log << "Used Memory is " << 8 * memory << " [Byte]" << "\n\n";
    this->getMemory();
//   }
}

DfXcpotfitting::~DfXcpotfitting()
{
    delete[] Alpha;
    delete[] tBeta;
    delete[] tAlpha;
    delete[] tmpE2B;
    delete[] tmpE2A;
    delete[] tmpvectorB;
    delete[] tmpvectorA;
    delete[] densPpqB;
    delete[] densPpqA;
    delete[] workvector2;
    delete[] workvector1;
    delete[] n12;
    delete[] index0;
    delete[] Intval;
    delete[] index2;
    delete[] index1;
    delete[] fDeltaB;
    delete[] gSigmaB;
    delete[] fDeltaA;
    delete[] gSigmaA;
    delete[] deltaX;
    delete[] X;
    delete[] indx;
    delete[] tmpvector2;
    delete[] colinv;
    for (int i = 0; i < 2 * this->naux; ++i) {
        delete[] *(InvJacobi + i);
    }
    delete[] InvJacobi;
    for (int i = 0; i < 2 * this->naux; ++i) {
        delete[] *(Jacobi + i);
    }
    delete[] Jacobi;
}

// int DfXcpotfitting::dfXcpEstmemval() {
//   return  MAXNPQA * 3 + MAXNA * 2 + this->naux * 25 + MaxEleNum * 7 + MaxEleAux * 3 + this->naux * this->naux * 8 + this->TotalPpq * 2;
// }

int DfXcpotfitting::getMemory()
{
    Jacobi = new double*[2 * this->naux + 1];
    for (int i = 0; i < (2 * this->naux); ++i) {
        *(Jacobi + i) = new double[2 * this->naux];
    }
    InvJacobi = new double*[2 * this->naux + 1];
    for (int j = 0; j < (2 * this->naux); ++j) {
        *(InvJacobi + j) = new double[ 2 * this->naux + 1 ];
    }

    colinv      = new double[ 2 * this->naux ];
    tmpvector2  = new double[ 2 * this->naux ];
    indx        = new int  [ 2 * this->naux ];
    X           = new double[ 2 * this->naux ];
    deltaX      = new double[ 2 * this->naux ];
    gSigmaA     = new double[ this->naux ];
    fDeltaA     = new double[ this->naux ];
    gSigmaB     = new double[ this->naux ];
    fDeltaB     = new double[ this->naux ];

    index1      = new int  [ MAXNPQA ];
    index2      = new int  [ MAXNPQA ];
    Intval      = new double[ MAXNPQA ];
    index0      = new int  [ MAXNA ];
    n12         = new int  [ MAXNA ];

    workvector1 = new int  [ MaxEleNum ];
    workvector2 = new int  [ MaxEleNum ];
    densPpqA    = new double[this->TotalPpq];
    densPpqB    = new double[this->TotalPpq];

//   MyuGammaA   = new double[ this->naux ];
//   MyuGammaB   = new double[ this->naux ];
//   NyuGammaA   = new double[ this->naux ];
//   NyuGammaB   = new double[ this->naux ];
    tmpvectorA  = new double[ this->naux ];
    tmpvectorB  = new double[ this->naux ];
    tmpE2A      = new double[ this->naux ];
    tmpE2B      = new double[ this->naux ];
    tAlpha      = new double[ this->naux ];
    tBeta       = new double[ this->naux ];
    Alpha       = new double[ this->naux ];

//   rowA        = new int  [ MaxEleNum ];
//   colA        = new int  [ MaxEleNum ];
//   PpqA        = new double[ MaxEleNum ];
//   rowB        = new int  [ MaxEleNum ];
//   colB        = new int  [ MaxEleNum ];
//   PpqB        = new double[ MaxEleNum ];
//   rowS        = new int  [ MaxEleAux ];
//   colS        = new int  [ MaxEleAux ];
//   Sgd         = new double[ MaxEleAux ];

    return 0;
}

// Note: Main calculation parts of calcE1, calcE2 and calcE3 are essentially
//   same (Matrix * Vector).
//   In the c++ style these member functions should be governed, but
//   in order to emphasize the procedures, we would make different
//   functions.
int DfXcpotfitting::dfXcpMain()
{
    TlLogX& Log = TlLogX::getInstance();
//   Log.precision(8);
//   Log.setf(ios::scientific, ios::floatfield);

    Log << "=== START ==========\n";
    Log << "START : " << TlTime::getNow() << "\n";

    this->initAtbl();
    Log << "initAtbl : " << TlTime::getNow() << "\n";
    Log.flush();

    this->readPpq();
    Log << "readPpq  : " << TlTime::getNow() << "\n";
    Log.flush();

    this->readSgd();
    Log << "readSgd  : " << TlTime::getNow() << "\n";
    Log.flush();

    this->bfscalMyu();
    Log << "bfscalMyu : " << TlTime::getNow() << "\n";
    Log.flush();

    int   flag = -1;      // if flag = 0, converged
    int   mem  = NRiter+1;    // memory of NRiter
    while (NRiter * flag) {
        Log << TlUtils::format("N-R iter = %3d", mem-NRiter);

        flag = this->newton(mem-NRiter);
        NRiter--;
    }
    Log << "newton  : " << TlTime::getNow() << "\n";
    Log.flush();

    if (flag) {
        this->calcE1();
        Log << "calcE1 : " << TlTime::getNow() << "\n";
        Log.flush();

        this->calcE2();
        Log << "calcE2 : " << TlTime::getNow() << "\n";
        Log.flush();

        this->calcE3();
        Log << "calcE3 : " << TlTime::getNow() << "\n";
        Log.flush();

        this->calcMN();
        Log << "calcMN : " << TlTime::getNow() << "\n";
        Log.flush();
    }

    this->afscalMyu();
    Log << "afscalMyu : " << TlTime::getNow() << "\n";
    Log.flush();

    Log << "END : " << TlTime::getNow() << "\n";
    Log << "=== END ============\n";
    Log.flush();

    return 0;
}

// Make Alpha Table
// 1. get option about Alpha value from Fl_Inputdata
// 2. analize the option
// 3. All Alpha values of each atom are determined, referring TlAlphatbl
//    which contain some Alpha table sets.
//        a. construct TlAlphatbl giving option
//           default or user defined or some optional set
//        b. get Alpha value of each atom by double getalphaval(int) in
//           TlAlphatbl
// 4. make atomnumtbl from int getAuxtoAnum(int) in DfInitialguess
//    each auxiliary function beints to what kind of atoms
// 5. make Alpha Table of each auxiliary function
int DfXcpotfitting::initAtbl()
{
    double coef = - alphaval * 3.0 / 2.0 * pow(3.0 / M_PI, F13);

    // if this class must scaled, the following sentences are needed.
    if (this->scftype == "sp") {
        coef *= R3_2;
    }

    for (int i = 0; i < naux; ++i) {
        this->Alpha[i] = coef;
    }

    return 0;
}

// Read Ppq matrix and put into row, col and Ppq.
// Because Ppq is a sparse matrix, all elements can be put in the main memory
// and be read at one time.
// Read MyuGamma & NyuGamma
int DfXcpotfitting::readPpq()
{
    int niteration  = this->number_iteration % this->number_rotation - 1;
    if (niteration == 0) {
        niteration = this->number_rotation;
    }
    if (niteration == -1) {
        niteration = this->number_rotation - 1;
    }

    this->niteration1 = this->number_iteration % this->number_rotation;
    if (this->niteration1 == 0) {
        this->niteration1 = this->number_rotation;
    }
    // rotational numbers are
    // 1, 2, ..., nrot

    if (this->scftype == "nsp") {
//     Fl_Mtr_PpqD  ppq( RSFD,niteration, "r" );
//     ppq.open( "fl_Work", ppq.getfilename(), "read");
//     ppq.getelmnum( &nPpqA );
//     ppq.readc( nPpqA, rowA, colA, PpqA );
//     ppq.close();
        this->PpqA.load("fl_Work/fl_Mtr_PpqDr" + TlUtils::xtos(niteration));

//     int totalmA;
//     Fl_Vct_Myu myu(niteration);
//     myu.open("fl_Work", myu.getfilename(), "read");
//     myu.getelemnum(&totalmA);
//     myu.read(totalmA, MyuGammaA);
//     myu.close();
        this->MyuGammaA.load("fl_Work/fl_Vct_Myu" + TlUtils::xtos(niteration));

//     int totalnA;
//     Fl_Vct_Nyu   nyu( niteration );
//     nyu.open( "fl_Work", nyu.getfilename(), "read" );
//     nyu.getelemnum( &totalnA );
//     nyu.read( totalnA, NyuGammaA );
//     nyu.close();
        this->NyuGammaA.load("fl_Work/fl_Vct_Nyu" + TlUtils::xtos(niteration));

//     if (this->outlevel < -3){
//       Log << "*** PRINT OUT nPpq in readPpq ***\n";
//       for (int i = 0; i < nPpqA; ++i) {
//  Log << "row[" << i << "] = " << rowA[i] << ",  ";
//  Log << "col[" << i << "] = " << colA[i] << ",  ";
//  Log << "Ppq[" << i << "] = " << PpqA[i] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT MyuGamma in readPpq ***\n";
//       for (int j = 0; j < totalmA; ++j){
//  Log << "MyuGamma[" << j << "] = " << MyuGammaA[j] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT NyuGamma in readPpq ***\n";
//       for (int k = 0; k < totalnA; ++k){
//  Log << "NyuGamma[" << k << "] = " << NyuGammaA[k] << "\n";
//       }
//       Log << "\n";
//     }
    } else {
//     Fl_Mtr_PpqD  ppqA(RSFD,niteration, "a", "r");
//     ppqA.open("fl_Work", ppqA.getfilename(), "read");
//     ppqA.getelmnum(&nPpqA);
//     ppqA.readc(nPpqA, rowA, colA, PpqA);
//     ppqA.close();
        this->PpqA.load("fl_Work/fl_Mtr_PpqDra" + TlUtils::xtos(niteration));

//     Fl_Mtr_PpqD  ppqB(RSFD,niteration, "b", "r");
//     ppqB.open("fl_Work", ppqB.getfilename(), "read");
//     ppqB.getelmnum( &nPpqB );
//     ppqB.readc( nPpqB, rowB, colB, PpqB );
//     ppqB.close();
        this->PpqB.load("fl_Work/fl_Mtr_PpqDrb" + TlUtils::xtos(niteration));

//     int totalmA;
//     Fl_Vct_Myu   myuA( niteration,"a" );
//     myuA.open( "fl_Work", myuA.getfilename(), "read" );
//     myuA.getelemnum( &totalmA );
//     myuA.read( totalmA, MyuGammaA );
//     myuA.close();
        this->MyuGammaA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(niteration));

//     int totalmB;
//     Fl_Vct_Myu   myuB( niteration,"b" );
//     myuB.open( "fl_Work", myuB.getfilename(), "read" );
//     myuB.getelemnum( &totalmB );
//     myuB.read( totalmB, MyuGammaB );
//     myuB.close();
        this->MyuGammaB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(niteration));

//     int totalnA;
//     Fl_Vct_Nyu   nyuA( niteration, "a" );
//     nyuA.open( "fl_Work", nyuA.getfilename(), "read" );
//     nyuA.getelemnum( &totalnA );
//     nyuA.read( totalnA, NyuGammaA );
//     nyuA.close();
        this->NyuGammaA.load("fl_Work/fl_Vct_Nyua" + TlUtils::xtos(niteration));

//     int totalnB;
//     Fl_Vct_Nyu   nyuB( niteration, "b" );
//     nyuB.open( "fl_Work", nyuB.getfilename(), "read" );
//     nyuB.getelemnum( &totalnB );
//     nyuB.read( totalnB, NyuGammaB );
//     nyuB.close();
        this->NyuGammaB.load("fl_Work/fl_Vct_Nyub" + TlUtils::xtos(niteration));

//     if (this->outlevel < -3){
//       Log << "*** PRINT OUT nPpqA in readPpq ***\n";
//       for (int i = 0; i < nPpqA; ++i){
//  Log << "rowA[" << i << "] = " << rowA[i] << ",  ";
//  Log << "colA[" << i << "] = " << colA[i] << ",  ";
//  Log << "PpqA[" << i << "] = " << PpqA[i] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT nPpqB in readPpq ***\n";
//       for (int j = 0; j < nPpqB; ++j){
//  Log << "rowB[" << j << "] = " << rowB[j] << ",  ";
//  Log << "colB[" << j << "] = " << colB[j] << ",  ";
//  Log << "PpqB[" << j << "] = " << PpqB[j] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT MyuGammaA in readPpq ***\n";
//       for (int k = 0; k < totalmA; ++k){
//  Log << "MyuGammaA[" << k << "] = " << MyuGammaA[k] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT MyuGammaB in readPpq ***\n";
//       for (int l = 0; l < totalmB; ++l){
//  Log << "MyuGammaB[" << l << "] = " << MyuGammaB[l] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT NyuGammaA in readPpq ***\n";
//       for (int m = 0; m < totalnA; ++m){
//  Log << "NyuGammaA[" << m << "] = " << NyuGammaA[m] << "\n";
//       }
//       Log << "\n";

//       Log << "*** PRINT OUT NyuGammaB in readPpq ***\n";
//       for (int n = 0; n < totalnB; ++n){
//  Log << "NyuGammaB[" << n << "] = " << NyuGammaB[n] << "\n";
//       }
//       Log << "\n";
//     }
    }

    return 0;
}

// Read Sgd matrix.
// Because Sgd is a sparse matrix, all elements can be put in the main memory
// and be read at one time, though index g and d is about 2.5 times greater
// than p or q.
int DfXcpotfitting::readSgd()
{
//   Fl_Mtr_Sgdt sgd(RSFD);
//   sgd.open( "fl_Work", sgd.getfilename(), "read" );
//   sgd.getelmnum( &nSgd );
//   sgd.readc( nSgd, rowS, colS, Sgd );
//   sgd.close();

    this->Sgd.load("fl_Work/fl_Mtr_Sgdt");

//   if (this->outlevel < -5){
//     Log << "*** PRINT OUT Sgd in readSgd ***\n";
//     for (int i = 0; i < nSgd; ++i) {
//       Log << "row[" << i << "] = " << rowS[i] << ",  ";
//       Log << "col[" << i << "] = " << colS[i] << ",  ";
//       Log << "Sgd[" << i << "] = " << Sgd[i]  << "\n";
//     }
//     Log << "\n";
//   }

    return 0;
}

// This is a temporal function.
// The following case must be added.
// 1. NSP or SP         -> clear
// 2. Disk or Direct
// 3. Workfull or not
// 4. 4-divided inversed matrix routine
// Now we have Disk and not Workfull case.
int DfXcpotfitting::newton(int number)
{
    TlLogX& Log = TlLogX::getInstance();
    int   naux2 = 2 * naux;

    // Prepare the densPpq only once

    if (number == 1) {
        for (int i = 0; i < TotalPpq; ++i) {
            densPpqA[i] = 0.0;
        }

//     for (int i = 0; i < nPpqA; ++i) {
//       const int thrindex = ( rowA[i] + 1 ) * rowA[i] / 2 + colA[i];
//       if (rowA[i] == colA[i]){
//  PpqA[i] *= 0.5; // same as /= 2.0
//       }
//       densPpqA[thrindex] = PpqA[i];
//     }
        for (int i = 0; i < this->PpqA.getNumOfRows(); i++) {
            for (int j = 0; j <= i; i++) {
                const int thrindex = (i +1) * i / 2 + j;
                if (i == j) {
                    PpqA(i, j) *= 0.5;
                }
                densPpqA[thrindex] = PpqA(i, j);
            }
        }

    }

    // make Jacobian ( -J )
    // calc gSigma and fDelta
    for (int i = 0; i < naux2; ++i) {
        for (int j = 0; j < naux2; ++j) {
            Jacobi[i][j]    = 0.0;
            InvJacobi[i][j] = 0.0;
        }
    }

    for (int i = 0; i < naux; ++i) {
        gSigmaA[i] = fDeltaA[i] = 0.0;
    }

    // A1 Matrix
//   for (int k = 0; k < nSgd; ++k) {
//     Jacobi[ rowS[k] ][ colS[k] ] = Sgd[k];
//     Jacobi[ colS[k] ][ rowS[k] ] = Sgd[k];
//   }
    for (int i = 0; i < this->Sgd.getNumOfRows(); i++) {
        for (int j = 0; j <= i; i++) {
            Jacobi[i][j] = Sgd(i, j);
            Jacobi[j][i] = Sgd(i, j);
        }
    }

    // calc fDelta
    // Sum Ppq [ pqGamma ] is calculated once
    // Note that the contents of tAlpha are kept during the arival of newton.
    // tAlpha = Sum PpqA [ pqGamma ]

    if (number == 1) {
        for (int i = 0; i < naux; ++i) {
            tAlpha[i] = 0.0;
        }

        Fl_Int_Pqg fip;
        //fip.open( "fl_Temp", fip.getfilename(), "read" );

        int icont;
        do {
            int npqG;
            int nG;

            fip.read(&icont, &npqG, &nG, index0, n12, index1, index2, Intval);

            for (int i = 0; i < MaxEleNum; ++i) {
                workvector1[i] = 0;
            }

            spmtrprd(densPpqA, npqG, nG, index0, n12, index1, index2, Intval, naux, workvector1, tmpvectorA);

            for (int j = 0; j < naux; ++j) {
                tAlpha[j] += tmpvectorA[j];
            }
        } while (icont);
        //fip.close();
    }

    for (int j = 0; j < naux; ++j) {
        fDeltaA[j] = tAlpha[j];
    }

    // Read fl_Int_Gds one time
    Fl_Int_Gds fig;
    //fig.open( "fl_Temp", fig.getfilename(), "read" );

    int icont;
    do {
        int nGDS;
        int nS;

        fig.read(&icont, &nGDS, &nS, index0, n12, index1, index2, Intval);

        // A2 & A3 Matrix
        putA23(MyuGammaA, nGDS, nS, index0, n12, index1, index2, Intval, naux, Jacobi);

        // A4 Matrix
        putA4(NyuGammaA, nGDS, nS, index0, n12, index1, index2, Intval, naux, Jacobi);

        // calc gSigma
        for (int k = 0; k < naux; ++k) {
            tmpvectorA[k] = tmpvectorB[k] = 0.0;
        }

        spvctprd(MyuGammaA, nGDS, nS, index0, n12, index1, index2, Intval, naux, tmpvectorA);

        for (int j = 0; j < naux; ++j) {
            gSigmaA[j] += tmpvectorA[j];
        }

        // calc fDelta
        spvctprd3(this->MyuGammaA, this->NyuGammaA, nGDS, nS, index0, n12, index1, index2, Intval, naux, tmpvectorB);

        for (int k = 0; k < naux; ++k) {
            fDeltaA[k] -= tmpvectorB[k];
        }

        if (this->outlevel < -15) {
            Log << "*** PRINT OUT GDSigma in newton ***" << "\n";
            Log << "icont = " << icont << "\n";
            Log << "nGDS  = " << nGDS   << "\n";
            Log << "nS    = " << nS     << "\n";
            Log << "\n";
            for (int i = 0; i < nS; ++i) {
                Log << "indexS [" << i << "] = " << index0[i]  << "\n";
            }
            Log << "\n";
            for (int j = 0; j < nS; ++j) {
                Log << "nGD    [" << j << "] = " << n12[j]     << "\n";
            }
            Log << "\n";
            for (int k = 0; k < nGDS; ++k) {
                Log << "indexG [" << k << "] = " << index1[k]  << "\n";
            }
            for (int m = 0; m < nGDS; ++m) {
                Log << "indexD [" << m << "] = " << index2[m]  << "\n";
            }
            Log << "\n";
            for (int l = 0; l < nGDS; ++l) {
                Log << "GDS    [" << l << "] = " << Intval[l]  << "\n";
            }
            Log << "\n";
        }

    } while (icont);
    //fig.close();

    // keep tmpvector for calcE2
    for (int j = 0; j < naux; ++j) {
        tmpE2A[j] = gSigmaA[j];
    }

    // calc gSigma (continue)
    for (int j = 0; j < naux; ++j) {
        tmpvectorA[j] = 0.0;
    }

    spvctprd4(NyuGammaA, Sgd, naux, tmpvectorA);

    for (int j = 0; j < naux; ++j) {
        gSigmaA[j] -= tmpvectorA[j];
    }
    // calc gSigma end

    if (this->outlevel < -10) {
        Log << "*** PRINT OUT JacobianA in newton ***" << "\n";
        Log << "*** A1 Part ***" << "\n";
        for (int i = 0; i < naux; ++i) {
            Log << "J[" << i << "][j] =" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << " " << Jacobi[i][j];
                if ((j+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }
        Log << "*** A2 Part ***" << "\n";

        for (int i2 = naux; i2 < naux2; ++i2) {
            Log << "J[" << i2 << "][j] =" << "\n";
            for (int j2 = 0; j2 < naux; ++j2) {
                Log << " " << Jacobi[i2][j2];
                if ((j2+1)%6 == 0) Log << "\n";
            }
            Log << "\n";
        }

        Log << "*** A3 Part ***" << "\n";
        for (int i3 = 0; i3 < naux; ++i3) {
            Log << "J[" << i3 << "][j] =" << "\n";
            for (int j3 = naux; j3 < naux2; ++j3) {
                Log << " " << Jacobi[i3][j3];
                if ((j3+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }

        Log << "*** A4 Part ***" << "\n";
        for (int i4 = naux; i4 < naux2; ++i4) {
            Log << "J[" << i4 << "][j] =" << "\n";
            for (int j4 = naux; j4 < naux2; ++j4) {
                Log << " " << Jacobi[i4][j4];
                if ((j4+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }
        Log << "\n";
    }

    if (this->outlevel < -4) {
        Log << "*** PRINT OUT gSigmaA in newton ***" << "\n";
        for (int i = 0; i < naux; ++i) {
            Log << "gSigmaA[" << i << "] = " << gSigmaA[i] << "\n";
        }
        Log << "\n";
        Log << "*** PRINT OUT fDeltaA in newton ***" << "\n";
        for (int j = 0; j < naux; ++j) {
            Log << "fDeltaA[" << j << "] = " << fDeltaA[j] << "\n";
        }
        Log << "\n";
    }

    // make inversed  Jacobian
    for (int ii = 0; ii < naux2; ++ii) {
        tmpvector2[ii] = 0.0;
        indx[ii]       = 0;
    }

    double d;
    LUdcmp(Jacobi, naux2, indx, &d, tmpvector2);

    for (int j = 0; j < naux2; ++j) {
        for (int i = 0; i < naux2; ++i) {
            colinv[i] = 0.0;
        }
        colinv[j] = 1.0;
        LUbksb(Jacobi, naux2, indx, colinv);
        for (int k = 0; k < naux2; ++k) {
            InvJacobi[k][j] = colinv[k];
        }
    }

    if (this->outlevel < -10) {
        Log << "*** PRINT OUT Inversed JacobianA in newton ***" << "\n";
        Log << "*** X1 Part ***" << "\n";
        for (int i = 0; i < naux; ++i) {
            Log << "J-1[" << i << "][j] =" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << " " << InvJacobi[i][j];
                if ((j+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }

        Log << "*** X2 Part ***" << "\n";
        for (int i2 = naux; i2 < naux2; ++i2) {
            Log << "J-1[" << i2 << "][j] =" << "\n";
            for (int j2 = 0; j2 < naux; ++j2) {
                Log << " " << InvJacobi[i2][j2];
                if ((j2+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }

        Log << "*** X3 Part ***" << "\n";
        for (int i3 = 0; i3 < naux; ++i3) {
            Log << "J-1[" << i3 << "][j] =" << "\n";
            for (int j3 = naux; j3 < naux2; ++j3) {
                Log << " " << InvJacobi[i3][j3];
                if ((j3+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }

        Log << "*** X4 Part ***" << "\n";
        for (int i4 = naux; i4 < naux2; ++i4) {
            Log << "J-1[" << i4 << "][j] =" << "\n";
            for (int j4 = naux; j4 < naux2; ++j4) {
                Log << " " << InvJacobi[i4][j4];
                if ((j4+1) % 6 == 0) {
                    Log << "\n";
                }
            }
            Log << "\n";
        }
        Log << "\n";
    }

    // calc deltaMyu and deltaNyu
    for (int i = 0; i < naux; ++i) {
        X[i]      = gSigmaA[i];
        X[i+naux] = fDeltaA[i];
    }

    for (int j = 0; j < naux2; ++j) {
        deltaX[j] = 0.0;
        for (int k = 0; k < naux2; ++k) {
            deltaX[j] += InvJacobi[j][k] * X[k];
        }
    }

    for (int j = 0; j < naux2; ++j) {
        deltaX[j] *= DFactor;
    }

    double stotal = 0.0;
    for (int in = 0; in < naux2; ++in) {
        stotal += deltaX[in] * deltaX[in];
    }

    stotal = sqrt(stotal);

    double fukuescaleval = 2.0 * sqrt((double)naux2);
    if (stotal > fukuescaleval) {
        for (int i = 0; i < naux2; ++i) {
            deltaX[i] = fukuescaleval * deltaX[i] / stotal;
        }
    }

    if (this->outlevel < -4) {
        Log << "*** PRINT OUT deltaMyuA in newton ***" << "\n";
        for (int i = 0; i < naux; ++i) {
            Log << "deltaMyuA[" << i << "] = " << deltaX[i+naux] << "\n";
        }
        Log << "\n";
        Log << "*** PRINT OUT deltaNyuA in newton ***" << "\n";
        for (int j = 0; j < naux; ++j) {
            Log << "deltaNyuA[" << j << "] = " << deltaX[j] << "\n";
        }
        Log << "\n";
    }

    // renew MyuGamma and NyuGamma
    for (int il = 0; il < naux; ++il) {
        NyuGammaA[il] += deltaX[il];
        MyuGammaA[il] += deltaX[il+naux];
    }

    if (this->outlevel < -3) {
        Log << "*** PRINT OUT MyuGammaA in newton ***" << "\n";
        for (int i = 0; i < naux; ++i) {
            Log << "MyuGammaA[" << i << "] = " << MyuGammaA[i] << "\n";
        }
        Log << "\n";
        Log << "*** PRINT OUT NyuGammaA in newton ***" << "\n";
        for (int j = 0; j < naux; ++j) {
            Log << "NyuGammaA[" << j << "] = " << NyuGammaA[j] << "\n";
        }
        Log << "\n";
    }

    double dtotal = 0.0;
    for (int im = 0; im < naux2; ++im) {
        dtotal += deltaX[im] * deltaX[im];
    }

    if (this->scftype == "nsp") {
        // Non-spin-polarized case
        dtotal /= (DFactor * DFactor);
        dtotal  = sqrt(dtotal / naux2);

        if (outlevel <=-1) {
            Log << TlUtils::format("  dX(stotal)=%10.4le  update-dX(dtotal)=%10.4le", stotal, dtotal);
            Log << TlUtils::format("  f-s-v=%7.4lf DFactor=%4.2lf naux2=%4ld\n", fukuescaleval, DFactor, naux2);
            Log.flush();
        }

        if (dtotal < NREPS) {
            return 0;
        } else {
            return -1;
        }
    } else {
        // Spin-polarized case

        // Prepare the densPpq only once
        if (number == 1) {
            for (int i = 0; i < TotalPpq; ++i) {
                densPpqB[i] = 0.0;
            }

//       for (int i = 0; i < nPpqB; ++i) {
//  int thrindex = ( rowB[i] + 1 ) * rowB[i] / 2 + colB[i];
//  if (rowB[i] == colB[i]){
//    PpqB[i] /= 2.0;
//  }
//  densPpqB[thrindex] = PpqB[i];
//       }
            for (int i = 0; i < this->PpqB.getNumOfRows(); i++) {
                for (int j = 0; j <= i; j++) {
                    const int thrindex = (i +1) * i / 2 + j;
                    if (i == j) {
                        PpqB(i, j) *= 0.5;
                    }
                    densPpqB[thrindex] = this->PpqB(i, j);
                }
            }

        }

        // make Jacobian ( -J )
        // calc gSigma and fDelta
        for (int i = 0; i < naux2; ++i) {
            for (int j = 0; j < naux2; ++j) {
                Jacobi[i][j]    = 0.0;
                InvJacobi[i][j] = 0.0;
            }
        }

        for (int i = 0; i < naux; ++i) {
            gSigmaB[i] = fDeltaB[i] = 0.0;
        }

        // A1 Matrix
        // This procedure is not needed.
//     for (int k = 0; k < nSgd; ++k) {
//       Jacobi[ rowS[k] ][ colS[k] ] = Sgd[k];
//       Jacobi[ colS[k] ][ rowS[k] ] = Sgd[k];
//     }
        for (int i = 0; i < this->Sgd.getNumOfRows(); i++) {
            for (int j = 0; i <= j; i++) {
                Jacobi[i][j] = Sgd(i, j);
                Jacobi[j][i] = Sgd(i, j);
            }
        }

        // calc fDelta
        // Sum Ppq [ pqGamma ] is calculated once
        // Note that the contents of tBeta are kept during the arival of newton.
        // tBeta = Sum PpqB [ pqGamma ]

        if (number == 1) {
            for (int i = 0; i < naux; ++i) {
                tBeta[i] = 0.0;
            }

            Fl_Int_Pqg fip;
            //fip.open( "fl_Temp", fip.getfilename(), "read" );

            do {
                int icont;
                int npqG;
                int nG;
                fip.read(&icont, &npqG, &nG, index0, n12, index1, index2, Intval);

                for (int i = 0; i < MaxEleNum; ++i) {
                    workvector2[i] = 0;
                }

                spmtrprd(densPpqB, npqG, nG, index0, n12, index1, index2, Intval, naux, workvector2, tmpvectorB);

                for (int j = 0; j < naux; ++j) {
                    tBeta[j] += tmpvectorB[j];
                }
            } while (icont);
        }

        for (int j = 0; j < naux; ++j) {
            fDeltaB[j] = tBeta[j];
        }

        // Read fl_Int_Gds one time
        Fl_Int_Gds fig;
        //fig.open( "fl_Temp", fig.getfilename(), "read" );

        do {
            int icont;
            int nGDS;
            int nS;
            fig.read(&icont, &nGDS, &nS, index0, n12, index1, index2, Intval);

            // A2 & A3 Matrix
            putA23(MyuGammaB, nGDS, nS, index0, n12, index1, index2, Intval, naux, Jacobi);

            // A4 Matrix
            putA4(NyuGammaB, nGDS, nS, index0, n12, index1, index2, Intval, naux, Jacobi);

            // calc gSigma
            for (int k = 0; k < naux; ++k) {
                tmpvectorA[k] = tmpvectorB[k] = 0.0;
            }

            spvctprd(MyuGammaB, nGDS, nS, index0, n12, index1, index2, Intval, naux, tmpvectorA);

            for (int j = 0; j < naux; ++j) {
                gSigmaB[j] += tmpvectorA[j];
            }

            // calc fDelta
            spvctprd3(MyuGammaB, NyuGammaB, nGDS, nS, index0, n12, index1, index2, Intval, naux, tmpvectorB);

            for (int k = 0; k < naux; ++k) {
                fDeltaB[k] -= tmpvectorB[k];
            }
        } while (icont);
        //fig.close();

        // keep tmpvector for calcE2
        for (int j = 0; j < naux; ++j) {
            tmpE2B[j] = gSigmaB[j];
        }

        // calc gSigma (continue)
        for (int j = 0; j < naux; ++j) {
            tmpvectorB[j] = 0.0;
        }

        spvctprd4(NyuGammaB, Sgd, naux, tmpvectorB);

        for (int j = 0; j < naux; ++j) {
            gSigmaB[j] -= tmpvectorB[j];
        }
        // calc gSigma end

        if (this->outlevel < -10) {
            Log << "*** PRINT OUT JacobianB in newton ***" << "\n";
            Log << "*** A1 Part ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "J[" << i << "][j] =" << "\n";
                for (int j = 0; j < naux; ++j) {
                    Log << " " << Jacobi[i][j];
                }
                Log << "\n";
            }

            Log << "*** A2 Part ***" << "\n";
            for (int i2 = naux; i2 < naux2; ++i2) {
                Log << "J[" << i2 << "][j] =" << "\n";
                for (int j2 = 0; j2 < naux; ++j2) {
                    Log << " " << Jacobi[i2][j2];
                }
                Log << "\n";
            }

            Log << "*** A3 Part ***" << "\n";
            for (int i3 = 0; i3 < naux; ++i3) {
                Log << "J[" << i3 << "][j] =" << "\n";
                for (int j3 = naux; j3 < naux2; ++j3)
                    Log << " " << Jacobi[i3][j3];
                Log << "\n";
            }

            Log << "*** A4 Part ***" << "\n";
            for (int i4 = naux; i4 < naux2; ++i4) {
                Log << "J[" << i4 << "][j] =" << "\n";
                for (int j4 = naux; j4 < naux2; ++j4) {
                    Log << " " << Jacobi[i4][j4];
                }
                Log << "\n";
            }
            Log << "\n";
        }

        // make inversed  Jacobian
        for (int ii = 0; ii < naux2; ++ii) {
            tmpvector2[ii] = 0.0;
            indx[ii]       = 0;
        }

        double d;
        LUdcmp(Jacobi, naux2, indx, &d, tmpvector2);

        for (int j = 0; j < naux2; ++j) {
            for (int i = 0; i < naux2; ++i) {
                colinv[i] = 0.0;
            }
            colinv[j] = 1.0;
            LUbksb(Jacobi, naux2, indx, colinv);
            for (int k = 0; k < naux2; ++k) {
                InvJacobi[k][j] = colinv[k];
            }
        }

        if (this->outlevel < -10) {
            Log << "*** PRINT OUT Inversed JacobianB in newton ***" << "\n";
            Log << "*** X1 Part ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "J-1[" << i << "][j] =" << "\n";
                for (int j = 0; j < naux; ++j) {
                    Log << " " << InvJacobi[i][j];
                }
                Log << "\n";
            }

            Log << "*** X2 Part ***" << "\n";
            for (int i2 = naux; i2 < naux2; ++i2) {
                Log << "J-1[" << i2 << "][j] =" << "\n";
                for (int j2 = 0; j2 < naux; ++j2) {
                    Log << " " << InvJacobi[i2][j2];
                }
                Log << "\n";
            }

            Log << "*** X3 Part ***" << "\n";
            for (int i3 = 0; i3 < naux; ++i3) {
                Log << "J-1[" << i3 << "][j] =" << "\n";
                for (int j3 = naux; j3 < naux2; ++j3) {
                    Log << " " << InvJacobi[i3][j3];
                }
                Log << "\n";
            }

            Log << "*** X4 Part ***" << "\n";
            for (int i4 = naux; i4 < naux2; ++i4) {
                Log << "J-1[" << i4 << "][j] =" << "\n";
                for (int j4 = naux; j4 < naux2; ++j4) {
                    Log << " " << InvJacobi[i4][j4];
                }
                Log << "\n";
            }
            Log << "\n";
        }

        if (this->outlevel < -4) {
            Log << "*** PRINT OUT gSigmaB in newton ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "gSigmaB[" << i << "] = " << gSigmaB[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT fDeltaB in newton ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "fDeltaB[" << j << "] = " << fDeltaB[j] << "\n";
            }
            Log << "\n";
        }

        // calc deltaMyu and deltaNyu
        for (int i = 0; i < naux; ++i) {
            X[i]      = gSigmaB[i];
            X[i+naux] = fDeltaB[i];
        }

        for (int j = 0; j < naux2; ++j) {
            deltaX[j] = 0.0;
            for (int k = 0; k < naux2; ++k) {
                deltaX[j] += InvJacobi[j][k] * X[k];
            }
        }

        for (int j = 0; j < naux2; ++j) {
            deltaX[j] *= DFactor;
        }

        stotal = 0.0;
        for (int in = 0; in < naux2; ++in) {
            stotal += deltaX[in] * deltaX[in];
        }

        stotal = sqrt(stotal);

        if (this->outlevel < -2) {
            Log << "*** PRINT OUT stotal in newton ***" << "\n";
            Log << "stotal = " << stotal << "\n";
            Log << "\n";
        }

        double  fukuescaleval = 2.0 * sqrt((double)naux2) ;
        if (stotal > fukuescaleval) {
            for (int i = 0; i < naux2; ++i) {
                deltaX[i] = fukuescaleval * deltaX[i] / stotal;
            }
        }

        if (this->outlevel < -4) {
            Log << "*** PRINT OUT deltaMyuB in newton ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "deltaMyuB[" << i << "] = " << deltaX[i+naux] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT deltaNyuB in newton ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "deltaNyuB[" << j << "] = " << deltaX[j] << "\n";
            }
            Log << "\n";
        }

        // renew MyuGamma and NyuGamma
        for (int il = 0; il < naux; ++il) {
            NyuGammaB[il] += deltaX[il];
            MyuGammaB[il] += deltaX[il+naux];
        }

        if (this->outlevel < -3) {
            Log << "*** PRINT OUT MyuGammaB in newton ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGammaB[" << i << "] = " << MyuGammaB[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT NyuGammaB in newton ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "NyuGammaB[" << j << "] = " << NyuGammaB[j] << "\n";
            }
            Log << "\n";
        }

        for (int im = 0; im < naux2; ++im) {
            dtotal += deltaX[im] * deltaX[im];
        }

        dtotal /= (DFactor * DFactor);
        dtotal = sqrt(dtotal / naux2);

        if (this->outlevel < -2) {
            Log << "*** PRINT OUT dtotal in newton ***" << "\n";
            Log << "dtotal = " << dtotal << "\n";
            Log << "\n";
        }

        if (dtotal < 2*NREPS) {
            return 0;
        } else {
            return -1;
        }
    }             // end of Spin-pol if
}

// Calculate E1
// the following procedures are already done at calc-fDelta in newton
// 1. put pqGamma into each workarea
// 2. calclate partial E1
// 3. goto 1 while end
// 4. delete Ppq memory
int DfXcpotfitting::calcE1()
{
    TlLogX& Log = TlLogX::getInstance();
    E1A = E1B = 0.0;

    if (this->scftype == "nsp") {
        for (int i = 0; i < naux; ++i) {
            E1A += MyuGammaA[i] * tAlpha[i];
        }

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E1 in calcE1 ***" << "\n";
            Log << "E1 = " << E1A << "\n";
            Log << "\n";
        }

    } else {
        for (int i = 0; i < naux; ++i) {
            E1A += MyuGammaA[i] * tAlpha[i];
            E1B += MyuGammaB[i] * tBeta[i];
        }

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E1A in calcE1 ***" << "\n";
            Log << "E1A = " << E1A << "\n";
            Log << "\n";
            Log << "*** PRINT OUT E1B in calcE1 ***" << "\n";
            Log << "E1B = " << E1B << "\n";
            Log << "\n";
        }
    }

    return 0;
}

// Calculate E2
// we already calculated tmpE2 with spvctprd
int DfXcpotfitting::calcE2()
{
    TlLogX& Log = TlLogX::getInstance();
    E2A = E2B = 0.0;

    if (this->scftype == "nsp") {
        for (int i = 0; i < naux; ++i) {
            E2A += NyuGammaA[i] * tmpE2A[i];
        }

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E2 in calcE2 ***" << "\n";
            Log << "E2 = " << E2A << "\n";
            Log << "\n";
        }

    } else {
        for (int i = 0; i < naux; ++i) {
            E2A += NyuGammaA[i] * tmpE2A[i];
            E2B += NyuGammaB[i] * tmpE2B[i];
        }

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E2A in calcE2 ***" << "\n";
            Log << "E2A = " << E2A << "\n";
            Log << "\n";
            Log << "*** PRINT OUT E2B in calcE2 ***" << "\n";
            Log << "E2B = " << E2B << "\n";
            Log << "\n";
        }
    }

    return 0;
}


// Calculate E3
// 1. put Sgd into workarea at one time
// 2. catlculate E3 from NyuGamma and Sgd
int DfXcpotfitting::calcE3()
{
    TlLogX& Log = TlLogX::getInstance();
    E3A = E3B = 0.0;

    if (this->scftype == "nsp") {
        // Non-spin-polarized case

        E3A = spvctprd2(NyuGammaA, Sgd);

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E3 in calcE3 ***" << "\n";
            Log << "E3 = " << E3A << "\n";
            Log << "\n";
        }
    } else {
        // Spin-polarized case

        E3A = spvctprd2(NyuGammaA, Sgd);
        E3B = spvctprd2(NyuGammaB, Sgd);

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT E3A in calcE3 ***" << "\n";
            Log << "E3A = " << E3A << "\n";
            Log << "\n";
            Log << "*** PRINT OUT E3B in calcE3 ***" << "\n";
            Log << "E3B = " << E3B << "\n";
            Log << "\n";
        }
    }

    return 0;
}

// calculate Myu & Nyu
int DfXcpotfitting::calcMN()
{
    TlLogX& Log = TlLogX::getInstance();
    if (this->scftype == "nsp") {
        // Non-spin-polarized case

        SFactorMA = pow(E1A * E2A / (E3A * E3A), F13);
        SFactorNA = pow(E1A * E1A / (E2A * E3A), F13);

        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = SFactorMA * MyuGammaA[i];
            NyuGammaA[i] = SFactorNA * NyuGammaA[i];
        }

        if (this->outlevel < -2) {
            Log << "*** PRINT OUT SFactorM in calcMN ***" << "\n";
            Log << "SFactorM = " << SFactorMA << "\n";
            Log << "\n";
            Log << "*** PRINT OUT SFactorN in calcMN ***" << "\n";
            Log << "SFactorN = " << SFactorNA << "\n";
            Log << "\n";
        }

        if (this->outlevel < -3) {
            Log << "*** PRINT OUT MyuGamma in calcMN ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGamma[" << i << "] = " << MyuGammaA[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT NyuGamma in calcMN ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "NyuGamma[" << j << "] = " << NyuGammaA[j] << "\n";
            }
            Log << "\n";
        }
    } else {
        // Spin-polarized case

        SFactorMA = pow(E1A * E2A / (E3A * E3A), F13);
        SFactorNA = pow(E1A * E1A / (E2A * E3A), F13);
        SFactorMB = pow(E1B * E2B / (E3B * E3B), F13);
        SFactorNB = pow(E1B * E1B / (E2B * E3B), F13);

        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = SFactorMA * MyuGammaA[i];
            NyuGammaA[i] = SFactorNA * NyuGammaA[i];
            MyuGammaB[i] = SFactorMB * MyuGammaB[i];
            NyuGammaB[i] = SFactorNB * NyuGammaB[i];
        }

        if (this->outlevel < -2) {
            Log << "*** PRINT OUT SFactorMA in calcMN ***" << "\n";
            Log << "SFactorMA = " << SFactorMA << "\n";
            Log << "\n";
            Log << "*** PRINT OUT SFactorMB in calcMN ***" << "\n";
            Log << "SFactorMB = " << SFactorMB << "\n";
            Log << "\n";
            Log << "*** PRINT OUT SFactorNA in calcMN ***" << "\n";
            Log << "SFactorNA = " << SFactorNA << "\n";
            Log << "\n";
            Log << "*** PRINT OUT SFactorNB in calcMN ***" << "\n";
            Log << "SFactorNB = " << SFactorNB << "\n";
            Log << "\n";
        }

        if (this->outlevel < -3) {
            Log << "*** PRINT OUT MyuGammaA in calcMN ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGammaA[" << i << "] = " << MyuGammaA[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT NyuGammaA in calcMN ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "NyuGammaA[" << j << "] = " << NyuGammaA[j] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT MyuGammaB in calcMN ***" << "\n";
            for (int k = 0; k < naux; ++k) {
                Log << "MyuGammaB[" << k << "] = " << MyuGammaB[k] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT NyuGammaB in calcMN ***" << "\n";
            for (int l = 0; l < naux; ++l) {
                Log << "NyuGammaB[" << l << "] = " << NyuGammaB[l] << "\n";
            }
            Log << "\n";
        }
    }

    return 0;
}

int DfXcpotfitting::bfscalMyu()
{
    TlLogX& Log = TlLogX::getInstance();
    // Scale Myu before calculation
    if (this->scftype == "nsp") {
        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = MyuGammaA[i] / Alpha[i];
        }

        if (this->outlevel < -5) {
            Log << "*** PRINT OUT MyuGamma in bfscalMyu ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGamma[" << i << "] = " << MyuGammaA[i] << "\n";
            }
            Log << "\n";
        }
    } else {
        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = MyuGammaA[i] / Alpha[i];
            MyuGammaB[i] = MyuGammaB[i] / Alpha[i];
        }

        if (this->outlevel < -5) {
            Log << "*** PRINT OUT MyuGammaA in bfscalMyu ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGammaA[" << i << "] = " << MyuGammaA[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT MyuGammaB in bfscalMyu ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "MyuGammaB[" << j << "] = " << MyuGammaB[j] << "\n";
            }
            Log << "\n";
        }
    }

    return 0;
}

// Scale Myu after calculation
// File Log in fl_Vct_Myu and fl_Vct_Nyu
// if SCF type is non-spin-polarized case, file name is fl_Vct_Myu(int),
// while in spin-polarized case, file names are fl_Vct_Myua(int) and fl_Vct_Myub(int).
int DfXcpotfitting::afscalMyu()
{
    TlLogX& Log = TlLogX::getInstance();
    if (this->scftype == "nsp") {
        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = MyuGammaA[i] * Alpha[i];
        }

//     Fl_Vct_Myu myu( niteration1 );
//     myu.open( "fl_Work", myu.getfilename(), "write" );
//     myu.putelemnum( &naux );
//     myu.write( naux, MyuGammaA );
//     myu.close();
        this->MyuGammaA.save("fl_Work/fl_Vct_Myu" + TlUtils::xtos(niteration1));

//     Fl_Vct_Nyu   nyu( niteration1 );
//     nyu.open( "fl_Work", nyu.getfilename(), "write" );
//     nyu.putelemnum( &naux );
//     nyu.write( naux, NyuGammaA );
//     nyu.close();
        this->NyuGammaA.save("fl_Work/fl_Vct_Nyu" + TlUtils::xtos(niteration1));

        if (this->outlevel <= -2) {
            Log << "*** PRINT OUT MyuGamma in afscalMyu ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGamma[" << i << "] = " << MyuGammaA[i] << "\n";
            }
            Log << "\n";
        }
    } else {
        for (int i = 0; i < naux; ++i) {
            MyuGammaA[i] = MyuGammaA[i] * Alpha[i];
            MyuGammaB[i] = MyuGammaB[i] * Alpha[i];
        }

//     Fl_Vct_Myu   myua( niteration1, "a" );
//     myua.open( "fl_Work", myua.getfilename(), "write" );
//     myua.putelemnum( &naux );
//     myua.write( naux, MyuGammaA );
//     myua.close();
        this->MyuGammaA.save("fl_Work/fl_Vct_Myua" + TlUtils::xtos(niteration1));

//     Fl_Vct_Myu   myub( niteration1, "b" );
//     myub.open( "fl_Work", myub.getfilename(), "write" );
//     myub.putelemnum( &naux );
//     myub.write( naux, MyuGammaB );
//     myub.close();
        this->MyuGammaB.save("fl_Work/fl_Vct_Myub" + TlUtils::xtos(niteration1));

//     Fl_Vct_Nyu   nyua( niteration1, "a" );
//     nyua.open( "fl_Work", nyua.getfilename(), "write" );
//     nyua.putelemnum( &naux );
//     nyua.write( naux, NyuGammaA );
//     nyua.close();
        this->NyuGammaA.save("fl_Work/fl_Vct_Nyua" + TlUtils::xtos(niteration1));

//     Fl_Vct_Nyu   nyub( niteration1, "b" );
//     nyub.open( "fl_Work", nyub.getfilename(), "write" );
//     nyub.putelemnum( &naux );
//     nyub.write( naux, NyuGammaB );
//     nyub.close();
        this->NyuGammaB.save("fl_Work/fl_Vct_Nyub" + TlUtils::xtos(niteration1));

        if (this->outlevel < 0) {
            Log << "*** PRINT OUT MyuGammaA in afscalMyu ***" << "\n";
            for (int i = 0; i < naux; ++i) {
                Log << "MyuGammaA[" << i << "] = " << this->MyuGammaA[i] << "\n";
            }
            Log << "\n";
            Log << "*** PRINT OUT MyuGammaB in afscalMyu ***" << "\n";
            for (int j = 0; j < naux; ++j) {
                Log << "MyuGammaB[" << j << "] = " << this->MyuGammaB[j] << "\n";
            }
            Log << "\n";
        }
    }

    return 0;
}

void DfXcpotfitting::Error(const char* string)
{
    std::cout << "*** Error in DfXcpotfitting:" << std::endl;
    std::cout << string << std::endl;
    std::cout << std::endl;
}

void DfXcpotfitting::Print(int number, double* vector)
{
    TlLogX& Log = TlLogX::getInstance();
    for (int i = 0; i < number; ++i) {
        Log << "vector[" << i << "] = " << vector[i] << "\n";
    }
}

void DfXcpotfitting::Print(int number, int* row, int* col, double* matrix)
{
    TlLogX& Log = TlLogX::getInstance();
    for (int i = 0; i < number; ++i) {
        Log << "matrix[" << row[i] << ", " << col[i] << "] = " << matrix[i] << "\n";
    }
}

