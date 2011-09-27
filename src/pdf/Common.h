#ifndef COMMON_H
#define COMMON_H

// // The followings are the contents of "INTSIZ".
// //********************* PARAMETER ***************************************
// //     PARAMETER (MAXSHL=3000,MAXCON=200,MAXATM=800)
// //     PARAMETER (MAXSH2=5000,MAXCO2=400)
// //     PARAMETER (MAXNA=300,MXNPQA=4462)
// //     PARAMETER (IBLKSZ=40)
// //     PARAMETER (MAXDEN=9000)
// //     PARAMETER (MXORBF=5000,MXDENF=9000,MXXCP1=9000,MXXCP2=9000)
// //***********************************************************************
// // INTSIZ END

#define MAXNA           300             // max nA.
#define MAXNPQA         4462            // max npqA.

// #define MaxFNLen        15
// #define MaxAtmNum       2000
// #define MaxAuxNum       18000            // this statement is put into common.h
// // calclation can be carried out more than
// // MaxAuxNum.  if naux is less than MaxAuxNum,
// // SAlpahBeta_inv can be put into workarea
// // at one time.
// //#define MaxEleNum       6750000         // 4500 * 4500 /3
#define MaxEleNum       60500000        // 11000 * 11000 /2
// #define MaxEleAux       200000000        // 20000 * 20000 /2
// #define MaxMemory       2000000000       // set max memory is 2000 MB.


// #define F13             0.3333333333333333
// #define F23             0.6666666666666667
// #define F43             1.3333333333333333
// #define R3_2            1.2599210498948732
// #define PAI             3.14159265358979323846264338327905
// #define BOHR            0.52917706

// // for Parallel Computing

// // #define KSHELL_BLOCK 20
// #define KSHELL_BLOCK 4

#endif // COMMON_H
