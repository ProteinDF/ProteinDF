#include <fstream>
#include <string>
#include "DfXmatrixRou.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"

DfXmatrixRou::DfXmatrixRou()
{
}


DfXmatrixRou::~DfXmatrixRou()
{
}


void DfXmatrixRou::main()
{
    TlLogX& Log = TlLogX::getInstance();

//   Fl_Globalinput Gbi(">>>>SCF");
//   int pout = Gbi.printlevel();

    // read S matrix
    TlSymmetricMatrix Sab2;
    {
        Sab2.load("fl_Work/fl_Mtr_Sab2.matrix");

//     if (-1 > pout){
//       Out << "@@@@ -1>pout\n";
//       Out << "Sab2 (Srou,rou) matrix\n";
//       Sab2.print(Out);
//     }
    }

    // diagonalize S matrix
    {
        Log << " @@@@ diagonalization with JAMOL4 routine\n";

        TlVector eigVal;
        TlMatrix eigVec;

        Sab2.diagonal(&eigVal, &eigVec);

//     if( -1>pout ){
//       Out << "@@@@ -1>pout\n";
//       Out << "eigenvector of Sab2 matrix\n";
//       eigVec.print(Out);
//     }

        Log << "eigenvalue  of Sab2 matrix\n";
        eigVal.print(Log);
    }

    return;
}
