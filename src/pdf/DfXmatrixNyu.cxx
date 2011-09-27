#include <fstream>
#include <string>
#include "DfXmatrixNyu.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"

DfXmatrixNyu::DfXmatrixNyu()
{
}


DfXmatrixNyu::~DfXmatrixNyu()
{
}


void DfXmatrixNyu::main()
{
    TlLogX& Log = TlLogX::getInstance();

//   Fl_Globalinput Gbi(">>>>SCF");
//   int  pout = Gbi.printlevel();

    // read S matrix
    TlSymmetricMatrix Sgd;
    {
        Sgd.load("fl_Work/fl_Mtr_Sgd.matrix");

//     if( -1>pout ){
//       Out << "@@@@ -1>pout\n";
//       Out << "Sgd (Snyu,nyu) matrix\n";
//       Sgd.print(Out);
//     }
    }

    // diagonalize S matrix
    {

        Log  << " @@@@ diagonalization with JAMOL4 routine\n";

        TlVector eigVal;
        TlMatrix eigVec;

        Sgd.diagonal(&eigVal, &eigVec);

//     if( -1>pout ){
//       Out << "@@@@ -1>pout\n";
//       Out << "eigenvector of Sgd matrix\n";
//       eigVec.print(Out);
//     }

        Log << "eigenvalue  of Sgd matrix\n";
        eigVal.print(Log);
    }

    return;
}
