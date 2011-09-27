#include <fstream>
#include <string>
#include "DfXmatrixMyu.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"

DfXmatrixMyu::DfXmatrixMyu()
{
}


DfXmatrixMyu::~DfXmatrixMyu()
{
}


void DfXmatrixMyu::main()
{
    TlLogX& Log = TlLogX::getInstance();

    //Fl_Globalinput  Gbi(">>>>SCF");
    //int  pout = Gbi.printlevel();

    // read S matrix
    TlSymmetricMatrix Sgd2;
    {
        Sgd2.load("fl_Work/fl_Mtr_Sgd2.matrix");

//     if( -1>pout ){
//       Out << "@@@@ -1>pout\n";
//       Out << "Sgd2 (Smyu,myu) matrix\n";
//       Sgd2.print(Out);
//     }
    }

    // diagonalize S matrix
    {
        Log  << " @@@@ diagonalization with JAMOL4 routine\n";

        TlVector eigVal;
        TlMatrix eigVec;
        Sgd2.diagonal(&eigVal, &eigVec);

//     if( -1>pout ){
//       Out << "@@@@ -1>pout" << "\n";
//       Out << "eigenvector of Sgd2 matrix" << "\n";
//       eigVec.print(Out);
//     }
        Log << "eigenvalue  of Sgd2 matrix\n";
        eigVal.print(Log);
    }

    return;
}
