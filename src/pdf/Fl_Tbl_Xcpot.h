#ifndef FL_TBL_XCPOT_H
#define FL_TBL_XCPOT_H

#include "Fl_Geometry.h"
#include "Fl_Gto.h"
#include "Fl_Gto_Xcpot.h"
#include "TlPrdctbl.h"

class Fl_Tbl_Xcpot {
private:
    static int  flag1;         // flag for Constructor.
    static int  flag2;         // flag for Destructor.
    static int cGtoTotalNum;  // s*1, p*3, d*5.

// ---------- Common variable in Fl_Tbl_Xcpot ----------

    enum { Maxfilename  =  80 };               // MaxNumber of FileName.
    enum { MaxAtomName  =  20 };               // MaxNumber of AtomName.
    enum { MaxLabelname =  20 };               // MaxNumber of Label.
    enum { MaxFlBasisTblmemory = 100000000 };  // toriaezu.

// ---------- Private Data member ----------

    int  AtomNum;         // Total Atom Number.
    int   AtomKindNum;     // Number of Atomtype by AtomKind Label2.
//  int  cGtoTotalNum;    // CGTO Total Number.
//   char* TblFile;         // Pointer to OutTableFile.
    // This Class make a Table : "XcpotTable" .

// ---------- Peculiar structure for an Atom. ----------

    struct BasisTable {     // This structure has a line_data in XcpotTable.
        int   cGtonum;       // Number in Fl_Gto_Xcpot.
//  int    contraction;   // Number of Contraction.
        int   InpAtomnumber; // Nucleus Number.
        char*  Atom;          // Atomtype.
        char*  shell;         // Shelltype { s , p , d }.
        char*  dumy1;         // dumy1. : details of shelltype.
        char*  dumy2;         // dumy2. : Label2.
    };

    BasisTable* BsTbl;      // Declaration for Structure_Object.

// ---------- Private Member Function ----------

    int   prepare();
    int  Estmemval();
    int   getMemory();
    int   makeTable();      // Table : XcpotTable .
    int   setData();        // SetFunction  to Private Data Member.


// ========== Public Member Function ==========

public:

    Fl_Tbl_Xcpot();
    ~Fl_Tbl_Xcpot();

    int    getcGtoTotalNum(void)     {
        return cGtoTotalNum;
    };

    int    getcgtonum(int num)      {
        return BsTbl[num].cGtonum;
    };
    int    getInpAtomnum(int num)   {
        return BsTbl[num].InpAtomnumber;
    };
    char*   getInpAtomtype(int num)  {
        return BsTbl[num].Atom;
    };
    char*   getShelltype(int num)    {
        return BsTbl[num].shell;
    };
    char*   getNote1(int num)        {
        return BsTbl[num].dumy1;
    };
    char*   getNote2(int num)        {
        return BsTbl[num].dumy2;
    };
};
#endif
/* end of file(Fl_Tbl_Xcpot.h) */
