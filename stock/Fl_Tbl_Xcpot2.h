#ifndef FL_TBL_XCPOT2_H
#define FL_TBL_XCPOT2_H

#include "Fl_Geometry.h"
#include "Fl_Gto.h"
#include "Fl_Gto_Xcpot2.h"
#include "TlPrdctbl.h"
#include "Fl_Geometry.h"

class Fl_Tbl_Xcpot2 {

public:
    Fl_Tbl_Xcpot2(const Fl_Geometry& flGeom);
    ~Fl_Tbl_Xcpot2();

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

private:
    const Fl_Geometry& flGeom_;
    static int  flag1;          // flag for Constructor.
    static int  flag2;          // flag for Destructor.
    static int cGtoTotalNum;   // s*1, p*3, d*5.

// ---------- Common variable in Fl_Tbl_Xcpot2 ----------

    enum { Maxfilename  =  80 };               // MaxNumber of FileName.
    enum { MaxAtomName  =  20 };               // MaxNumber of AtomName.
    enum { MaxLabelname =  20 };               // MaxNumber of Label.
    enum { MaxFlBasisTblmemory = 100000000 };  // toriaezu.

// ---------- Private Data member ----------

    int  AtomNum;         // Total Atom Number.
    int   AtomKindNum;     // Number of Atomtype by AtomKind Label2.
//  int  cGtoTotalNum;    // CGTO Total Number.
    std::string TblFile;         // Pointer to OutTableFile.
    // This Class make a Table : "Xcpot2Table" .

// ---------- Peculiar structure for an Atom. ----------

    struct BasisTable {     // This structure has a line_data in Xcpot2Table.
        int   cGtonum;       // Number in Fl_Gto_Xcpot2.
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
    int   makeTable();      // Table : Xcpot2Table .
    int   setData();        // SetFunction  to Private Data Member.


// ========== Public Member Function ==========


};
#endif
/* end of file(Fl_Tbl_Xcpot2.h) */
