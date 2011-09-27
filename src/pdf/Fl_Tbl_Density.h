#ifndef FL_TBL_DENSITY_H
#define FL_TBL_DENSITY_H

#include <vector>
#include <string>

#include "Fl_Geometry.h"
#include "Fl_Gto.h"
#include "Fl_Gto_Density.h"
#include "TlPrdctbl.h"

class Fl_Tbl_Density {
public:
    Fl_Tbl_Density();
    ~Fl_Tbl_Density();

    int getcGtoTotalNum() const {
        return this->BsTbl.size() -1;
    }

    int getcgtonum(int num) const {
        //std::cout << TlUtils::format("getcgtonum(%d) = %d", num, BsTbl[num].cGtonum) << std::endl;
        return BsTbl[num].cGtonum;
    }

    int getInpAtomnum(int num) const {
        return BsTbl[num].InpAtomnumber;
    }

    std::string getInpAtomtype(int num) const {
        return BsTbl[num].Atom;
    }

    char getShelltype(int num) const {
        return BsTbl[num].shell;
    }

    std::string getNote1(int num) const {
        return BsTbl[num].dumy1;
    }

    std::string getNote2(int num) const {
        return BsTbl[num].dumy2;
    }

private:
    int prepare();
    int makeTable();      // Table : DensityTable .
    int setData();        // SetFunction  to Private Data Member.

private:
    static int flag1;         // flag for Constructor.
    static int flag2;         // flag for Destructor.

    struct BasisTable {      // This structure has a line_data in DensityTable.
        int cGtonum;       // Number in Fl_Gto_Density.
        int InpAtomnumber; // Nucleus Number.
        std::string Atom;          // Atomtype.
        char shell;         // Shelltype { s , p , d }.
        std::string dumy1;         // dumy1. : details of shelltype.
        std::string dumy2;         // dumy2. : Label2.
    };

    std::vector<BasisTable> BsTbl;

};

#endif

