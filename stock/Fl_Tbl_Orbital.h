#ifndef FL_TBL_ORBITAL_H
#define FL_TBL_ORBITAL_H

#include <vector>
#include <string>

#include "Fl_Geometry.h"
#include "Fl_Gto.h"
//#include "Fl_Gto_Orbital.h"
#include "TlPrdctbl.h"

class Fl_Tbl_Orbital {
public:
    Fl_Tbl_Orbital(const Fl_Geometry& flGeom,
                   const std::string& path = Fl_Tbl_Orbital::getDefaultFileName());
    ~Fl_Tbl_Orbital();

public:
    static std::string getDefaultFileName() {
        return "fl_Table/OrbitalTable";
    }

public:
    int getcGtoTotalNum() const {
        return this->cGtoTotalNum;
    };

    int getcgtonum(const int num) const {
        return this->m_BasisTable[num].cGtonum;
    };

    int getInpAtomnum(const int num) const {
        return this->m_BasisTable[num].InpAtomnumber;
    };

    std::string getInpAtomtype(const int num) const {
        return this->m_BasisTable[num].Atom;
    };

    char getShelltype(int num) const {
        return this->m_BasisTable[num].shell;
    };

    std::string getNote1(int num) const {
        return this->m_BasisTable[num].dumy1;
    };

    std::string getNote2(int num) const {
        return this->m_BasisTable[num].dumy2;
    };

private:
    // ---------- Private Member Function ----------
    void prepare();
    void getMemory();
    void makeTable();      // Table : OrbitalTable .
    void setData();        // SetFunction  to Private Data Member.

private:
    const Fl_Geometry& flGeom_;

    static int  flag1;         // flag for Constructor.
    static int cGtoTotalNum;  // s*1, p*3, d*5.

    int  AtomNum;         // Total Atom Number.
    int   AtomKindNum;     // Number of Atomtype by AtomKind Label2.
    std::string TblFile2;     // Pointer to OutTableFile.
    // This Class make a Table : "OrbitalTable" .

    const std::string m_sTableFilePath;

    static const std::string S_ORB_STR[1];
    static const std::string P_ORB_STR[3];
    static const std::string D_ORB_STR[5];

    struct BasisTable {     // This structure has a line_data in OrbitalTable.
        int   cGtonum;       // Number in Fl_Gto_Orbital.
        int   InpAtomnumber; // Nucleus Number.
        std::string Atom;          // Atomtype.
        char shell;         // Shelltype { s , p , d }.
        std::string dumy1;         // dumy1. : details of shelltype.
        std::string dumy2;         // dumy2. : Label2.
    };

    std::vector<BasisTable> m_BasisTable;

};

#endif // FL_TBL_ORBITAL_H


