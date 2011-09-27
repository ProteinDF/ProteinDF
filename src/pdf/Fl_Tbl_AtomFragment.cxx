#include <fstream>
#include <ios>
#include "Fl_Tbl_AtomFragment.h"
#include "TlLogX.h"

Fl_Tbl_AtomFragment::Fl_Tbl_AtomFragment()
{
    this->m_pAtomFragTbl = NULL;
    this->number_atom = 0;
}

Fl_Tbl_AtomFragment :: ~Fl_Tbl_AtomFragment()
{
    if (this->m_pAtomFragTbl != NULL) {
        delete[](this->m_pAtomFragTbl);
    }
}


void Fl_Tbl_AtomFragment::makeTable()
{
    TlLogX& Log = TlLogX::getInstance();

    std::ofstream fo;
    const std::string sTblFile = "fl_Table/AtomFragmentTable";
    fo.open(sTblFile.c_str(), std::ios::out | std::ios::trunc);
    if (!fo) {
        Log << " Cannot open " << sTblFile << "\n";
        CnErr.abort();
    }

    for (int i = 0; i < number_atom; i++) {
        fo << std::setw(8) << this->m_pAtomFragTbl[i].atom
        << std::setw(6) << this->m_pAtomFragTbl[i].fragment
        << std::endl;
    }
    fo.close();
}


void Fl_Tbl_AtomFragment::setData()
{
    TlLogX& Log = TlLogX::getInstance();

    std::ifstream fi;
    std::string sTblFile = "fl_Table/AtomFragmentTable";
    fi.open(sTblFile.c_str(), std::ios::in);
    if (!fi) {
        Log << " Cannot open " << sTblFile << "\n";
        CnErr.abort("Fl_Tbl_AtomFragment", "", "setData", "Cannot open AtomFragmentTable") ;
    }
    number_atom = 0;

    while (fi) {
        fi >> this->m_pAtomFragTbl[ number_atom ].atom
        >> this->m_pAtomFragTbl[ number_atom ].fragment;
        number_atom++;
    }
    fi.close();
}

void Fl_Tbl_AtomFragment::putFragment(int atomindex, int fragindex)
{
    AtomFragmentTable* new_AtomFragTbl = new AtomFragmentTable[number_atom+1];
    for (int i = 0; i < number_atom; i++) {
        new_AtomFragTbl[i] = this->m_pAtomFragTbl[i];
    }
    new_AtomFragTbl[ number_atom ].atom     = atomindex;
    new_AtomFragTbl[ number_atom ].fragment = fragindex;

    if (this->m_pAtomFragTbl != NULL) {
        delete this->m_pAtomFragTbl;
    }

    this->m_pAtomFragTbl = new_AtomFragTbl;
    number_atom++;
}

int Fl_Tbl_AtomFragment::getFragment(int atomindex)
{
    return this->m_pAtomFragTbl[atomindex].fragment;
}

