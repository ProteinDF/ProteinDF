#ifndef FL_TBL_ATOMFRAGMENT_H
#define FL_TBL_ATOMFRAGMENT_H

#include <iomanip>
#include "CnError.h"

class Fl_Tbl_AtomFragment {
public:
    Fl_Tbl_AtomFragment();
    ~Fl_Tbl_AtomFragment();

    void  makeTable();
    void  setData();
    void  putFragment(int atomindex, int fragindex);
    int  getFragment(int fragindex);

private:
    int number_atom;

    struct AtomFragmentTable {
        int atom;
        int fragment;
    };

    AtomFragmentTable* m_pAtomFragTbl; // Declaration for Structure_Object
};

#endif // FL_TBL_ATOMFRAGMENT_H

