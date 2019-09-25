// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef FL_TBL_ATOMFRAGMENT_H
#define FL_TBL_ATOMFRAGMENT_H

#include <iomanip>
#include "CnError.h"

class Fl_Tbl_AtomFragment {
   public:
    Fl_Tbl_AtomFragment();
    ~Fl_Tbl_AtomFragment();

    void makeTable();
    void setData();
    void putFragment(int atomindex, int fragindex);
    int getFragment(int fragindex);

   private:
    int number_atom;

    struct AtomFragmentTable {
        int atom;
        int fragment;
    };

    AtomFragmentTable* m_pAtomFragTbl;  // Declaration for Structure_Object
};

#endif  // FL_TBL_ATOMFRAGMENT_H
