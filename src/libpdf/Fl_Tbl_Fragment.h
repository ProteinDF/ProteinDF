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

#ifndef FL_TBL_FRAGMENT_H
#define FL_TBL_FRAGMENT_H

#include <iomanip>
#include <string>
#include "CnError.h"
#include "Fl_Fragment.h"

class Fl_Tbl_Fragment {
 public:
  Fl_Tbl_Fragment(const Fl_Geometry& flGeom);
  ~Fl_Tbl_Fragment();

  int getQclo(int fragindex, int fragqcloindex);
  int getQcloAlpha(int fragindex, int fragqcloindex);
  int getQcloBeta(int fragindex, int fragqcloindex);
  int getFragment(int qcloindex);
  int getFragmentAlpha(int qcloindex);
  int getFragmentBeta(int qcloindex);
  int getFragmentqclo(int qcloindex);
  int getFragmentqcloAlpha(int qcloindex);
  int getFragmentqcloBeta(int qcloindex);
  int getNumberFragmentqclo(int fragindex);

 private:
  static int flag;  // Flag for Constructor

  int number_qclo;           // Number of QCLOs
  int number_fragment;       // Number of Fragments
  int* number_fragmentqclo;  // Number of QCLOs in each fragment

  std::string TblFile;  // Name of the TableFile
  Fl_Fragment FlFrag;   // Pointer to Fl_Fragment object

  struct FragmentTable {     // Line_data in FragmentTable
    int qclo;                // Index of QCLO
    int fragment;            // Index of Fragment
    int fragment_alpha;      // Index of Fragment
    int fragment_beta;       // Index of Fragment
    int fragmentqclo;        // Index of QCLO in each fragment
    int fragmentqclo_alpha;  // Index of QCLO in each fragment
    int fragmentqclo_beta;   // Index of QCLO in each fragment
  };

  FragmentTable* FragTbl;  // Declaration for Structure_Object

  void prepare();  // set QCLO information
  void makeTable();
  void setData();
};

#endif  // FL_TBL_FRAGMENT_H
