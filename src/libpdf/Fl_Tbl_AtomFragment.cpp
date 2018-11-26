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

#include "Fl_Tbl_AtomFragment.h"
#include <fstream>
#include <ios>
#include "TlLogging.h"

Fl_Tbl_AtomFragment::Fl_Tbl_AtomFragment() {
  this->m_pAtomFragTbl = NULL;
  this->number_atom = 0;
}

Fl_Tbl_AtomFragment::~Fl_Tbl_AtomFragment() {
  if (this->m_pAtomFragTbl != NULL) {
    delete[](this->m_pAtomFragTbl);
  }
}

void Fl_Tbl_AtomFragment::makeTable() {
  TlLogging& log = TlLogging::getInstance();

  std::ofstream fo;
  const std::string sTblFile = "fl_Table/AtomFragmentTable";
  fo.open(sTblFile.c_str(), std::ios::out | std::ios::trunc);
  if (!fo) {
    log.error("Cannot open " + sTblFile);
    CnErr.abort();
  }

  for (int i = 0; i < number_atom; i++) {
    fo << std::setw(8) << this->m_pAtomFragTbl[i].atom << std::setw(6)
       << this->m_pAtomFragTbl[i].fragment << std::endl;
  }
  fo.close();
}

void Fl_Tbl_AtomFragment::setData() {
  TlLogging& log = TlLogging::getInstance();

  std::ifstream fi;
  std::string sTblFile = "fl_Table/AtomFragmentTable";
  fi.open(sTblFile.c_str(), std::ios::in);
  if (!fi) {
    log.error("Cannot open " + sTblFile);
    CnErr.abort("Fl_Tbl_AtomFragment", "", "setData",
                "Cannot open AtomFragmentTable");
  }
  number_atom = 0;

  while (fi) {
    fi >> this->m_pAtomFragTbl[number_atom].atom >>
        this->m_pAtomFragTbl[number_atom].fragment;
    number_atom++;
  }
  fi.close();
}

void Fl_Tbl_AtomFragment::putFragment(int atomindex, int fragindex) {
  AtomFragmentTable* new_AtomFragTbl = new AtomFragmentTable[number_atom + 1];
  for (int i = 0; i < number_atom; i++) {
    new_AtomFragTbl[i] = this->m_pAtomFragTbl[i];
  }
  new_AtomFragTbl[number_atom].atom = atomindex;
  new_AtomFragTbl[number_atom].fragment = fragindex;

  if (this->m_pAtomFragTbl != NULL) {
    delete this->m_pAtomFragTbl;
  }

  this->m_pAtomFragTbl = new_AtomFragTbl;
  number_atom++;
}

int Fl_Tbl_AtomFragment::getFragment(int atomindex) {
  return this->m_pAtomFragTbl[atomindex].fragment;
}
