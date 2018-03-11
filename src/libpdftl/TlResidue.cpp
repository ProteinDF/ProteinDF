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

#include "TlResidue.h"
#include <cassert>
#include "TlUtils.h"

const char* TlResidue::residueTable[] = {
    "ALA",   "ALAC-", "ALAN+", "ARG",   "ASNN+", "ASN",   "ASP",
    "CYS",   "GLN",   "GLU",   "GLY",   "GLYC-", "HID",   "HIP",
    "HIS",   "LEU",   "LEUN+", "LEUC-", "LYS",   "MET",   "PHE",
    "PHEC-", "PRO",   "SER",   "SERC-", "THR",   "THRN+", "THRC-",
    "TRP",   "TYR",   "TYRC-", "TYRN+", "VALC-", "VAL",   "HEM"};

int TlResidue::number() {
  return (sizeof(TlResidue::residueTable) / sizeof(char*));
}

// [TH] 仕様がおかしい。開始０に変更した方がよい。
std::string TlResidue::label(int p) {
  assert((1 <= p) && p <= this->number());

  //   if( p>number() ){
  //     CnErr.abort( "TlResidue", "", "label", TlUtils::format("illegal order
  //     of amino acid information table %d", p));
  //   }
  //    printf("%s\n",residue[p-1]);
  return std::string(TlResidue::residueTable[p - 1]);
}

TlResidue::TlResidue() {}

bool TlResidue::isResidue(const std::string& str) {
  bool bAnswer = false;

  int nNumOfResidues = sizeof(TlResidue::residueTable) / sizeof(char*);

  for (int i = 0; i < nNumOfResidues; ++i) {
    if (str == std::string(TlResidue::residueTable[i])) {
      bAnswer = true;
      break;
    }
  }

  return bAnswer;
}
