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

#include <stdio.h>
#include <stdlib.h>
#include "TlPrdctbl.h"
#include "TlUtils.h"

const char* TlPrdctbl::periodicTable[] = {
    "X",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr"
};

TlPrdctbl::TlPrdctbl()
{
}

TlPrdctbl::~TlPrdctbl()
{
}

int TlPrdctbl::getAtomicNumber(const std::string& str)
{
    int answer = -1;
    const int nMaxAtomIndex = sizeof(periodicTable) / sizeof(char*);
    for (int i = 0; i < nMaxAtomIndex; i++) {
        if (str.compare(periodicTable[i]) == 0) {
            answer = i;
            break;
        }
    }

    return answer;
}

std::string TlPrdctbl::getSymbol(int n)
{
    std::string answer = "";
    const int nMaxAtomIndex = sizeof(periodicTable) / sizeof(char*);
    if ((1 <= n) && (n <= nMaxAtomIndex)) {
        answer = std::string(periodicTable[n]);
    }

    return answer;
}

