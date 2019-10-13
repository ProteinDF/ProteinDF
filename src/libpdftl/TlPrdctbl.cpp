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

#include "TlPrdctbl.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "TlUtils.h"

const double TlPrdctbl::UNDEFINED = std::numeric_limits<double>::quiet_NaN();

const char* TlPrdctbl::periodicTable[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};

// https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const double TlPrdctbl::vdwRadii[] = {
    0.00,
    1.20,  // H
    1.40,  // He

    1.82,  // Li
    1.53,  // Be
    1.92,  // B
    1.70,  // C
    1.55,  // N
    1.52,  // O
    1.47,  // F
    1.54,  // Ne

    2.27,       // Na
    1.73,       // Mg
    1.84,       // Al
    2.10,       // Si
    1.80,       // P
    1.80,       // S
    1.75,       // Cl
    1.88,       // Ar
    2.75,       // K
    2.31,       // Ca
    2.11,       // Sc
    UNDEFINED,  // Ti
    UNDEFINED,  // V
    UNDEFINED,  // Cr
    UNDEFINED,  // Mn
    UNDEFINED,  // Fe
    UNDEFINED,  // Co
    1.63,       // Ni
    1.40,       // Cu
    1.39,       // Zn
    1.87,       // Ga
    2.11,       // Ge
    1.85,       // As
    1.90,       // Se
    1.85,       // Br
    2.02,       // Kr
    3.03,       // Rb
    2.49,       // Sr
    UNDEFINED,  // Y
    UNDEFINED,  // Zr
    UNDEFINED,  // Nb
    UNDEFINED,  // Mo
    UNDEFINED,  // Tc
    UNDEFINED,  // Ru
    UNDEFINED,  // Rh
    1.63,       // Pd
    1.72,       // Ag
    1.58,       // Cd
    1.93,       // Id
    2.17,       // Sn
    2.06,       // Sb
    2.06,       // Te
    1.98,       // I
    2.16,       // Xe
    3.43,       // Cs
    2.68,       // Ba
    UNDEFINED,  // La
    UNDEFINED,  // Ce
    UNDEFINED,  // Pr
    UNDEFINED,  // Nd
    UNDEFINED,  // Pm
    UNDEFINED,  // Sm
    UNDEFINED,  // Eu
    UNDEFINED,  // Gd
    UNDEFINED,  // Tb
    UNDEFINED,  // Dy
    UNDEFINED,  // Ho
    UNDEFINED,  // Er
    UNDEFINED,  // Tm
    UNDEFINED,  // Yb
    UNDEFINED,  // Lu
    UNDEFINED,  // Hf
    UNDEFINED,  // Ta
    UNDEFINED,  // W
    UNDEFINED,  // Re
    UNDEFINED,  // Os
    UNDEFINED,  // Ir
    1.75,       // Pt
    1.66,       // Au
    1.55,       // Hg
    1.96,       // Tl
    2.02,       // Pb
    2.07,       // Bi
    1.97,       // Po
    2.02,       // At
    2.20,       // Rn
    3.48,       // Fr
    2.83,       // Ra
    UNDEFINED,  // Ac
    UNDEFINED,  // Th
    UNDEFINED,  // Pa
    1.86,       // U
    UNDEFINED,  // Np
    UNDEFINED,  // Pu
    UNDEFINED,  // Am
    UNDEFINED,  // Cm
    UNDEFINED,  // Bk
    UNDEFINED,  // Cf
    UNDEFINED,  // Es
    UNDEFINED,  // Fm
    UNDEFINED,  // Md
    UNDEFINED,  // No
    UNDEFINED,  // Lr
};

// J. C. Slater, J. Chem. Phys., 41, 3199 (1964).
// DOI: 10.1063/1.1725697
const double TlPrdctbl::BraggSlaterRadii[] = {
    0.00,
    0.25,  // H
    0.50,  // He(undefined)

    1.45,  // Li
    1.05,  // Be
    0.85,  // B
    0.70,  // C
    0.65,  // N
    0.60,  // O
    0.50,  // F
    0.50,  // Ne(undefined)

    1.80,  // Na
    1.50,  // Mg
    1.25,  // Al
    1.10,  // Si
    1.00,  // P
    1.00,  // S
    1.00,  // Cl
    1.00,  // Ar(undefined)

    2.20,  // K
    1.80,  // Ca
    1.60,  // Sc
    1.40,  // Ti
    1.35,  // V
    1.40,  // Cr
    1.40,  // Mn
    1.40,  // Fe
    1.35,  // Co
    1.35,  // Ni
    1.35,  // Cu
    1.35,  // Zn
    1.30,  // Ga
    1.25,  // Ge
    1.15,  // As
    1.15,  // Se
    1.15,  // Br
    1.15,  // Kr(undefined)

    2.35,  // Rb
    2.00,  // Sr
    1.80,  // Y
    1.55,  // Zr
    1.45,  // Nb
    1.45,  // Mo
    1.35,  // Tc
    1.30,  // Ru
    1.35,  // Rh
    1.40,  // Pd
    1.60,  // Ag
    1.55,  // Cd
    1.55,  // In
    1.45,  // Sn
    1.45,  // Sb
    1.40,  // Te
    1.40,  // I
    1.40,  // Xe(undefined)

    2.60,  // Cs
    2.15,  // Ba
    1.95,  // La
    1.85,  // Ce
    1.85,  // Pr
    1.85,  // Nd
    1.85,  // Pm
    1.85,  // Sm
    1.85,  // Eu
    1.80,  // Gd
    1.75,  // Tb
    1.75,  // Dy
    1.75,  // Ho
    1.75,  // Er
    1.75,  // Tu
    1.75,  // Yb
    1.75,  // Lu
    1.55,  // Hf
    1.45,  // Ta
    1.35,  // W
    1.35,  // Re
    1.30,  // Os
    1.35,  // Ir
    1.35,  // Pt
    1.35,  // Au
    1.50,  // Hg
    1.90,  // Ti
    1.60,  // Bi
    1.90,  // Po
    1.90,  // At(undefined)
    1.90,  // Rn(undefined)

    2.15,  // Fr(undefined)
    2.15,  // Ra
    1.95,  // Ac
    1.80,  // Th
    1.80,  // Pa
    1.75,  // U
    1.75,  // Np
    1.75,  // Pu
    1.75   // Am
};

// B. Cordero, et.al., Dalton. Trans., 2831-2838, 2008.
// DOI: 10.1039/b801115j
const double TlPrdctbl::covalentRadii[] = {
    0.00,
    0.31,  // H
    0.28,  // He

    1.28,  // Li
    0.96,  // Be
    0.84,  // B
    0.76,  // C(sp3)
    0.71,  // N
    0.66,  // O
    0.57,  // F
    0.58,  // Ne

    1.66,  // Na
    1.41,  // Mg
    1.21,  // Al
    1.11,  // Si
    1.07,  // P
    1.05,  // S
    1.02,  // Cl
    1.06,  // Ar

    2.03,  // K
    1.76,  // Ca
    1.70,  // Sc
    1.60,  // Ti
    1.53,  // V
    1.39,  // Cr
    1.39,  // Mn(LS)
    1.32,  // Fe(LS)
    1.26,  // Co(LS)
    1.24,  // Ni
    1.32,  // Cu
    1.22,  // Zn
    1.22,  // Ga
    1.20,  // Ge
    1.19,  // As
    1.20,  // Se
    1.20,  // Br
    1.16,  // Kr

    2.20,  // Rb
    1.95,  // Sr
    1.90,  // Y
    1.75,  // Zr
    1.64,  // Nb
    1.54,  // Mo
    1.47,  // Tc
    1.46,  // Ru
    1.42,  // Rh
    1.39,  // Pd
    1.45,  // Ag
    1.44,  // Cd
    1.42,  // In
    1.39,  // Sn
    1.39,  // Sb
    1.38,  // Te
    1.39,  // I
    1.40,  // Xe

    2.44,  // Cs
    2.15,  // Ba
    2.07,  // La
    2.04,  // Ce
    2.03,  // Pr
    2.01,  // Nd
    1.99,  // Pm
    1.98,  // Sm
    1.98,  // Eu
    1.96,  // Gd
    1.94,  // Tb
    1.92,  // Dy
    1.92,  // Ho
    1.89,  // Er
    1.90,  // Tm
    1.87,  // Yb
    1.87,  // Lu
    1.75,  // Hf
    1.70,  // Ta
    1.62,  // W
    1.51,  // Re
    1.44,  // Os
    1.41,  // Ir
    1.36,  // Pt
    1.36,  // Au
    1.32,  // Hg
    1.45,  // Tl
    1.46,  // Pb
    1.48,  // Bi
    1.40,  // Po
    1.50,  // At
    1.50,  // Rn
    2.60,  // Fr
    2.21,  // Ra
    2.15,  // Ac
    2.06,  // Th
    2.00,  // Pa
    1.96,  // U
    1.90,  // Np
    1.87,  // Pu
    1.80,  // Am
    1.69,  // Cm
};

TlPrdctbl::TlPrdctbl() {}

TlPrdctbl::~TlPrdctbl() {}

int TlPrdctbl::getAtomicNumber(const std::string& str) {
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

std::string TlPrdctbl::getSymbol(int n) {
    std::string answer = "";
    const int maxAtomIndex = sizeof(periodicTable) / sizeof(char*);
    if ((1 <= n) && (n < maxAtomIndex)) {
        answer = std::string(periodicTable[n]);
    }

    return answer;
}

double TlPrdctbl::getVdwRadii(const int atomicNumber) {
    double answer = UNDEFINED;

    const int maxAtomicNumber = sizeof(TlPrdctbl::vdwRadii) / sizeof(double);
    if ((1 <= atomicNumber) && (atomicNumber < maxAtomicNumber)) {
        answer = TlPrdctbl::vdwRadii[atomicNumber];
    }

    return answer;
}

double TlPrdctbl::getBraggSlaterRadii(const int atomicNumber) {
    double answer = 0.0;

    const int maxAtomicNumber =
        sizeof(TlPrdctbl::BraggSlaterRadii) / sizeof(double);
    if ((1 <= atomicNumber) && (atomicNumber < maxAtomicNumber)) {
        answer = TlPrdctbl::BraggSlaterRadii[atomicNumber];
    }

    return answer;
}

double TlPrdctbl::getCovalentRadii(const int atomicNumber) {
    double answer = 0.0;

    const int maxAtomicNumber =
        sizeof(TlPrdctbl::covalentRadii) / sizeof(double);
    if ((1 <= atomicNumber) && (atomicNumber < maxAtomicNumber)) {
        answer = TlPrdctbl::covalentRadii[atomicNumber];
    }

    return answer;
}
