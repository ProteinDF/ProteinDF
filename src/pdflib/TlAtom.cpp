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

#include <iostream>
#include "TlAtom.h"
#include "TlPosition.h"

#define MAX_NUMBER_OF_ATOMS 110

const char* TlAtom::m_sSymbols[MAX_NUMBER_OF_ATOMS] = {
    "X",
    "H",  "He",
    "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra",
    "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt"
};

const double TlAtom::stdAtomicWeight_[MAX_NUMBER_OF_ATOMS] = {
    0.0,
    1.00794, 4.002602,
    6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762, 32.065, 35.453, 39.948,
    39.0983, 40.078, 44.955912, 47.867, 50.9415, 51.9961, 54.938045, 55.845, 58.933195, 58.6934, 63.546, 65.38, 69.723, 72.64, 74.92160, 78.96, 79.904, 83.798,
    85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.96, 98.0, 101.07, 102.90550, 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,
    132.9054519, 137.327,
    138.90547, 140.116, 140.90765, 144.242, 145.0, 150.36, 151.964, 157.25, 158.92535, 162.500, 164.93032, 167.259, 168.93421, 173.054, 174.9668,
    178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.59, 204.3833, 207.2, 208.98040, 209.0, 210.0, 222.0,
    223.0, 226.0,
    227.0, 232.03806, 231.03588, 238.02891, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0,
    265.0, 268.0, 271.0, 272.0, 270.0, 276.0
};

////////////////////////////////////////////////////////////////////////
// ** van der Waals radii of atom groups **
// hydrogens are ignored;
// according to:
// Linus Pauling "The Nature of the Chemical Bond, 3rd ed.",
//   Cornell University Press, Ithaca, N.Y., 1960.
//
// see p.529 of Gordon M. Barrow "Physical Chemistry, 4th ed.",
//   McGraw-Hill, 1979
// cf. p.26 pf Ooi, Oobatake, Nemethy, Scheraga (1987)
//   Proc Natl Acad Sci USA 84:3086-3090.
const double TlAtom::m_dVdw_radii[] = {
    1.7,  // unknown

    1.95, // Br_sp3_6

    1.8,  // Cl_sp3_6

    1.7,  // C_sp1_6
    1.7,  // C_sp1_7
    1.7,  // C_sp2_6
    1.7,  // C_sp2_7
    2.0,  // C_sp3_6
    2.0,  // C_sp3_7

    1.35, // F_sp3_6

    2.15, // I_sp3_6

    1.5,  // N_sp1_4
    1.5,  // N_sp1_7
    1.5,  // N_sp2_1
    1.5,  // N_sp2_2
    1.5,  // N_sp2_3
    1.5,  // N_sp2_4
    1.5,  // N_sp2_5
    1.5,  // N_sp2_7
    1.5,  // N_sp3_1
    1.5,  // N_sp3_2
    1.5,  // N_sp3_3
    1.5,  // N_sp3_5
    1.5,  // N_sp3_7

    1.4,  // O_sp2_1
    1.4,  // O_sp2_2
    1.4,  // O_sp2_4
    1.4,  // O_sp2_5
    1.4,  // O_sp3_1
    1.4,  // O_sp3_2
    1.4,  // O_sp3_4
    1.4,  // O_sp3_5

    1.9,  // P_sp2_7
    1.9,  // P_sp3_1
    1.9,  // P_sp3_7

    1.85, // S_sp1_7
    1.85, // S_sp2_2
    1.85, // S_sp2_4
    1.85, // S_sp2_6
    1.85, // S_sp2_7
    1.85, // S_sp3_1
    1.85, // S_sp3_2
    1.85, // S_sp3_6
    1.85  // S_sp3_7
};

//
TlAtom::TlAtom(const std::string& symbol, const TlPosition& pos, double charge)
    : charge_(charge), m_dAsa(0), m_position(pos), m_element(0) {
    this->m_properties.reset();
    this->m_sPcClass = TlAtom::UNDEFINED;

    this->setElement(symbol);
}

TlAtom::~TlAtom()
{
}

void TlAtom::setElement(unsigned int n)
{
    this->m_element = n;
}

void TlAtom::setElement(const std::string& symbol)
{
    this->m_element = TlAtom::getElementNumber(symbol);
}

void TlAtom::setCharge(double c)
{
    this->charge_ = c;
}

double TlAtom::getCharge() const
{
    return this->charge_;
}

// TODO: rename atomic number
int TlAtom::getElementNumber(const std::string& symbol)
{
    for (int i = 1; i < MAX_NUMBER_OF_ATOMS; ++i) {
        if (symbol == TlAtom::m_sSymbols[i]) {
            return i;
        }
    }

    return 0;
}


//
int TlAtom::getType() const
{
    bool b_sp1 = this->m_properties.test(sp);
    bool b_sp2 = this->m_properties.test(sp2);
    bool b_sp3 = this->m_properties.test(sp3);

    if (this->m_element == this->getElementNumber("Br")) {
        if (b_sp3 && this->m_sPcClass ==6) {
            return Br_sp3_6;
        }
    } else if (this->m_element == this->getElementNumber("Cl")) {
        if (b_sp3 && this->m_sPcClass == 6) {
            return Cl_sp3_6;
        }
    } else if (this->m_element == this->getElementNumber("C")) {
        if (b_sp1 && this->m_sPcClass == 6) {
            return C_sp1_6;
        } else if (b_sp1 && this->m_sPcClass == 7) {
            return C_sp1_7;
        } else if (b_sp2 && this->m_sPcClass == 6) {
            return C_sp2_6;
        } else if (b_sp2 && this->m_sPcClass == 7) {
            return C_sp2_7;
        } else if (b_sp3 && this->m_sPcClass == 6) {
            return C_sp3_6;
        } else if (b_sp3 && this->m_sPcClass == 7) {
            return C_sp3_7;
        }
    } else if (this->m_element == this->getElementNumber("F")) {
        if (b_sp3 && this->m_sPcClass == 6) {
            return F_sp3_6;
        }
    } else if (this->m_element == this->getElementNumber("I")) {
        if (b_sp3 && this->m_sPcClass == 6) {
            return I_sp3_6;
        }
    } else if (this->m_element == this->getElementNumber("N")) {
        if (b_sp1 && this->m_sPcClass == 4) {
            return N_sp1_4;
        } else if (b_sp1 && this->m_sPcClass == 7) {
            return N_sp1_7;
        } else if (b_sp2 && this->m_sPcClass == 1) {
            return N_sp2_1;
        } else if (b_sp2 && this->m_sPcClass == 2) {
            return N_sp2_2;
        } else if (b_sp2 && this->m_sPcClass == 3) {
            return N_sp2_3;
        } else if (b_sp2 && this->m_sPcClass == 4) {
            return N_sp2_4;
        } else if (b_sp2 && this->m_sPcClass == 5) {
            return N_sp2_5;
        } else if (b_sp2 && this->m_sPcClass == 7) {
            return N_sp2_7;
        } else if (b_sp3 && this->m_sPcClass == 1) {
            return N_sp3_1;
        } else if (b_sp3 && this->m_sPcClass == 2) {
            return N_sp3_2;
        } else if (b_sp3 && this->m_sPcClass == 3) {
            return N_sp3_3;
        } else if (b_sp3 && this->m_sPcClass == 5) {
            return N_sp3_5;
        } else if (b_sp3 && this->m_sPcClass == 7) {
            return N_sp3_7;
        }
    } else if (this->m_element == this->getElementNumber("O")) {
        if (b_sp2 && this->m_sPcClass == 1) {
            return O_sp2_1;
        } else if (b_sp2 && this->m_sPcClass == 2) {
            return O_sp2_2;
        } else if (b_sp2 && this->m_sPcClass == 4) {
            return O_sp2_4;
        } else if (b_sp2 && this->m_sPcClass == 5) {
            return O_sp2_5;
        } else if (b_sp3 && this->m_sPcClass == 1) {
            return O_sp3_1;
        } else if (b_sp3 && this->m_sPcClass == 2) {
            return O_sp3_2;
        } else if (b_sp3 && this->m_sPcClass == 4) {
            return O_sp3_4;
        } else if (b_sp3 && this->m_sPcClass == 5) {
            return O_sp3_5;
        }
    } else if (this->m_element == this->getElementNumber("P")) {
        if (b_sp2 && this->m_sPcClass == 7) {
            return P_sp2_7;
        } else if (b_sp3 && this->m_sPcClass == 1) {
            return P_sp3_1;
        } else if (b_sp3 && this->m_sPcClass == 7) {
            return P_sp3_7;
        }
    } else if (this->m_element == this->getElementNumber("S")) {
        if (b_sp1 && this->m_sPcClass == 7) {
            return S_sp1_7;
        } else if (b_sp2 && this->m_sPcClass == 2) {
            return S_sp2_2;
        } else if (b_sp2 && this->m_sPcClass == 4) {
            return S_sp2_4;
        } else if (b_sp2 && this->m_sPcClass == 6) {
            return S_sp2_6;
        } else if (b_sp2 && this->m_sPcClass == 7) {
            return S_sp2_7;
        } else if (b_sp3 && this->m_sPcClass == 1) {
            return S_sp3_1;
        } else if (b_sp3 && this->m_sPcClass == 2) {
            return S_sp3_2;
        } else if (b_sp3 && this->m_sPcClass == 6) {
            return S_sp3_6;
        } else if (b_sp3 && this->m_sPcClass == 7) {
            return S_sp3_7;
        }
    }
    return UNKNOWN;
}


