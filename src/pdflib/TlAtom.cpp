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


