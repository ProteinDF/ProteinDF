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

