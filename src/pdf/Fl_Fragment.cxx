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
#include <fstream>

#include "Fl_Fragment.h"
#include "PdfUtils.h"
#include "TlUtils.h"

#include "CnError.h"
#include "Fl_Geometry.h"
#include "TlPrdctbl.h"


Fl_Fragment::Fl_Fragment(const Fl_Geometry& flGeom) 
    : flGeom_(flGeom), m_filePath("fl_Input/fl_Fragment")
{
}


Fl_Fragment:: ~Fl_Fragment()
{
}


std::string Fl_Fragment::getStr() const
{
    std::string str = "";

    const int numOfFragments = this->m_fragments.size();
    str += TlUtils::format("Frame = [%s]\n", this->m_frameName.c_str());
    str += TlUtils::format("number_of_fragments = %d\n", numOfFragments);

    int counter = 0;
    std::map<int, Fragment>::const_iterator pEnd = this->m_fragments.end();
    for (std::map<int, Fragment>::const_iterator p = this->m_fragments.begin(); p != pEnd; ++p) {
        str += TlUtils::format("// #%d fragment\n", counter);
        str += TlUtils::format("index = %d\n", p->first);
        str += TlUtils::format("name = [%s]\n", p->second.name.c_str());
        str += TlUtils::format("number_of_electrons = %d\n", p->second.numOfElectrons);
        str += TlUtils::format("number_of_alpha_electrons = %d\n", p->second.numOfAlphaElectrons);
        str += TlUtils::format("number_of_beta_electrons = %d\n", p->second.numOfBetaElectrons);
        str += TlUtils::format("number_of_orbitals = %d\n", p->second.numOfOrbitals);
        str += TlUtils::format("number_of_occupied_orbitals = %d\n", p->second.numOfOccupiedOrbitals);
        str += TlUtils::format("number_of_occupied_alpha_orbitals = %d\n", p->second.numOfOccupiedAlphaOrbitals);
        str += TlUtils::format("number_of_occupied_beta_orbitals = %d\n", p->second.numOfOccupiedBetaOrbitals);
        str += TlUtils::format("number_of_unoccupied_orbitals = %d\n", p->second.numOfUnoccupiedOrbitals);
        str += TlUtils::format("number_of_unoccupied_alpha_orbitals = %d\n", p->second.numOfUnoccupiedAlphaOrbitals);
        str += TlUtils::format("number_of_unoccupied_beta_orbitals = %d\n", p->second.numOfUnoccupiedBetaOrbitals);
        str += "\n";

        ++counter;
    }

    return str;
}


void Fl_Fragment::save() const
{
    std::ofstream fs(this->m_filePath.c_str(), std::ios::out | std::ios::trunc);
    fs << this->getStr();
    fs.close();
}


void Fl_Fragment::load()
{
    std::ifstream ifs(this->m_filePath.c_str(), std::ios::in);

    int numOfFragments = 0;
    int currentFragment = 0;
    while (!ifs.eof()) {
        std::string line = "";
        std::getline(ifs, line);

        if (PdfUtils::isComment(line) == true) {
            continue;
        }
        TlUtils::trim_ws(line);

        std::string term1 = TlUtils::getPdfParam(line);
        std::string term2 = TlUtils::getPdfParam(line);
        std::string term3 = TlUtils::getPdfParam(line);

        if ((term1 == "Frame") && (term2 == "=")) {
            this->m_frameName = term3;
            continue;
        } else if ((term1 == "number_of_fragments") && (term2 == "=")) {
            numOfFragments = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "index") && (term2 == "=")) {
            currentFragment = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "name") && (term2 == "=")) {
            this->m_fragments[currentFragment].name = term3;
            continue;
        } else if ((term1 == "number_of_electrons") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfElectrons = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_alpha_electrons") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfAlphaElectrons = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_beta_electrons") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfBetaElectrons = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_occupied_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfOccupiedOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_occupied_alpha_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfOccupiedAlphaOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_occupied_beta_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfOccupiedBetaOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_unoccupied_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfUnoccupiedOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_unoccupied_alpha_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfUnoccupiedAlphaOrbitals = std::atoi(term3.c_str());
            continue;
        } else if ((term1 == "number_of_unoccupied_beta_orbitals") && (term2 == "=")) {
            this->m_fragments[currentFragment].numOfUnoccupiedBetaOrbitals = std::atoi(term3.c_str());
            continue;
        } else {
            std::cerr << TlUtils::format("unknown parameter: %s %s %s\n",
                                         term1.c_str(), term2.c_str(), term3.c_str());
        }
    }

    assert(static_cast<std::size_t>(numOfFragments) == this->m_fragments.size());

    ifs.close();
}


int Fl_Fragment::getNumOfFragments() const
{
    return this->m_fragments.size();
}


int Fl_Fragment::getNumOfOrbitals(const int index) const
{
    int answer = -1;
    std::map<int, Fragment>::const_iterator p = this->m_fragments.find(index);
    if (p != this->m_fragments.end()) {
        answer = p->second.numOfOrbitals;
    }

    return answer;
}


int Fl_Fragment::getNumOfOccupiedOrbitals(const int index) const
{
    int answer = -1;
    std::map<int, Fragment>::const_iterator p = this->m_fragments.find(index);
    if (p != this->m_fragments.end()) {
        answer = p->second.numOfOccupiedOrbitals;
    }

    return answer;
}


int Fl_Fragment::getNumOfOccupiedAlphaOrbitals(const int index) const
{
    int answer = -1;
    std::map<int, Fragment>::const_iterator p = this->m_fragments.find(index);
    if (p != this->m_fragments.end()) {
        answer = p->second.numOfOccupiedAlphaOrbitals;
    }

    return answer;
}


int Fl_Fragment::getNumOfOccupiedBetaOrbitals(const int index) const
{
    int answer = -1;
    std::map<int, Fragment>::const_iterator p = this->m_fragments.find(index);
    if (p != this->m_fragments.end()) {
        answer = p->second.numOfOccupiedBetaOrbitals;
    }

    return answer;
}


void Fl_Fragment::calcDefault(const std::vector<int>& norbcut)
{
    //const Fl_Geometry FlGeom((*this->pPdfParam_)["coordinates"]);
    const TlPrdctbl Prdctbl;

    std::vector<int> electrons(this->m_fragments.size(), 0);
    const int natom = this->flGeom_.getNumOfAtoms();

    // read AtomFragmentTable
    {
        std::vector<int> atomFragment(natom);
        {
            std::ifstream fi;
            int atom, fragment;
            fi.open("fl_Table/AtomFragmentTable");
            while (fi) {
                fi >> atom >> fragment;
                atomFragment[atom] = fragment;
            }
            fi.close();
        }

        // count electrons in each fragment
        for (int i = 0;  i < natom; ++i) {
            const std::string element = this->flGeom_.getAtomSymbol(i);
            const int fragindex = atomFragment[i];
            electrons[fragindex] += Prdctbl.getAtomicNumber(element);
        }
    }

    // set default values
    std::map<int, Fragment>::iterator pEnd = this->m_fragments.end();
    for (std::map<int, Fragment>::iterator p = this->m_fragments.begin(); p != pEnd; ++p) {
        const int index = p->first;
        if (p->second.numOfElectrons == -1) {
            p->second.numOfElectrons = electrons[index];
        }
        if (p->second.numOfAlphaElectrons == -1) {
            p->second.numOfAlphaElectrons = (electrons[index] +1) / 2;
        }
        if (p->second.numOfBetaElectrons == -1) {
            p->second.numOfBetaElectrons = electrons[index] / 2;
        }

        if (p->second.numOfOrbitals == -1) {
            if ((p->second.numOfOccupiedOrbitals != -1) &&
                    (p->second.numOfUnoccupiedOrbitals != -1)) {
                p->second.numOfOrbitals =
                    p->second.numOfOccupiedOrbitals
                    + p->second.numOfUnoccupiedOrbitals;
            } else {
                p->second.numOfOrbitals = norbcut[index];
            }
        }
        if (p->second.numOfOccupiedOrbitals == -1) {
            p->second.numOfOccupiedOrbitals =
                (p->second.numOfElectrons +1) / 2;
        }
        if (p->second.numOfOccupiedAlphaOrbitals == -1) {
            p->second.numOfOccupiedAlphaOrbitals =
                p->second.numOfAlphaElectrons;
        }
        if (p->second.numOfOccupiedBetaOrbitals == -1) {
            p->second.numOfOccupiedBetaOrbitals =
                p->second.numOfBetaElectrons;
        }
        if (p->second.numOfUnoccupiedOrbitals == -1) {
            p->second.numOfUnoccupiedOrbitals =
                norbcut[index] - p->second.numOfOccupiedOrbitals;
        }
        if (p->second.numOfUnoccupiedAlphaOrbitals == -1) {
            p->second.numOfUnoccupiedAlphaOrbitals =
                norbcut[index] - p->second.numOfOccupiedAlphaOrbitals;
        }
        if (p->second.numOfUnoccupiedBetaOrbitals == -1) {
            p->second.numOfUnoccupiedBetaOrbitals =
                norbcut[index] - p->second.numOfOccupiedBetaOrbitals;
        }
    }
}


