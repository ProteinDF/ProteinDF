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
#include <string>
#include <set>
#include <fstream>

#include "Fl_Geometry.h"
#include "CnError.h"
#include "TlPrdctbl.h"
#include "TlMath.h"

#include "TlUtils.h"
#include "TlFile.h"
#include "PdfUtils.h"
#include "TlPseudoYaml.h"

// Fl_Geometry::Fl_Geometry() : isUpdate_(false)
// {
//     this->clear();
// }


Fl_Geometry::Fl_Geometry(const TlSerializeData& geomData) : isUpdate_(false)
{
    this->setup(geomData);
}


Fl_Geometry::Fl_Geometry(const Fl_Geometry& rhs)
        : filePath_(rhs.filePath_), isUpdate_(false), atoms_(rhs.atoms_)
{
}


Fl_Geometry::~Fl_Geometry()
{
    if ((this->filePath_ != "") &&
        (this->isUpdate_ == true)) {
        this->save();
    }
}

void Fl_Geometry::clear()
{
    this->atoms_.clear();

    this->isUpdate_ = true;
}

void Fl_Geometry::pushBack(const Fl_Geometry::AtomData& atomData)
{
    this->atoms_.push_back(atomData);

    this->isUpdate_ = true;
}


TlSerializeData Fl_Geometry::getSerializeData() const
{
    TlSerializeData data;
    const int size = this->atoms_.size();
    for (int i = 0; i < size; ++i) {
        const std::string symbol = this->atoms_[i].atom.getSymbol();
        const double charge = this->atoms_[i].atom.getCharge();
        const std::string label = this->atoms_[i].label;
        const double x = this->atoms_[i].atom.getPosition().x();
        const double y = this->atoms_[i].atom.getPosition().y();
        const double z = this->atoms_[i].atom.getPosition().z();

        TlSerializeData atom;
        atom["symbol"] = symbol;
        atom["label"] = label;
        TlSerializeData position;
        position.pushBack(x);
        position.pushBack(y);
        position.pushBack(z);
        atom["xyz"] = position;
        atom["charge"] = charge;
        data["coordinates"]["global"].pushBack(atom);
    }

    return data;
}


void Fl_Geometry::save() const
{
    std::ofstream fs(this->filePath_.c_str(), std::ios::out | std::ios::trunc);
    if (fs.is_open() != true) {
        CnErr.abort(TlUtils::format("could not open file \"%s\" at Fl_Geometry::write()", this->filePath_.c_str()));
    }

    fs << " -- fl_Geometry --\n";
    fs << " numAtom       :  " << this->atoms_.size() << "\n";

    fs << " -- Atom Charge Label coordinate\n";
    const int size = this->atoms_.size();
    for (int i = 0; i < size; ++i) {
        const std::string symbol = this->atoms_[i].atom.getSymbol();
        const double charge = this->atoms_[i].atom.getCharge();
        const std::string label = this->atoms_[i].label;
        const double x = this->atoms_[i].atom.getPosition().x();
        const double y = this->atoms_[i].atom.getPosition().y();
        const double z = this->atoms_[i].atom.getPosition().z();

        fs << TlUtils::format("%-2s %+16.10e [%s] %+16.10e %+16.10e %+16.10e\n",
                              symbol.c_str(), charge, label.c_str(),
                              x, y, z);
    }

    fs << " -- end of file --\n";
    fs.close();
}

void Fl_Geometry::load()
{
    //std::cerr << "Fl_Geometry::load() -> " << this->filePath_ << std::endl;
    this->atoms_.clear();

    std::ifstream fs(this->filePath_.c_str(), std::ios::in);
    if (fs.fail()) {
        fs.clear(); // if stream at eof, must clear error.
    }
    fs.seekg(0, std::ios::beg);

    std::string line;
    int atomCount = 0;
    while (std::getline(fs, line)) {
        if (PdfUtils::isComment(line) == true) {
            continue;
        }
        TlUtils::trim_ws(line);

        if (line.substr(0, 7) == "numAtom") {
            std::string term1 = TlUtils::getPdfParam(line);
            std::string term2 = TlUtils::getPdfParam(line);
            std::string term3 = TlUtils::getPdfParam(line);
            int numOfAtoms = std::atoi(term3.c_str());
            this->atoms_.resize(numOfAtoms);
            continue;
        }
        if ((line.substr(0, 13) == "Molecularname") ||
                (line.substr(0, 7) == "numInfo")) {
            // for old format
            continue;
        } else {
            // new format
            std::string term1 = TlUtils::getPdfParam(line);
            std::string term2 = TlUtils::getPdfParam(line);
            std::string term3 = TlUtils::getPdfParam(line);
            std::string term4 = TlUtils::getPdfParam(line);
            std::string term5 = TlUtils::getPdfParam(line);
            std::string term6 = TlUtils::getPdfParam(line);
            std::string term7 = TlUtils::getPdfParam(line);
            std::string term8 = TlUtils::getPdfParam(line);

            std::string atom = "";
            double charge = 0.0;
            std::string label = "";
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;

            if (term8 != "") {
                // obsolete
                atom = term1;
                charge = std::atof(term3.c_str());
                label = term5;
                x = std::atof(term6.c_str());
                y = std::atof(term7.c_str());
                z = std::atof(term8.c_str());
            } else if (term7 != "") {
                // obsolete
                atom = term1;
                charge = std::atof(term2.c_str());
                label = term4;
                x = std::atof(term5.c_str());
                y = std::atof(term6.c_str());
                z = std::atof(term7.c_str());
            } else {
                atom = term1;
                charge = std::atof(term2.c_str());
                label = term3;
                x = std::atof(term4.c_str());
                y = std::atof(term5.c_str());
                z = std::atof(term6.c_str());
            }

            this->atoms_[atomCount].atom.setElement(atom);
            this->atoms_[atomCount].atom.setCharge(charge);
            this->atoms_[atomCount].atom.moveTo(x, y, z);
            this->atoms_[atomCount].label = label;

            ++atomCount;
        }
    }
}


void Fl_Geometry::setup(const TlSerializeData& geomData)
{
    this->atoms_.clear();
    TlSerializeData::MapConstIterator groupEnd = geomData.endMap();
    for (TlSerializeData::MapConstIterator group = geomData.beginMap(); group != groupEnd; ++group) {
        TlSerializeData::ArrayConstIterator atomEnd = group->second.endArray();
        for (TlSerializeData::ArrayConstIterator atom = group->second.beginArray(); atom != atomEnd; ++atom) {

            AtomData ad;
            ad.atom.setElement((*atom)["symbol"].getStr());
            ad.atom.setCharge((*atom)["charge"].getDouble());
            const double x = (*atom)["xyz"].getAt(0).getDouble();
            const double y = (*atom)["xyz"].getAt(1).getDouble();
            const double z = (*atom)["xyz"].getAt(2).getDouble();
            ad.atom.moveTo(x, y, z);
            ad.label = (*atom)["label"].getStr();

//             std::cerr << TlUtils::format("%s (%e, %e, %e)",
//                                          ad.atom.getSymbol().c_str(),
//                                          ad.atom.getPosition().x(),
//                                          ad.atom.getPosition().y(),
//                                          ad.atom.getPosition().z())
//                       << std::endl;
            this->atoms_.push_back(ad);
        }
    }
}


int Fl_Geometry::getAtomKindNumber() const
{
    std::set<std::string> atomKind;
    std::vector<AtomData>::const_iterator pEnd = this->atoms_.end();
    for (std::vector<AtomData>::const_iterator p = this->atoms_.begin(); p != pEnd; ++p) {
        atomKind.insert(p->atom.getSymbol() + p->label);
    }

    return atomKind.size();
}

TlAtom Fl_Geometry::getAtom(int i) const
{
    assert(0 <= i);
    assert(i < static_cast<int>(this->atoms_.size()));

    return this->atoms_[i].atom;
}

std::string Fl_Geometry::getAtomSymbol(int i) const
{
    assert(0 <= i);
    assert(i < static_cast<int>(this->atoms_.size()));

    return this->atoms_[i].atom.getSymbol();
}

double Fl_Geometry::getCharge(int i) const
{
    assert(0 <= i);
    assert(i < static_cast<int>(this->atoms_.size()));

    return this->atoms_[i].atom.getCharge();
}

std::string Fl_Geometry::getLabel(int i) const
{
    assert(0 <= i);
    assert(i < static_cast<int>(this->atoms_.size()));

    return this->atoms_[i].label;
}

TlPosition Fl_Geometry::getCoordinate(int i) const
{
    assert(0 <= i);
    assert(i < static_cast<int>(this->atoms_.size()));

    return this->atoms_[i].atom.getPosition();
}

int Fl_Geometry::getNumOfAtoms() const
{
    return this->atoms_.size();
}

int Fl_Geometry::getNumOfDummyAtoms() const
{
    int count = 0;
    const int numOfAtoms = this->getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        if ("X" == this->getAtomSymbol(i)) {
            ++count;
        }
    }
    return count;
}

