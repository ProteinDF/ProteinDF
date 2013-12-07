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

#ifndef FL_GEOMETRY_H
#define FL_GEOMETRY_H

#include <vector>
#include "TlAtom.h"
#include "TlSerializeData.h"

// supplies geometry data
//   Bohr radius, 0.529177249 is used in au2an(). In fl_Geometry file,
//   dimention of length is used in a.u.
class Fl_Geometry { //: public File {
public:
    struct AtomData {
    public:
        AtomData() : atom(), label("") {
        }
    public:
        TlAtom atom;
        std::string label;
    };

public:
    //Fl_Geometry();
    Fl_Geometry(const TlSerializeData& data); // for data["coordinates"]
    Fl_Geometry(const Fl_Geometry& rhs);
    ~Fl_Geometry();

public:
    // static std::string getDefaultFileName() {
    //     return "fl_Input/fl_Geometry";
    // }

public:
    void clear();
    void pushBack(const AtomData& atomData);

    /// 原子記号を返す
    std::string getAtom(int i) const;

    /// 電荷を返す
    double getCharge(int i) const;

    /// return distinct label of atom for basis set.
    std::string getLabel(int i) const;

    TlPosition getCoordinate(int i) const;

    // return number of atom in data.
    int getNumOfAtoms() const;

    /// 原子の種類の数を返す(Xを含む)
    int getAtomKindNumber() const;

    /// ダミー原子(X)の数を返す。
    int getNumOfDummyAtoms() const;

    TlSerializeData getSerializeData() const;
    
private:
    void load();
    void save() const;

    void setup(const TlSerializeData& geomData);
    
private:
    // translate value[a.u.] to value[angstrom]
    double au2an() {
        return 0.529177249; // BOHR = 0.529177249
        //return 0.5291772108; // BOHR = 0.529177249
    }

    // translate value[angstrom] to value[a.u.]
    double an2au() {
        return 1.0/au2an(); // inverse of BOHR
    }

private:
    std::string filePath_;
    bool isUpdate_;
    std::vector<AtomData> atoms_;
};

#endif // FL_GEOMETRY_H

