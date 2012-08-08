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

