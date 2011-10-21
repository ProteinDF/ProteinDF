#ifndef FL_VCT_TABLE_H
#define FL_VCT_TABLE_H

#include <cassert>
#include <vector>
#include <string>

class Fl_Vct_Table {
public:
    struct Info {
public:
        explicit Info(int nAtomBelong =0, int nBasisType =0, double dPreExp =0.0, double dExpAlpha =0.0) :
                atomBelong(nAtomBelong), basisType(nBasisType), preExp(dPreExp), expAlpha(dExpAlpha) {
        }
        Info(const Info& rhs)
                : atomBelong(rhs.atomBelong), basisType(rhs.basisType), preExp(rhs.preExp), expAlpha(rhs.expAlpha) {
        }

        int atomBelong;
        int  basisType;
        double preExp;
        double expAlpha;
    };

public:
    Fl_Vct_Table(const std::string& sFilePath);
    Fl_Vct_Table(const Fl_Vct_Table& rhs);
    virtual ~Fl_Vct_Table();

public:
    inline int getSize() const;

    void resize(int size);
    void addData(const Info& info);

    inline Info  operator[](int index) const;
    inline Info& operator[](int index);

public:
    bool save() const;
    bool load();

protected:
    bool save(std::ofstream& ofs) const;
    bool load(std::ifstream& ifs);

private:
    std::string m_sFilePath;
    std::vector<Info> m_data;
};

////////////////////////////////////////////////////////////////////////
// inline

inline int Fl_Vct_Table::getSize() const
{
    return this->m_data.size();
}

inline Fl_Vct_Table::Info Fl_Vct_Table::operator[](int index) const
{
    assert(0 <= index && index < this->getSize());

    return this->m_data[index];
}

inline Fl_Vct_Table::Info& Fl_Vct_Table::operator[](int index)
{
    assert(0 <= index && index < this->getSize());

    return this->m_data[index];
}

////////////////////////////////////////////////////////////////////////
// derivative

class Fl_Vct_RhoTable : public Fl_Vct_Table {
public:
    Fl_Vct_RhoTable() : Fl_Vct_Table("fl_Work/fl_Vct_Rtable") {
    };
    virtual ~Fl_Vct_RhoTable() {
    }

};

class Fl_Vct_MyuTable : public Fl_Vct_Table {
public:
    Fl_Vct_MyuTable() : Fl_Vct_Table("fl_Work/fl_Vct_Mtable") {
    };
    virtual ~Fl_Vct_MyuTable() {
    }

};

class Fl_Vct_NyuTable : public Fl_Vct_Table {
public:
    Fl_Vct_NyuTable() : Fl_Vct_Table("fl_Work/fl_Vct_Ntable") {
    };
    virtual ~Fl_Vct_NyuTable() {
    }

};

#endif // FL_VCT_TABLE_H
