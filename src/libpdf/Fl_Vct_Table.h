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

#ifndef FL_VCT_TABLE_H
#define FL_VCT_TABLE_H

#include <cassert>
#include <string>
#include <vector>

class Fl_Vct_Table {
 public:
  struct Info {
   public:
    explicit Info(int nAtomBelong = 0, int nBasisType = 0, double dPreExp = 0.0,
                  double dExpAlpha = 0.0)
        : atomBelong(nAtomBelong),
          basisType(nBasisType),
          preExp(dPreExp),
          expAlpha(dExpAlpha) {}
    Info(const Info& rhs)
        : atomBelong(rhs.atomBelong),
          basisType(rhs.basisType),
          preExp(rhs.preExp),
          expAlpha(rhs.expAlpha) {}

    int atomBelong;
    int basisType;
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

  inline Info operator[](int index) const;
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

inline int Fl_Vct_Table::getSize() const { return this->m_data.size(); }

inline Fl_Vct_Table::Info Fl_Vct_Table::operator[](int index) const {
  assert(0 <= index && index < this->getSize());

  return this->m_data[index];
}

inline Fl_Vct_Table::Info& Fl_Vct_Table::operator[](int index) {
  assert(0 <= index && index < this->getSize());

  return this->m_data[index];
}

////////////////////////////////////////////////////////////////////////
// derivative

class Fl_Vct_RhoTable : public Fl_Vct_Table {
 public:
  Fl_Vct_RhoTable() : Fl_Vct_Table("fl_Work/fl_Vct_Rtable"){};
  virtual ~Fl_Vct_RhoTable() {}
};

class Fl_Vct_MyuTable : public Fl_Vct_Table {
 public:
  Fl_Vct_MyuTable() : Fl_Vct_Table("fl_Work/fl_Vct_Mtable"){};
  virtual ~Fl_Vct_MyuTable() {}
};

class Fl_Vct_NyuTable : public Fl_Vct_Table {
 public:
  Fl_Vct_NyuTable() : Fl_Vct_Table("fl_Work/fl_Vct_Ntable"){};
  virtual ~Fl_Vct_NyuTable() {}
};

#endif  // FL_VCT_TABLE_H
