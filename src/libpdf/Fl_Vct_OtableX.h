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

#ifndef FL_VCT_OTABLEX_H
#define FL_VCT_OTABLEX_H

#include <string>
#include <vector>

class Fl_Vct_OtableX {
   public:
    struct Info {
       public:
        explicit Info(int nAtomBelong = 0, int nBasisType = 0,
                      int nContract = 0, double dPreExp = 0.0,
                      double dExpAlpha = 0.0)
            : atomBelong(nAtomBelong),
              basisType(nBasisType),
              contract(nContract),
              preExp(dPreExp),
              expAlpha(dExpAlpha) {}

        int atomBelong;
        int basisType;
        int contract;
        double preExp;
        double expAlpha;
    };

   public:
    Fl_Vct_OtableX();
    ~Fl_Vct_OtableX();

   public:
    int getSize() const;

    void resize(int size);
    void addData(const Info& info);

    Info operator[](int index) const;
    Info& operator[](int index);

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

#endif  // FL_VCT_OTABLEX_H
