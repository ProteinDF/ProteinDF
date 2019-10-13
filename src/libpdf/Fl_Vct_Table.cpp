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

#include "Fl_Vct_Table.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////////
// Fl_Vct_Table
//

Fl_Vct_Table::Fl_Vct_Table(const std::string& sFilePath)
    : m_sFilePath(sFilePath) {}

Fl_Vct_Table::Fl_Vct_Table(const Fl_Vct_Table& rhs)
    : m_sFilePath(rhs.m_sFilePath), m_data(rhs.m_data) {}

Fl_Vct_Table::~Fl_Vct_Table() {}

void Fl_Vct_Table::resize(int size) {
    assert(0 < size);

    this->m_data.resize(size);
}

void Fl_Vct_Table::addData(const Fl_Vct_Table::Info& info) {
    this->m_data.push_back(info);
}

bool Fl_Vct_Table::save() const {
    bool bAnswer = true;
    std::ofstream ofs;
    ofs.open(this->m_sFilePath.c_str(),
             std::ofstream::out | std::ofstream::binary);

    bAnswer = this->save(ofs);

    ofs.close();

    return bAnswer;
}

bool Fl_Vct_Table::load() {
    bool bAnswer = false;

    std::ifstream ifs;
    ifs.open(this->m_sFilePath.c_str());

    if (ifs.fail()) {
        std::cerr << "[error] could not open file. " << this->m_sFilePath
                  << std::endl;
        std::abort();
    }

    bAnswer = this->load(ifs);

    if (bAnswer != true) {
        std::cerr << "Fl_Vct_Table::load() failed.: " << this->m_sFilePath
                  << std::endl;
    }

    ifs.close();

    return bAnswer;
}

bool Fl_Vct_Table::save(std::ofstream& ofs) const {
    bool bAnswer = true;

    const int nSize = this->getSize();
    ofs.write(reinterpret_cast<const char*>(&nSize), sizeof(int));

    // 下位互換性のため
    //   for (int i = 0; i < nSize; ++i){
    //     const int tmp = this->m_data[i].atomBeint;
    //     ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(int));
    //   }
    //   for (int i = 0; i < nSize; ++i){
    //     const int tmp = this->m_data[i].basisType;
    //     ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(int));
    //   }
    //   for (int i = 0; i < nSize; ++i){
    //     const double tmp = this->m_data[i].preExp;
    //     ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
    //   }
    //   for (int i = 0; i < nSize; ++i){
    //     const double tmp = this->m_data[i].expAlpha;
    //     ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
    //   }

    ofs.write(reinterpret_cast<const char*>(&(this->m_data[0])),
              sizeof(Fl_Vct_Table::Info) * nSize);

    return bAnswer;
}

bool Fl_Vct_Table::load(std::ifstream& ifs) {
    bool bAnswer = false;

    int nSize = 0;
    ifs.read(reinterpret_cast<char*>(&nSize), sizeof(int));

    this->m_data.clear();
    this->resize(nSize);

    //   for (int i = 0; i < nSize; i++){
    //     int tmp;
    //     ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int));
    //     this->m_data[i].atomBeint = tmp;
    //   }
    //   for (int i = 0; i < nSize; i++){
    //     int tmp;
    //     ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int));
    //     this->m_data[i].basisType = tmp;
    //   }
    //   for (int i = 0; i < nSize; i++){
    //     double tmp;
    //     ifs.read(reinterpret_cast<char*>(&tmp), sizeof(double));
    //     this->m_data[i].preExp = tmp;
    //   }
    //   for (int i = 0; i < nSize; i++){
    //     double tmp;
    //     ifs.read(reinterpret_cast<char*>(&tmp), sizeof(double));
    //     this->m_data[i].expAlpha = tmp;
    //   }

    ifs.read(reinterpret_cast<char*>(&(this->m_data[0])),
             sizeof(Fl_Vct_Table::Info) * nSize);

    bAnswer = true;

    return bAnswer;
}
