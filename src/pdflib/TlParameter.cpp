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

#include <cassert>
#include <iostream>
#include "TlParameter.h"

////////////////////////////////////////////////////////////////////////
// TlParameter::Param

TlParameter::Param::Param()
{
}

TlParameter::Param::Param(const TlParameter::Param& rhs) : m_data(rhs.m_data)
{
}

TlParameter::Param::~Param()
{
}

TlParameter::Param& TlParameter::Param::operator=(const TlParameter::Param& rhs)
{
    if (this != & rhs) {
        this->m_data = rhs.m_data;
    }

    return *this;
}

std::string TlParameter::Param::operator[](const std::string& sKey) const
{
    std::string sAnswer = "";

    ParamData::const_iterator p = this->m_data.find(sKey);
    if (p != this->m_data.end()) {
        sAnswer = p->second;
    }

    return sAnswer;
}

std::string& TlParameter::Param::operator[](const std::string& sKey)
{
    return this->m_data[sKey];
}

int TlParameter::Param::size() const
{
    return this->m_data.size();
}

std::vector<std::string> TlParameter::Param::getKeywords() const
{
    std::vector<std::string> answer;
    answer.reserve(this->size());

    ParamData::const_iterator p_end = this->m_data.end();
    for (ParamData::const_iterator p = this->m_data.begin();
            p != p_end; ++p) {
        answer.push_back(p->first);
    }

    return answer;
}


void TlParameter::Param::merge(const TlParameter::Param& rhs)
{
    for (ParamData::const_iterator p = rhs.m_data.begin(); p != rhs.m_data.end(); p++) {
        this->m_data[p->first] = p->second;
    }
}

TlParameter::Param& TlParameter::Param::operator+=(const TlParameter::Param& rhs)
{
    this->merge(rhs);

    return *this;
}

////////////////////////////////////////////////////////////////////////
// TlParameter

TlParameter::TlParameter() : m_sFilePath("")
{
}

TlParameter::TlParameter(const TlParameter& rhs) : m_sFilePath(rhs.m_sFilePath), m_data(rhs.m_data)
{
}

TlParameter::~TlParameter()
{
}

int TlParameter::size() const
{
    return this->m_data.size();
}

// TlParameter::const_iterator TlParameter::begin() const {
//   return this->m_data.begin();
// }

// TlParameter::const_iterator TlParameter::end() const {
//   return this->m_data.end();
// }

// TlParameter::iterator TlParameter::begin() {
//   return this->m_data.begin();
// }

// TlParameter::iterator TlParameter::end() {
//   return this->m_data.end();
// }

std::vector<std::string> TlParameter::getGroups() const
{
    std::vector<std::string> answer;
    answer.reserve(this->size());

    Group::const_iterator p_end = this->m_data.end();
    for (Group::const_iterator p = this->m_data.begin();
            p != p_end; ++p) {
        answer.push_back(p->first);
    }

    return answer;
}

void TlParameter::setFilePath(const std::string& sFilePath)
{
    this->m_sFilePath = sFilePath;
}

std::string TlParameter::getFilePath() const
{
    return this->m_sFilePath;
}

TlParameter& TlParameter::operator=(const TlParameter& rhs)
{
    if (this != &rhs) {
        this->m_sFilePath   = rhs.m_sFilePath;
        this->m_data        = rhs.m_data;
    }

    return *this;
}

void TlParameter::merge(const TlParameter& rhs)
{
    for (Group::const_iterator p = rhs.m_data.begin(); p != rhs.m_data.end(); p++) {
        this->m_data[p->first].merge(p->second);
    }
}

TlParameter& TlParameter::operator+=(const TlParameter& rhs)
{
    this->merge(rhs);

    return *this;
}

TlParameter::Param TlParameter::operator[](const std::string& sGroupKey) const
{
    Param param;

    Group::const_iterator p = this->m_data.find(sGroupKey);
    if (p != this->m_data.end()) {
        param = p->second;
    }

    return param;
}

TlParameter::Param& TlParameter::operator[](const std::string& sGroupKey)
{
    return this->m_data[sGroupKey];
}

bool TlParameter::load(const std::string& sFilePath)
{
    this->m_sFilePath = sFilePath;

    std::ifstream ifs;
    ifs.open(sFilePath.c_str(), std::ios::in);

    bool bAnswer = this->load(ifs);

    ifs.close();

    return bAnswer;
}

bool TlParameter::load(std::ifstream& ifs)
{
    bool bAnswer = true;

    if (ifs.is_open() == false) {
    }

    if (bAnswer == true) {
        std::string sLine = "";
        std::string sGroupName = "";

        while (std::getline(ifs, sLine)) {
            while (!sLine.empty()) {
                // 先頭のホワイトスペースを除去
                TlUtils::trim_ws(sLine);

                if ((sLine.compare(0, 2, "//") == 0) || (sLine.compare(0, 2, "--") == 0)) {
                    // comment行
                    sLine = "";
                } else if (sLine.compare(0, 4, ">>>>") == 0) {
                    // group name
                    sLine = sLine.substr(4);
                    sGroupName = TlUtils::getPdfParam(sLine);
                } else {
                    // parameter行
                    const std::string keyword = TlUtils::getPdfParam(sLine);
                    const std::string value   = TlUtils::getPdfParam(sLine);

                    (*this)[sGroupName][keyword] = value;
                }
            }
        }
    }

    return bAnswer;
}

bool TlParameter::save() const
{
    return this->save(this->m_sFilePath);
}

bool TlParameter::save(const std::string& sFilePath) const
{
    assert(!sFilePath.empty());

    std::ofstream ofs;
    ofs.open(sFilePath.c_str(), std::ios::out | std::ios::trunc);

    bool bAnswer = this->save(ofs);

    ofs.close();

    return bAnswer;
}

bool TlParameter::save(std::ofstream& ofs) const
{
    bool bAnswer = true;

    if (ofs.is_open() == false) {
        bAnswer = false;
    }

    if (bAnswer == true) {
        this->print(ofs);
    }

    return bAnswer;
}

