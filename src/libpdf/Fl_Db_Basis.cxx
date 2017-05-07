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
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>

#include "Fl_Db_Basis.h"
#include "TlUtils.h"
#include "TlStringTokenizer.h"

Fl_Db_Basis::Fl_Db_Basis(const std::string& sBasisName) : m_sBasisName(sBasisName)
{
    this->setData();
}


Fl_Db_Basis::~Fl_Db_Basis()
{
}


void Fl_Db_Basis::setData()
{
    std::string DbFile = "basis2";
    const char* PdfHome = std::getenv("PDF_HOME");
    if (PdfHome != NULL) {
        DbFile = TlUtils::format("%s/data/basis2", PdfHome);
    }

    std::ifstream fi;
    fi.open(DbFile.c_str(), std::ios::in);
    if (!fi) {
        CnErr.abort(TlUtils::format("Cannot open %s", DbFile.c_str()));
    }

    const bool auxMode = (this->m_sBasisName.compare(0, 2, "A-") == 0);

    bool isNormalTerminate = false;
    bool isNameFound = false;
    bool isRead_J_part = false;
    static const int MAX_SHELL_TYPE = 5; // s, p, d, f, g
    std::vector<int> numOfCGTOs(MAX_SHELL_TYPE, 0); 
    while (fi) {
        if (isNormalTerminate == true) {
            break;
        }

        std::string line = "";
        std::getline(fi, line);
        
        TlUtils::trim_ws(line);
        TlUtils::rtrim_ws(line);
        if ((line.size() > 0) && (line[0] == '#')) {
            line = "";
        }
        if (line.empty() == true) {
            continue;
        }

        if (isNameFound != true) {
            if (this->m_sBasisName.compare(line) == 0) {
                isNameFound = true;
                continue;
            }
        } else {
            std::vector<std::string> parts = this->getWords(line);
            const int numOfParts = std::min<int>(parts.size(), MAX_SHELL_TYPE);
            for (int i = 0; i < numOfParts; ++i) {
                numOfCGTOs[i] = std::atoi(parts[i].c_str());
            }
            
            if (auxMode != true) {
                this->read_cGTOs(numOfCGTOs[0], "s", fi, &(this->cgto));
                this->read_cGTOs(numOfCGTOs[1], "p", fi, &(this->cgto));
                this->read_cGTOs(numOfCGTOs[2], "d", fi, &(this->cgto));
                this->read_cGTOs(numOfCGTOs[3], "f", fi, &(this->cgto));
                this->read_cGTOs(numOfCGTOs[4], "g", fi, &(this->cgto));
                
                isNormalTerminate = true;
            } else {
                if (isRead_J_part != true) {
                    // J
                    this->read_cGTOs(numOfCGTOs[0], "s", fi, &(this->rhoCgto));
                    this->read_cGTOs(numOfCGTOs[1], "p", fi, &(this->rhoCgto));
                    this->read_cGTOs(numOfCGTOs[2], "d", fi, &(this->rhoCgto));
                    this->read_cGTOs(numOfCGTOs[3], "f", fi, &(this->rhoCgto));
                    this->read_cGTOs(numOfCGTOs[4], "g", fi, &(this->rhoCgto));

                    isRead_J_part = true;
                } else {
                    // K
                    std::vector<std::string> parts = this->getWords(line);
                    const int numOfParts = std::min<int>(parts.size(), 4);
                    for (int i = 0; i < numOfParts; ++i) {
                        numOfCGTOs[i] = std::atoi(parts[i].c_str());
                    }
                    
                    this->read_cGTOs(numOfCGTOs[0], "s", fi, &(this->myuCgto));
                    this->read_cGTOs(numOfCGTOs[1], "p", fi, &(this->myuCgto));
                    this->read_cGTOs(numOfCGTOs[2], "d", fi, &(this->myuCgto));
                    this->read_cGTOs(numOfCGTOs[3], "f", fi, &(this->myuCgto));
                    this->read_cGTOs(numOfCGTOs[4], "g", fi, &(this->myuCgto));
                    
                    isNormalTerminate = true;
                }
            }
        }
    }

    fi.close();

    if (isNormalTerminate == false) {
        CnErr.abort(TlUtils::format("The orbital basis set \'%s\' is not found. Please check Basis set name.",
                                    this->m_sBasisName.c_str()));
    }
}


void Fl_Db_Basis::read_cGTOs(const int numOfCGTOs,
                             const std::string& shellType,
                             std::ifstream& fi,
                             std::vector<Cgto>* pCGTOs)
{
    std::string line = "";

    for (int i = 0; i < numOfCGTOs; ++i) {
        Cgto tmp_cgto;
        tmp_cgto.scalefactor = 1.0;
        tmp_cgto.normalizedfactor = 1.0;

        int GtoNum = 0;
        fi >> GtoNum;
        tmp_cgto.contraction.resize(GtoNum);
        for (int j = 0; j < GtoNum; ++j) {
            do {
                std::getline(fi, line);
                TlUtils::trim_ws(line);
                TlUtils::rtrim_ws(line);
                if ((line.size() > 0) && (line[0] == '#')) {
                    line = "";
                }
            } while (line.empty());

            std::vector<std::string> parts = this->getWords(line);
            tmp_cgto.contraction[j].exponent = std::atof(parts[0].c_str());
            if (parts.size() >= 2) {
                tmp_cgto.contraction[j].coefficient = std::atof(parts[1].c_str());
            }
        }

        tmp_cgto.shell = shellType[0];
        pCGTOs->push_back(tmp_cgto);
    }
}


std::vector<std::string> Fl_Db_Basis::getWords(const std::string& str)
{
    std::vector<std::string> answer;
    
    TlStringTokenizer st(str);
    while(st.hasMoreTokens()) {
        const std::string word = st.nextToken();
        answer.push_back(word);
    }

    return answer;
}
