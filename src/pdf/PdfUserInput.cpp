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

#include <fstream>
#include <iostream>
#include "PdfUserInput.h"
#include "Fl_Db_Basis.h"
#include "TlStringTokenizer.h"
#include "TlAtom.h"
#include "TlResidue.h"
#include "TlLogging.h"

#define BOHR_ANGSTROM (0.52917721092)

PdfUserInput::PdfUserInput(const std::string& filePath)
    : filePath_(filePath), data_(), param_(), log_(TlLogging::getInstance())
{
}

PdfUserInput::~PdfUserInput()
{
}

TlParameter PdfUserInput::getParameter() const
{
    return this->param_;
}

TlSerializeData PdfUserInput::getSerializeData() const
{
    return this->data_;
}


void PdfUserInput::load()
{
    this->load_conventional(); // fl_Userinput形式
}

// メモ
//
// o コメント行は"//"もしくは"--", "#"で開始します。
// o グループ行は">>>>"で始まり、空白もしくは改行までです。
void PdfUserInput::load_conventional()
{
    std::ifstream ifs;
    ifs.open(this->filePath_.c_str(), std::ios::in);

    if (!ifs.is_open()) {
        this->log_.critical(TlUtils::format("could not open file: %s. stop.",
                                            this->filePath_.c_str()));
        abort();
    }

    std::string sLine = "";
    std::string sGroup = "";
    std::string sKeyword = "";
    bool bNextValue = false; // keywordを読み込み後、"="を読み込んだらtrue。それ以外false
    int nNumOfLine = 0;

    while (!ifs.eof()) {
        // 1行読み込む
        std::string sTmp = "";
        std::getline(ifs, sTmp);
        sLine += sTmp;
        nNumOfLine++;

        // parse
        while (!sLine.empty()) {
            //std::cerr << "sline = " << sLine << std::endl;

            // 行頭のホワイトスペースを除去
            TlUtils::trim_ws(sLine);
            if (sLine.empty()) {
                break;
            }

            // コメント行は改行まで読み込みを破棄する
            if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 2, "--") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
                sLine = "";
                break;
            }

            if (sKeyword == "") {
                // keyword が設定されていない場合 ==============================

                // グループ行
                if (sLine.compare(0, 4, ">>>>") == 0) {
                    sLine = sLine.substr(4);
                    sGroup = TlUtils::getWord(sLine);
                    continue;
                }

                // キーワード行
                // "=" とホワイトスペースの前まで取得
                std::string::size_type nKeyword = sLine.find_first_of("= \f\n\r\t\v");
                sKeyword = sLine.substr(0, nKeyword);
                if (nKeyword == std::string::npos) {
                    sLine = "";
                } else {
                    sLine = sLine.substr(nKeyword);
                }
            } else {
                // keyword が設定されている場合 ================================
                if (bNextValue == false) {
                    // "=" がまだ読み込まれていない。
                    // "=" を読み込む。"="で無ければエラー
                    if (sLine[0] == '=') {
                        sLine = sLine.substr(1);
                        bNextValue = true;
                        continue;
                    } else {
                        // error
                        std::cerr << "could not found '='. line =" << nNumOfLine << ". stop." << std::endl;
                        //std::cerr << "sKeyword = " << sKeyword << std::endl;
                        //std::cerr << "sLine = " << sLine << std::endl;
                        abort();
                    }
                } else {
                    // "=" が読み込まれている。

                    // value を読み込む
                    // 括弧の種類によって分類
                    const char bracketType = sLine[0];

                    std::string sValue = "";
                    switch (bracketType) {
                    case '[':
                        {
                            // "[" で始まるvalue は"]" まで格納する(改行を含まない)。
                            const std::string::size_type nEndBracket = sLine.find_first_of(']');
                            if (nEndBracket != std::string::npos) {
                                sValue = sLine.substr(1, nEndBracket -1);
                                //this->m_param[sGroup][sKeyword] = sValue;
                                sLine = sLine.substr(nEndBracket +1);
                                
                                // 次のキーワード、グループを読み込むための後始末
                                //sKeyword = "";
                                //bNextValue = false;
                            } else {
                                // error
                                std::cerr << "could not found ']'. line =" << nNumOfLine << ". stop." << std::endl;
                                abort();
                            }
                        }
                        break;

                    case '{':
                        {
                            const int nStartBracketLine = nNumOfLine;
                            // "{" で始まるvalue は"}"もしくは"}end" まで格納する(改行を含む)。
                            bool bFoundEndBracket = false;
                            do {
                                std::string::size_type nEndBracket = sLine.find_first_of("}");
                                if (nEndBracket != std::string::npos) {
                                    // 終端記号を発見した
                                    sValue = sLine.substr(1, nEndBracket -1);
                                    //this->m_param[sGroup][sKeyword] = sValue;
                                    sLine = sLine.substr(nEndBracket +1);
                                    if (sLine.compare(0, 3, "end") == 0) {
                                        // "}end"の場合も終端記号としてOKにする<-互換のため
                                        sLine = sLine.substr(3);
                                    }
                                    
                                    // 次のキーワード、グループを読み込むための後始末
                                    //sKeyword = "";
                                    bFoundEndBracket = true;
                                    
                                    break; // while を抜ける
                                }

                                // 終端記号が見つからなかったので、次の行を読み込む
                                std::string sTmp = "";
                                std::getline(ifs, sTmp);
                                sLine += ("\n" + sTmp);
                                nNumOfLine++;
                            } while (!ifs.eof());
                            
                            if (bFoundEndBracket == false) {
                                // error
                                // 終端記号が見つからないままファイルの終わりまで来てしまった
                                std::cerr << "could not found '}'. line =" << nStartBracketLine << ". stop." << std::endl;
                                abort();
                            }
                        }
                        break;

                    default:
                        {
                            // 括弧で囲まれていないとき
                            // ホワイトスペースまでが値
                            std::string::size_type nValue = sLine.find_first_of(" \f\n\r\t\v");
                            //std::string sValue;
                            if (nValue != std::string::npos) {
                                sValue = sLine.substr(0, nValue);
                                sLine = sLine.substr(nValue);
                            } else {
                                sValue = sLine;
                                sLine = "";
                            }
                            //this->m_param[sGroup][sKeyword] = sValue;
                            //sKeyword = "";
                            //bNextValue = false;
                        }
                        break;
                    }
                    
                    // store value
                    //std::cerr << "g = " << sGroup << ", k = " << sKeyword << ", v = " << sValue << std::endl;
                    this->param_[sGroup][sKeyword] = sValue;
                    this->data_[sKeyword] = sValue;

                    sKeyword = "";
                    bNextValue = false;
                }
            }
        }
    }

    ifs.close();

    // make table ========================================================
    if (this->data_["geometry/cartesian/input"].getStr() != "") {
        const std::string str = this->data_["geometry/cartesian/input"].getStr();
        this->molecule_geometry_cartesian_input(str);

        this->data_["geometry/cartesian/input"] = "";
        this->param_["MOLECULE"]["geometry/cartesian/input"] = "stored";
    }

    if (this->data_["basis-set/orbital"].getStr() != "") {
        const std::string str = this->data_["basis-set/orbital"].getStr();
        this->moleculeBasisSetOrbital(str);

        this->data_["basis-set/orbital"] = "";
        this->param_["MOLECULE"]["basis-set/orbital"] = "stored";
    }

    if (this->data_["basis-set/density-auxiliary"].getStr() != "") {
        const std::string str = this->data_["basis-set/density-auxiliary"].getStr();
        this->moleculeBasisSetDensityAuxiliary(str);
        
        this->data_["basis-set/density-auxiliary"] = "s";
        this->param_["MOLECULE"]["basis-set/density-auxiliary"] = "stored";
    }

    if (this->data_["basis-set/exchange-auxiliary"].getStr() != "") {
        const std::string str = this->data_["basis-set/exchange-auxiliary"].getStr();
        this->moleculeBasisSetExchangeAuxiliary(str);
        
        this->data_["basis-set/exchange-auxiliary"] = "";
        this->param_["MOLECULE"]["basis-set/exchange-auxiliary"] = "stored";
    }

    if (this->data_["basis-set/gridfree"].getStr() != "") {
        const std::string str = this->data_["basis-set/gridfree"].getStr();
        this->moleculeBasisSetGridFree(str);
        
        this->data_["basis-set/gridfree"] = "";
        this->param_["MOLECULE"]["basis-set/gridfree"] = "stored";
    } else {
        this->log_.info(" use orbital basis-set for gridfree.");
        this->data_["basis_sets_GF"] = this->data_["basis_sets"];

        this->data_["basis-set/gridfree"] = "";
        this->param_["MOLECULE"]["basis-set/gridfree"] = "stored";
    }

    //this->m_param.print(std::cout);

    this->alias();
    this->check();
}


void PdfUserInput::molecule_geometry_cartesian_input(const std::string& str)
{
    std::istringstream in(str);

    while (in) {
        std::string sLine;
        std::getline(in, sLine);

        // 先頭のホワイトスペースを取り除く
        TlUtils::trim_ws(sLine);
        if (sLine.empty()) {
            continue;
        }

        // コメント行を除く
        if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
            continue;
        }

        // ホワイトスペースをデリミタにして、読み込み
        std::string sAtom    = "";
        std::string sLabel1  = "";
        std::string sLabel2  = "";
        TlPosition position;
        double charge = 0.0;

        TlStringTokenizer st(sLine);
        std::string tmp = st.nextToken();

        // Atom
        if ((tmp == "X") || (TlAtom::getElementNumber(tmp) != 0)) {
            sAtom = tmp;
            tmp = st.nextToken();
        } else {
            this->log_.critical(TlUtils::format("atom \"%s\" is not defined. stop.", tmp.c_str()));
            abort();
        }

        // Label2
        if (tmp[0] == '@') {
            sLabel2 = tmp.substr(1);
            tmp = st.nextToken();
        }

        // 座標
        {
            double x = std::atof(tmp.c_str());
            tmp = st.nextToken();
            double y = std::atof(tmp.c_str());
            tmp = st.nextToken();
            double z = std::atof(tmp.c_str());
            tmp = st.nextToken();

            position.moveTo(x, y, z);

            // angstrom -> a.u.
            if (this->data_["geometry/cartesian/unit"].getStr() == "angstrom") {
                //position /= 0.529177249;
                position /= BOHR_ANGSTROM;
            }
        }

        // 電荷
        if ((tmp.compare(0, 2, "//") == 0) ||
                (tmp.compare(0, 1, "#") == 0)) {
            tmp = "";
        }
        if (!tmp.empty()) {
            charge = std::atof(tmp.c_str());
        } else {
            charge = TlAtom::getElementNumber(sAtom.c_str());
        }

        // serializeDataに格納
        TlSerializeData atom;
        atom["symbol"] = sAtom;
        atom["charge"] = charge;
        TlSerializeData xyz;
        xyz.pushBack(position.x());
        xyz.pushBack(position.y());
        xyz.pushBack(position.z());
        atom["xyz"] = xyz;
        atom["label"] = sLabel2;
        this->data_["coordinates"]["_"].pushBack(atom);
    }
}

void PdfUserInput::moleculeBasisSetOrbital(const std::string& str)
{
    std::istringstream in(str);

    int dNumOfCgto = 0;

    while (in) {
        std::string sLine;
        std::getline(in, sLine);

        // 先頭のホワイトスペースを取り除く
        TlUtils::trim_ws(sLine);
        if (sLine.empty()) {
            continue;
        }

        // コメント行を除く
        if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
            continue;
        }

        // ホワイトスペースをデリミタにして、読み込み
        std::string sLabel = "";
        std::string sAlias = "";

        // Atom
        std::string sAtom = "";
        {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            if ((tmp == "X") || (TlAtom::getElementNumber(tmp) != 0)) {
                sAtom = tmp;
            } else {
                this->log_.critical(TlUtils::format("atom \"%s\" is not defined. stop.", tmp.c_str()));
                abort();
            }
        }

        // label2
        if (sLine[0] == '@') {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            tmp = tmp.substr(1);
            sLabel = tmp;
        }

        // equal
        if (sLine[0] == '=') {
            sLine = sLine.substr(1);
            TlUtils::trim_ws(sLine);
        } else {
            std::cerr << "equal is not found. stop." << std::endl;
            abort();
        }

        // Name
        std::string sName = "";
        if (sLine[0] == '"') {
            sLine = sLine.substr(1);
            std::string::size_type nNameEnd = sLine.find_first_of('"');
            if (nNameEnd != std::string::npos) {
                sName = sLine.substr(0, nNameEnd);
                sLine = (nNameEnd +1 < sLine.length()) ? sLine.substr(nNameEnd +1) : "";
            } else {
                std::cerr << "double quotation is not closed. stop." << std::endl;
                abort();
            }
        } else {
            sName = TlUtils::getWord(sLine);
        }

        // store to serializeData
        std::string key = sAtom;
        if (sLabel != "") {
            key += "@" + sLabel;
        }
        const TlSerializeData basisset = this->getBasisInfo(sName);
        this->data_["basis_sets"][key] = basisset;
    }
}


TlSerializeData PdfUserInput::getBasisInfo(const std::string& basisName)
{
    TlSerializeData data;
    data["name"] = basisName;
    
    Fl_Db_Basis flDbBasis(basisName);
    const int numOfCGTOs = flDbBasis.getTotalcgto();
    for (int i = 0; i < numOfCGTOs; ++i) {
        TlSerializeData cGTO;
        cGTO["shell_type"] = TlUtils::format("%c", flDbBasis.getShell(i));
        cGTO["scale_factor"] = flDbBasis.getScalefactor(i);
        const int contractions = flDbBasis.getContraction(i);

        for (int j = 0; j < contractions; ++j) {
            TlSerializeData pGTO;
            pGTO["exp"] = flDbBasis.getExpornent(i, j);
            pGTO["coef"] = flDbBasis.getCoefficient(i ,j);
            cGTO["pGTOs"].pushBack(pGTO);
        }

        data["cGTOs"].pushBack(cGTO);
    }

    return data;
}


void PdfUserInput::moleculeBasisSetDensityAuxiliary(const std::string& str)
{
    std::istringstream in(str);

    int dNumOfCgto = 0;

    while (in) {
        std::string sLine;
        std::getline(in, sLine);

        // 先頭のホワイトスペースを取り除く
        TlUtils::trim_ws(sLine);
        if (sLine.empty()) {
            continue;
        }

        // コメント行を除く
        if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
            continue;
        }

        // ホワイトスペースをデリミタにして、読み込み
        std::string sLabel = "";
        //std::string sResidue = "";
        std::string sAlias   = "";

        // Atom
        std::string sAtom = "";
        {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            if ((tmp == "X") || (TlAtom::getElementNumber(tmp) != 0)) {
                sAtom = tmp;
            } else {
                this->log_.critical(TlUtils::format("atom \"%s\" is not defined. stop.", tmp.c_str()));
                abort();
            }
        }

        // label2
        if (sLine[0] == '@') {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            tmp = tmp.substr(1);
            sLabel = tmp;
        }

        // equal
        if (sLine[0] == '=') {
            sLine = sLine.substr(1);
            TlUtils::trim_ws(sLine);
        } else {
            std::cerr << "equal is not found. stop." << std::endl;
            abort();
        }

        // Name
        std::string sName = "";
        {
            if (sLine[0] == '"') {
                sLine = sLine.substr(1);
                std::string::size_type nNameEnd = sLine.find_first_of('"');
                if (nNameEnd != std::string::npos) {
                    sName = sLine.substr(0, nNameEnd);
                    sLine = (nNameEnd +1 < sLine.length()) ? sLine.substr(nNameEnd +1) : "";
                } else {
                    std::cerr << "double quotation is not closed. stop." << std::endl;
                    abort();
                }
            } else {
                sName = TlUtils::getWord(sLine);
            }
        }

        // 格納 ============================================================
        Fl_Db_Basis flDbBasis(sName);

        // store to serializeData
        std::string key = sAtom;
        if (sLabel != "") {
            key += "@" + sLabel;
        }
        const int numOfCGTOs = flDbBasis.getrouTotalnum();
        for (int i = 0; i < numOfCGTOs; ++i) {
            TlSerializeData cGTO;
            cGTO["shell_type"] = TlUtils::format("%c", flDbBasis.getrouShell(i));
            cGTO["scale_factor"] = 1.0;

            const int contraction = flDbBasis.getrouContraction(i);
            for (int j = 0; j < contraction; j++) {
                TlSerializeData pGTO;
                pGTO["exp"] = flDbBasis.getrouExpornent(i, j);
                pGTO["coef"] = 1.0;
                cGTO["pGTOs"].pushBack(pGTO);
            }

            this->data_["basis_sets_j"][key]["cGTOs"].pushBack(cGTO);
        }
        this->data_["basis_sets_j"][key]["name"] = sName;
    }
}

void PdfUserInput::moleculeBasisSetExchangeAuxiliary(const std::string& str)
{
    std::istringstream in(str);

    int dNumOfCgto = 0;

    while (in) {
        std::string sLine;
        std::getline(in, sLine);

        // 先頭のホワイトスペースを取り除く
        TlUtils::trim_ws(sLine);
        if (sLine.empty()) {
            continue;
        }

        // コメント行を除く
        if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
            continue;
        }

        // ホワイトスペースをデリミタにして、読み込み
        std::string sLabel = "";
        //std::string sResidue = "";
        std::string sAlias   = "";

        // Atom
        std::string sAtom = "";
        {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            if ((tmp == "X") || (TlAtom::getElementNumber(tmp) != 0)) {
                sAtom = tmp;
            } else {
                this->log_.critical(TlUtils::format("atom \"%s\" is not defined. stop.", tmp.c_str()));
                abort();
            }
        }

        // label2
        if (sLine[0] == '@') {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            tmp = tmp.substr(1);
            sLabel = tmp;
        }

        // equal
        if (sLine[0] == '=') {
            sLine = sLine.substr(1);
            TlUtils::trim_ws(sLine);
        } else {
            std::cerr << "equal is not found. stop." << std::endl;
            abort();
        }

        // Name
        std::string sName = "";
        {
            if (sLine[0] == '"') {
                sLine = sLine.substr(1);
                std::string::size_type nNameEnd = sLine.find_first_of('"');
                if (nNameEnd != std::string::npos) {
                    sName = sLine.substr(0, nNameEnd);
                    sLine = (nNameEnd +1 < sLine.length()) ? sLine.substr(nNameEnd +1) : "";
                } else {
                    std::cerr << "double quotation is not closed. stop." << std::endl;
                    abort();
                }
            } else {
                sName = TlUtils::getWord(sLine);
            }
        }

        // 格納 ============================================================
        Fl_Db_Basis flDbBasis(sName);

        // store to serializeData
        std::string key = sAtom;
        if (sLabel != "") {
            key += "@" + sLabel;
        }
        const int numOfCGTOs = flDbBasis.getmyuTotalnum();
        for (int i = 0; i < numOfCGTOs; ++i) {
            TlSerializeData cGTO;
            cGTO["shell_type"] = TlUtils::format("%c", flDbBasis.getmyuShell(i));
            cGTO["scale_factor"] = 1.0;
            const int contractions = flDbBasis.getmyuContraction(i);

            for (int j = 0; j < contractions; ++j) {
                TlSerializeData pGTO;
                pGTO["exp"] = flDbBasis.getmyuExpornent(i, j);
                pGTO["coef"] = 1.0;
                cGTO["pGTOs"].pushBack(pGTO);
            }

            this->data_["basis_sets_k"][key]["cGTOs"].pushBack(cGTO);
        }
        this->data_["basis_sets_k"][key]["name"] = sName;
    }
}

void PdfUserInput::moleculeBasisSetGridFree(const std::string& str)
{
    // Fl_Gto_Orbital Bas;
    std::istringstream in(str);
    int dNumOfCgto = 0;

    while (in) {
        std::string sLine;
        std::getline(in, sLine);

        // 先頭のホワイトスペースを取り除く
        TlUtils::trim_ws(sLine);
        if (sLine.empty()) {
            continue;
        }

        // コメント行を除く
        if ((sLine.compare(0, 2, "//") == 0) ||
                (sLine.compare(0, 1, "#") == 0)) {
            continue;
        }

        // ホワイトスペースをデリミタにして、読み込み
        std::string sLabel = "";
        std::string sAlias = "";

        // Atom
        std::string sAtom = "";
        {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            if ((tmp == "X") || (TlAtom::getElementNumber(tmp) != 0)) {
                sAtom = tmp;
            } else {
                this->log_.critical(TlUtils::format("atom \"%s\" is not defined. stop.", tmp.c_str()));
                abort();
            }
        }

        // label2
        if (sLine[0] == '@') {
            std::string tmp = TlUtils::getWord(sLine);
            TlUtils::trim_ws(sLine);
            tmp = tmp.substr(1);
            sLabel = tmp;
        }

        // equal
        if (sLine[0] == '=') {
            sLine = sLine.substr(1);
            TlUtils::trim_ws(sLine);
        } else {
            std::cerr << "equal is not found. stop." << std::endl;
            abort();
        }

        // Name
        std::string sName = "";
        if (sLine[0] == '"') {
            sLine = sLine.substr(1);
            std::string::size_type nNameEnd = sLine.find_first_of('"');
            if (nNameEnd != std::string::npos) {
                sName = sLine.substr(0, nNameEnd);
                sLine = (nNameEnd +1 < sLine.length()) ? sLine.substr(nNameEnd +1) : "";
            } else {
                std::cerr << "double quotation is not closed. stop." << std::endl;
                abort();
            }
        } else {
            sName = TlUtils::getWord(sLine);
        }

        // 格納 ============================================================
        // Fl_Db_Basis flDbBasis(sName);
        // for (int i = 0; i < flDbBasis.getTotalcgto(); i++) {
        //     Fl_Gto_Orbital::Cgto cgto;
        //     cgto.basisName = sName;
        //     cgto.Snum = 0;
        //     cgto.Pnum = 0;
        //     cgto.Dnum = 0;
        //     cgto.atom = sAtom;
        //     cgto.label = sLabel;
        //     cgto.shell = flDbBasis.getShell(i);
        //     cgto.scalefactor = flDbBasis.getScalefactor(i);

        //     const int contraction = flDbBasis.getContraction(i);
        //     cgto.pgto.resize(contraction);
        //     for (int j = 0; j < contraction; j++) {
        //         cgto.pgto[j].exponent = flDbBasis.getExpornent(i, j);
        //         cgto.pgto[j].coefficient = flDbBasis.getCoefficient(i ,j);
        //     }

        //     Bas.set(dNumOfCgto, cgto);

        //     dNumOfCgto++;
        // }

        // store to serializeData
        std::string key = sAtom;
        if (sLabel != "") {
            key += "@" + sLabel;
        }
        const TlSerializeData basisset = this->getBasisInfo(sName);
        this->data_["basis_sets_GF"][key] = basisset;
    }
}

void PdfUserInput::alias()
{
    // // xc-potential
    // {
    //     std::string sXcPotential = this->data_["xc-potential"].getStr();
    //     std::string sTilde = "";
    //     if (sXcPotential[sXcPotential.length() -1] == '~') {
    //         sTilde = "~";
    //         sXcPotential = sXcPotential.substr(0, sXcPotential.length() -1);
    //     }

    //     if (TlUtils::toUpper(sXcPotential) == "VWN") {
    //         sXcPotential = "svwn";
    //     }

    //     this->param_["SCF"]["xc-potential"] = sXcPotential + sTilde;
    //     this->data_["xc-potential"] = sXcPotential + sTilde;
    // }
}

bool PdfUserInput::check()
{
    bool bAnswer = true;

    // xc-poteintial
    {
        std::string sXcPotential = this->data_["xc_functional"].getStr();
        std::string sTilde = "";
        if (sXcPotential[sXcPotential.length() -1] == '~') {
            sTilde = "~";
            sXcPotential = sXcPotential.substr(0, sXcPotential.length() -1);
        }

        if (sTilde == "") {
            if ((TlUtils::toUpper(this->data_["scf-memory-saving"].getStr()) == "NO")) {
                std::cout << " 'scf-memory-saving = yes' overridded." << std::endl;
            }
            this->param_["SCF"]["scf-memory-saving"] = "yes";
            this->data_["scf-memory-saving"] = "yes";
        }
    }

    return bAnswer;
}



