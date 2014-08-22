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

#ifndef FL_GTO_H
#define FL_GTO_H

#include <cassert>
#include <vector>
#include <fstream>
#include <string>
#include "TlSerializeData.h"

class Fl_Gto {
public:
    struct Cgto {
    public:
        Cgto() : Snum(0), Pnum(0), Dnum(0),
                 basisName(""), atom(""), label(""), shellname(""), shell(' '),
                 scalefactor(0.0), pgto() {
        }

        Cgto(const Cgto& rhs) : Snum(rhs.Snum), Pnum(rhs.Pnum), Dnum(rhs.Dnum),
                                basisName(rhs.basisName), atom(rhs.atom),
                                label(rhs.label), shellname(rhs.shellname),
                                shell(rhs.shell), scalefactor(rhs.scalefactor),
                                pgto(rhs.pgto) {
        }

        int contraction() const {
            return this->pgto.size();
        }
        
    public:
        struct Pgto {
        public:
            Pgto() : exponent(0.0), coefficient(0.0) {
            }

            Pgto(const Pgto& rhs) : exponent(rhs.exponent), coefficient(rhs.coefficient) {
            }
            
        public:
            /// orbital exponent of pGTO
            double exponent;
            /// contraction coefficient of pGTO
            double coefficient;
        };
        
        int Snum;
        int Pnum;
        int Dnum;
        std::string basisName;
        std::string atom;  // atom name to which CGTO belongs.
        std::string label; // label2 for atom to which CGTO belongs. this label starts in '@' character.
        std::string shellname;  // shell name, such as 1S, 2S1, 2P and so on.
        char shell;             // shell type, such as s, p, d and so on.
        double scalefactor;     // scale factor for CGTO.

        std::vector<Pgto> pgto;
    };

public:
    explicit Fl_Gto(const std::string& str = "");
    Fl_Gto(const TlSerializeData& data);
    Fl_Gto(const Fl_Gto& rhs);
    
    virtual ~Fl_Gto();

    Fl_Gto operator=(const Fl_Gto& rhs);
    
public:
    // getter ============================================================
    int getSnum(int i) const;
    int getPnum(int i) const;
    int getDnum(int i) const;

    /// 全原子種のCGTOの総数を返す
    int getNumOfCGTOs() const;

    /// CGTOの数を返す
    ///
    /// @param [in] atomSymbol 原子記号
    /// @param [in] label ラベル
//     int getNumOfCGTOs(const std::string& atomSymbol,
//                       const std::string& label) const;

    
//     void getNumOfPGTO(const std::string& atomSymbol,
//                       const std::string& label,
//                       int pGTO_index,
//                       double* pCoef, double* pExp) const;
    
    
    std::string getAtom(int s) const;
    std::string getLabel(int s) const;
    std::string getShellname(int s) const;
    char getShell(int s) const;
    double getScalefactor(int s) const;

    double getNormalizedfactor(int s, int l, int m, int n) const;

    int getContraction(int s) const;
    double getExponent(int s, int p) const;
    double getCoefficient(int s, int p) const;

    double getNormalized(int s, int p, int l, int m, int n) const;

    // y番目が属している基底関数の俗称を与える。
    std::string getBasisName(int y) const;

    // 原子種(ダミー原子含む)の数を返す
    int getBasiskindnumber() const;

    // i番目の原子種(ダミー原子含む)に属するCGTOの数を返す
    int getTermnumber(int i) const;

    // i番目の原子種(ダミー原子含む)に属すCGTOがcgtoの配列の何番目から始まっているかの要素番号を返す
    int getStartposition(int i) const;

    // setter ============================================================
    void set(int i, const Cgto& data);
    void push_back(const Cgto& data);

public:
    void load();
    void save(const std::string& savePath = "");

    void show() const;
    void show(const std::string& title) const;
    std::string getStr_AMOSS() const;
    std::string getStr_GAMESS() const;

public:
    // 規格化定数その２ <1|1/(r1-r2)|2>の規格化タイプ
    double getCoulombnormalized(int k, int p, int l, int m, int n) const;

protected:
    void setup(const TlSerializeData& data);
    
protected:
    //TlSerializeData data_;
    std::vector<Cgto> cgto;
    std::string m_sFileName;
    bool m_isUpdate;
};

#endif // FL_GTO_H
