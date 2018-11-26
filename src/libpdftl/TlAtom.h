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

#ifndef TLATOM_H
#define TLATOM_H

#include <bitset>
#include <string>
#include "TlPosition.h"
#include "tl_dense_general_matrix_lapack.h"

/// 原子情報を保持するクラス
class TlAtom {
 public:
  /// デフォルトコンストラクタ
  explicit TlAtom(const std::string& symbol = "",
                  const TlPosition& pos = TlPosition(), double charge = 0.0);

  /// デストラクタ
  ~TlAtom();

 public:
  /// 原子を指定する
  void setElement(unsigned int n);
  void setElement(const std::string& symbol);

  /// 指定された文字列の原子番号を返す
  static int getElementNumber(const std::string& symbol);

  /// 原子記号を返す
  inline std::string getSymbol() const;

  double getStdWeight() const;

  /// Van der Waals 半径を返す
  double vdwr() const;

  /// accession
  int getType() const;  // atom type of Miller et al. (1999)

  TlPosition getPosition() const;

  void moveTo(const TlPosition& position);
  void moveTo(const double& x, const double& y, const double& z);
  void shiftBy(const TlPosition& position);
  void shiftBy(const double& x, const double& y, const double& z);
  void rotateBy(const TlDenseGeneralMatrix_Lapack& rot);

  std::string getName() const;
  void setName(const std::string& name);

  void setCharge(double c);
  double getCharge() const;

 public:
  bool operator==(const TlAtom& rhs) const {
    bool answer = false;
    if ((this->m_element == rhs.m_element) &&
        (this->getPosition() == rhs.getPosition()) &&
        (std::fabs(this->getCharge() - rhs.getCharge()) < 1.0E-5)) {
      answer = true;
    }

    return answer;
  }

  bool operator!=(const TlAtom& rhs) const { return !(this->operator==(rhs)); }

 private:
  double charge_;  /// 電荷を保持する

  // physco-chemical properties of Bush & Sheridan (1993)
  std::bitset<16> m_properties;  // 16 = number of properties

  // physco-chemical class of Bush & Sheridan (1993)
  unsigned char m_sPcClass;

  // solvent-accessible surface area
  // double sasa_;

  // atomic symbols in order of element numbers
  static const char* m_sSymbols[];

  // atomic standard weight
  static const double stdAtomicWeight_[];

  // Van der Waals radii of atoms according to their types
  static const double m_dVdw_radii[];

  // position (coordinates) in 3D space
  TlPosition m_position;

  // atom name as defined by the data source;
  // e.g. "CA", "CB", ... , from PDB
  std::string m_name;

  // element number (H=1, He=2, ..., C=6, N=7, O=8, ...)
  int m_element;

  ////////////////////////////////////////////////////////////////////
  // physco-chemical properties od non-hydrogen atoms
  //  using Bush & Sheridan method (1993):
  //  Bush BL, Sheridan RP (1993) J Chem Inf Comput Sci 33: 756-762.
  // uses the following sixteen physco-chemical properties:
  //  # on hybridization
  //    0=sp, 1=sp2, 2=sp3, 3=conj, 4=res,
  //  # on number of neighbors (bonded atoms)
  //    5=x0, 6=x1, 7=x2, 8=x3, 9=x4,
  //  # on aromaticity (5-, 6-membered)
  //    10=ar, 11=ar5, 12=ar6
  //  # electronegativity
  //    13=neg,
  //  # for amide nitrogen
  //    14=namide,
  //  # for carboxylate carbon
  //    15=cx
  //
  //  each atom has 1 or 0 for each of the 16 properties in 'Properties'
  enum Properties {
    sp = 0,
    sp2 = 1,
    sp3 = 2,
    conj = 3,
    res = 4,
    x0 = 5,
    x1 = 6,
    x2 = 7,
    x3 = 8,
    x4 = 9,
    ar = 10,
    ar5 = 11,
    ar6 = 12,
    neg = 13,
    namide = 14,
    cx = 15,
    number_of_properties
  };

  ////////////////////////////////////////////////////////////////////
  // physcochemical classes of atoms:
  //    1=cation, 2=anion, 3=donor, 4=acceptor, 5=polar,
  //    6=hydrophobic, 7=none
  // according to:
  // Bush BL, Sheridan RP (1993) J. Chem. Inf. Comput. Sci 33:756-762.
  enum Classes {
    UNDEFINED = 0,
    CATION = 1,
    ANION = 2,
    DONOR = 3,
    ACCEPTOR = 4,
    POLAR = 5,
    HYDROPHOBIC = 6,
    NONE = 7
  };

  ////////////////////////////////////////////////////////////////////
  // atom types:
  // according to:
  // Miller MD, Sheridan RP, Kearsley SK (1999) J. Med. Chem. 42:1505-1514.
  //
  // 44 = 43 types + 1 (unknown)
  enum Types {
    UNKNOWN = 0,
    Br_sp3_6 = 1,
    Cl_sp3_6 = 2,
    C_sp1_6 = 3,
    C_sp1_7 = 4,
    C_sp2_6 = 5,
    C_sp2_7 = 6,
    C_sp3_6 = 7,
    C_sp3_7 = 8,
    F_sp3_6 = 9,
    I_sp3_6 = 10,
    N_sp1_4 = 11,
    N_sp1_7 = 12,
    N_sp2_1 = 13,
    N_sp2_2 = 14,
    N_sp2_3 = 15,
    N_sp2_4 = 16,
    N_sp2_5 = 17,
    N_sp2_7 = 18,
    N_sp3_1 = 19,
    N_sp3_2 = 20,
    N_sp3_3 = 21,
    N_sp3_5 = 22,
    N_sp3_7 = 23,
    O_sp2_1 = 24,
    O_sp2_2 = 25,
    O_sp2_4 = 26,
    O_sp2_5 = 27,
    O_sp3_1 = 28,
    O_sp3_2 = 29,
    O_sp3_4 = 30,
    O_sp3_5 = 31,
    P_sp2_7 = 32,
    P_sp3_1 = 33,
    P_sp3_7 = 34,
    S_sp1_7 = 35,
    S_sp2_2 = 36,
    S_sp2_4 = 37,
    S_sp2_6 = 38,
    S_sp2_7 = 39,
    S_sp3_1 = 40,
    S_sp3_2 = 41,
    S_sp3_6 = 42,
    S_sp3_7 = 43
  };
};

// inline member functions
inline std::string TlAtom::getSymbol() const {
  return std::string(TlAtom::m_sSymbols[this->m_element]);
}

inline double TlAtom::getStdWeight() const {
  return stdAtomicWeight_[this->m_element];
}

inline double TlAtom::vdwr() const { return m_dVdw_radii[this->getType()]; }

inline TlPosition TlAtom::getPosition() const { return m_position; }

inline void TlAtom::moveTo(const TlPosition& position) {
  this->m_position.moveTo(position);
}

inline void TlAtom::moveTo(const double& x, const double& y, const double& z) {
  this->m_position.moveTo(x, y, z);
}

inline void TlAtom::shiftBy(const TlPosition& position) {
  this->m_position.shiftBy(position);
}

inline void TlAtom::shiftBy(const double& x, const double& y, const double& z) {
  this->m_position.shiftBy(x, y, z);
}

inline void TlAtom::rotateBy(const TlDenseGeneralMatrix_Lapack& rot) {
  this->m_position.rotateBy(rot);
}

inline std::string TlAtom::getName() const { return this->m_name; }

inline void TlAtom::setName(const std::string& name) { this->m_name = name; }

#endif  // TLATOM_H
