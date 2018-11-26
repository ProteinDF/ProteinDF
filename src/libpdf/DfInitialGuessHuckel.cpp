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

#include "DfInitialGuessHuckel.h"
#include <iostream>
#include <limits>
#include "CnError.h"
#include "tl_dense_symmetric_matrix_lapack.h"

const double DfInitialGuessHuckel::EPS = std::numeric_limits<double>::epsilon();

DfInitialGuessHuckel::DfInitialGuessHuckel(TlSerializeData* pPdfParam)
    : DfInitialGuess(pPdfParam) {
  this->initialize();
}

DfInitialGuessHuckel::~DfInitialGuessHuckel() {}

void DfInitialGuessHuckel::initialize() {
  this->m_Hii.clear();

  this->m_Hii["H"][0] = -13.6;  // H
  this->m_Hii["C"][0] = -21.4;  // C
  this->m_Hii["C"][1] = -11.4;
  this->m_Hii["N"][0] = -26.0;  // N
  this->m_Hii["N"][1] = -13.4;
  this->m_Hii["O"][0] = -32.3;  // O
  this->m_Hii["O"][1] = -14.8;
  this->m_Hii["S"][0] = -20.0;  // S
  this->m_Hii["S"][1] = -11.0;
  this->m_Hii["S"][2] = -8.0;
}

void DfInitialGuessHuckel::createGuess() {
  this->createGuess<TlDenseSymmetricMatrix_Lapack,
                    TlDenseGeneralMatrix_Lapack>();
}

double DfInitialGuessHuckel::getHii(const std::string& sAtomName,
                                    const int nOrbitalType) {
  double dAnswer = 0.0;

  std::map<std::string, std::map<int, double> >::const_iterator p =
      this->m_Hii.find(sAtomName);
  if (p != this->m_Hii.end()) {
    std::map<int, double>::const_iterator q = p->second.find(nOrbitalType);
    if (q != p->second.end()) {
      dAnswer = q->second;
    }
  }

  return dAnswer;
}
