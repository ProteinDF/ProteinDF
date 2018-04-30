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

#include "TlRotateLCAO.h"
#include <iostream>

TlRotateLCAO::TlRotateLCAO(const TlOrbitalInfo& orbInfo) : orbInfo_(orbInfo) {}

TlRotateLCAO::~TlRotateLCAO() {}

TlMatrix TlRotateLCAO::exec(const TlMatrix& lcao, const TlMatrix& rot) {
  TlMatrix answer(lcao.getNumOfRows(), lcao.getNumOfCols());
  if ((rot.getNumOfRows() != 3) || (rot.getNumOfCols() != 3)) {
    return answer;
  }

  const int maxRow = lcao.getNumOfRows();
  const int maxCol = lcao.getNumOfCols();
  for (int row = 0; row < maxRow;) {
    const int shellType = this->orbInfo_.getShellType(row);
    switch (shellType) {
      case 0:  // s
        this->rotateLCAO_typeS(lcao, rot, row, maxCol, &answer);
        break;

      case 1:  // p
        this->rotateLCAO_typeP(lcao, rot, row, maxCol, &answer);
        break;

      case 2:  // d
        this->rotateLCAO_typeD(lcao, rot, row, maxCol, &answer);
        break;
    }
    row += (shellType * 2 + 1);
  }

  return answer;
}

void TlRotateLCAO::rotateLCAO_typeS(const TlMatrix& lcao, const TlMatrix& rot,
                                    const int row, const int maxCol,
                                    TlMatrix* ioMatrix) {
  for (int col = 0; col < maxCol; ++col) {
    ioMatrix->set(row, col, lcao.get(row, col));
  }
}

void TlRotateLCAO::rotateLCAO_typeP(const TlMatrix& lcao, const TlMatrix& rot,
                                    const int row, const int maxCol,
                                    TlMatrix* ioMatrix) {
  for (int col = 0; col < maxCol; ++col) {
    TlVector v(3);
    for (int i = 0; i < 3; ++i) {
      v[i] = lcao.get(row + i, col);
    }

    const TlVector u = v * rot;

    for (int i = 0; i < 3; ++i) {
      ioMatrix->set(row + i, col, u.get(i));
    }
  }
}

void TlRotateLCAO::rotateLCAO_typeD(const TlMatrix& lcao, const TlMatrix& rot,
                                    const int row, const int maxCol,
                                    TlMatrix* ioMatrix) {
  const double Cxx = rot.get(0, 0);
  const double Cxy = rot.get(0, 1);
  const double Cxz = rot.get(0, 2);
  const double Cyx = rot.get(1, 0);
  const double Cyy = rot.get(1, 1);
  const double Cyz = rot.get(2, 2);
  const double Czx = rot.get(3, 0);
  const double Czy = rot.get(3, 1);
  const double Czz = rot.get(3, 2);
  static const double INV_SQRT3 = 1.0 / std::sqrt(3.0);
  static const double INV_6 = 1.0 / 6.0;
  static const double COEF2 = 1.0 / 3.0 * std::sqrt(3.0);
  static const double COEF = COEF2 * 0.5;
  const double normFact = 1.0;  // to set

  for (int col = 0; col < maxCol; ++col) {
    const double dxy = lcao.get(row, col);
    const double dxz = lcao.get(row + 1, col);
    const double dyz = lcao.get(row + 2, col);
    const double dx2y2 = lcao.get(row + 3, col);
    const double dz2 = lcao.get(row + 4, col);

    const double dxy_new =
        (Cxx * Cyy + Cxy * Cyx) * dxy + (Cxx * Czy + Cxy * Czx) * dxz +
        (Cyx * Czy + Cyy * Czx) * dyz + (Cxx * Cxy + Cyx * Cyy) * dx2y2 -
        (Cxx * Cxy + Cyx * Cyy - 2.0 * Czx * Czy) * dz2 * COEF2;
    const double dxz_new =
        (Cxx * Cyz + Cxz * Cyx) * dxy + (Cxx * Czz + Cxz * Czx) * dxz +
        (Cyx * Czz + Cyz * Czx) * dyz + (Cxx * Cxz - Cyx * Cyz) * dx2y2 +
        (Cxx * Cxz + Cyx * Cyz - 2.0 * Czx * Czz) * dz2 * COEF2;
    const double dyz_new =
        (Cxy * Cyz + Cxz * Cyy) * dxy + (Cxy * Czz + Cxz * Czy) * dxz +
        (Cyy * Czz + Cyz * Czy) * dyz + (Cxy * Cxz + Cyy * Cyz) * dx2y2 +
        (Cxy * Cxz + Cyy * Cyz - 2.0 * Czy * Czz) * dz2 * COEF2;
    const double dx2y2_new =
        (Cxx * Cyx - Cxy * Cyy) * dxy + (Cxx * Czx - Cxy * Czy) * dxz +
        (Cyx * Czx - Cyy * Czy) * dyz * normFact +
        (Cxx * Cxx - Cyx * Cyx - Cxy * Cxy + Cyy * Cyy) * dx2y2 * 0.5 -
        (Cxx * Cxx + Cyx * Cyx - 2.0 * Czx * Czx - Cxy * Cxy - Cyy * Cyy +
         2.0 * Czy * Czy) *
            dz2 * COEF;
    const double dz2_new =
        (2.0 * Cxz * Cyz - Cxx * Cyx - Cxy * Cyy) * dxy * INV_SQRT3 +
        (2.0 * Cxz * Czz - Cxx * Czx - Cxy * Czy) * dxy * INV_SQRT3 +
        (2.0 * Cyz * Czz - Cxx * Czx - Cxy * Czy) * dxy * INV_SQRT3 +
        (2.0 * Cxz * Cxz - 2.0 * Cyz * Cyz - Cxx * Cxx + Cyx * Cyx - Cxy * Cxy +
         Cyy * Cyy) *
            dx2y2 * 0.5 * INV_SQRT3 -
        (2.0 * Cxz * Cxz + 2.0 * Cyz * Cyz - 4.0 * Czz * Czz - Cxx * Cxx -
         Cyx * Cyx + 2.0 * Czx * Czx - Cxy * Cxy - Cyy * Cyy +
         2.0 * Czy * Czy) *
            dz2 * INV_6;

    ioMatrix->set(row, col, dxy_new);
    ioMatrix->set(row + 1, col, dxz_new);
    ioMatrix->set(row + 2, col, dyz_new);
    ioMatrix->set(row + 3, col, dx2y2_new);
    ioMatrix->set(row + 4, col, dz2_new);
  }
}
