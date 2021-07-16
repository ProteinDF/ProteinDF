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

#ifndef DFINTEGRALS_H
#define DFINTEGRALS_H

#include "DfObject.h"
#include "TlSerializeData.h"

class DfHpq;
class DfOverlap;
class DfEri;
class DfXMatrix;
class DfInvMatrix;
class DfCD;
class DfGridFreeXC;
class DfGenerateGrid;
class TlDenseSymmetricMatrix_Lapack;

/// 1 電子ハミルトニアン, 2 中心積分、1 中心積分(Na)を計算するクラス
class DfIntegrals : public DfObject {
   protected:
    enum CalcState {
        Hpq = 1,
        Spq = 2,
        Sab2 = 4,
        Sgd = 8,
        Na = 16,
        Sab = 32,
        CD = 64,
        CDK = 128,
        X = 256,
        INV = 512,
        CHOLESKY_VECTORS_XC = 1024,
        GRID_FREE = 2048,
        GRID = 4096
    };

   public:
    explicit DfIntegrals(TlSerializeData* param = NULL);
    virtual ~DfIntegrals();

   public:
    void main();

   protected:
    virtual void saveParam();

   protected:
    virtual void createHpqMatrix();
    virtual void createOverlapMatrix();
    virtual void createERIMatrix();
    void createCholeskyVectors();
    void createCholeskyVectors_K();

    void createInverseMatrixes();
    void createXMatrix();
    void createCholeskyVectors_XC();
    void prepareGridFree();
    void createGrids();

   protected:
    virtual DfXMatrix* getDfXMatrixObject();
    virtual DfInvMatrix* getDfInvMatrixObject();
    virtual DfCD* getDfCDObject();
    virtual DfGridFreeXC* getDfGridFreeXCObject();
    virtual DfGenerateGrid* getDfGenerateGridObject();

   protected:
    virtual void saveInvSquareVMatrix(const TlDenseSymmetricMatrix_Lapack& v);

   protected:
    virtual void outputStartTitle(const std::string& stepName, const char lineChar = '-');
    virtual void outputEndTitle(const std::string& stepName = "", const char lineChar = '-');

   protected:
    std::string saveParamPath_;
};

#endif  // DFINTEGRALS_H
