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

#ifndef DFXMATRIX_PARALLEL_H
#define DFXMATRIX_PARALLEL_H

#include "DfXMatrix.h"

class TlDenseGeneralMatrix_Scalapack;
class TlDenseSymmetricMatrix_Scalapack;

class DfXMatrix_Parallel : public DfXMatrix {
   public:
    DfXMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix_Parallel();

   public:
    virtual void buildX();

    virtual void canonicalOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                        TlDenseGeneralMatrix_Lapack* pX,
                                        TlDenseGeneralMatrix_Lapack* pXinv,
                                        const std::string& eigvalFilePath = "");

    virtual void lowdinOrthogonalize(const TlDenseSymmetricMatrix_Lapack& S,
                                     TlDenseGeneralMatrix_Lapack* pX,
                                     TlDenseGeneralMatrix_Lapack* pXinv,
                                     const std::string& eigvalFilePath = "");

    void canonicalOrthogonalize(const TlDenseSymmetricMatrix_Scalapack& S,
                                TlDenseGeneralMatrix_Scalapack* pX,
                                TlDenseGeneralMatrix_Scalapack* pXinv,
                                const std::string& eigvalFilePath = "");

    void lowdinOrthogonalize(const TlDenseSymmetricMatrix_Scalapack& S,
                             TlDenseGeneralMatrix_Scalapack* pX,
                             TlDenseGeneralMatrix_Scalapack* pXinv,
                             const std::string& eigvalFilePath = "");

   protected:
    void buildX_LAPACK();
    void buildX_ScaLAPACK();
};

#endif  // DFXMATRIX_PARALLEL_H
