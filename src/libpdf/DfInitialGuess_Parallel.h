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

#ifndef DFINITIALGUESS_PARALLEL_H
#define DFINITIALGUESS_PARALLEL_H

#include "DfInitialGuess.h"
#include "tl_dense_general_matrix_scalapack.h"

class DfInitialGuess_Parallel : public DfInitialGuess {
    // ------------------------------------------------------------------------
    // constructor & destructor
    // ------------------------------------------------------------------------
public:
    DfInitialGuess_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess_Parallel();

    // ------------------------------------------------------------------------
    // state
    // ------------------------------------------------------------------------
protected:
    virtual unsigned int loadCalcState() const;
    virtual void saveCalcState(unsigned int cs);

    // ------------------------------------------------------------------------
    // create occupation
    // ------------------------------------------------------------------------
protected:
    virtual void createOccupation(RUN_TYPE runType);

    virtual void createOccupationByFile(const RUN_TYPE runType);

    // ------------------------------------------------------------------------
    // create initial guess (rho)
    // ------------------------------------------------------------------------
protected:
    // virtual void createRho();

    // ------------------------------------------------------------------------
    // create initial guess (Huckel/core)
    // ------------------------------------------------------------------------
protected:
    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();

    // ------------------------------------------------------------------------
    // create initial guess (Harris)
    // ------------------------------------------------------------------------
protected:
    virtual void createInitialGuessUsingHarris();

    // ------------------------------------------------------------------------
    // create initial guess (density matrix)
    // ------------------------------------------------------------------------
protected:

    // ------------------------------------------------------------------------
    // create initial guess (LCAO)
    // ------------------------------------------------------------------------
protected:
    // void createInitialGuessUsingLCAO(); // used in parent class
    virtual void createInitialGuessUsingLCAO(RUN_TYPE runType);

protected:
    void createInitialGuessUsingLCAO_onLAPACK(const RUN_TYPE runType);
    void createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType);

    // TlDenseGeneralMatrix_Scalapack getLCAO_onScaLAPACK(const RUN_TYPE runType);
    // TlDenseGeneralMatrix_Scalapack getLCAO_onScaLAPACK_txt(const RUN_TYPE runType);
    // TlDenseGeneralMatrix_Scalapack getLCAO_onScaLAPACK_bin(const RUN_TYPE runType);

protected:
    virtual DfDmatrix* getDfDmatrixObject(TlSerializeData* param);
};

#endif  // DFINITIALGUESS_PARALLEL_H
