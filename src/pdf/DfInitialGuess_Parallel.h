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
#include "TlDistributeMatrix.h"

class DfInitialGuess_Parallel : public DfInitialGuess {
public:
    DfInitialGuess_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess_Parallel();
    
protected:
    // virtual void createRho();
    
    virtual TlVector createOccupation(RUN_TYPE runType);

    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();
    virtual void createInitialGuessUsingHarris();

    virtual void createInitialGuessUsingLCAO(RUN_TYPE runType);

    /// 占有軌道情報を取得する
    virtual TlVector getOccupation(const RUN_TYPE runType);

    /// 占有軌道情報を保存する
    virtual void saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation);

protected:
    void createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType);
    void createInitialGuessUsingLCAO_onLAPACK(const RUN_TYPE runType);

    TlDistributeMatrix getLCAO_onScaLAPACK(const RUN_TYPE runType);

    virtual DfDmatrix* getDfDmatrixObject(TlSerializeData* param);
};

#endif // DFINITIALGUESS_PARALLEL_H
