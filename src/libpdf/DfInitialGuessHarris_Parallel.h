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

#ifndef DFINITIALGUESSHARRIS_PARALLEL_H
#define DFINITIALGUESSHARRIS_PARALLEL_H

#include "DfInitialGuessHarris.h"

/// Harrisの汎関数による初期値を作成する
class DfInitialGuessHarris_Parallel : public DfInitialGuessHarris {
public:
    DfInitialGuessHarris_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHarris_Parallel();

    virtual void main();

protected:
    void distributeHarrisDB();
};

#endif  // DFINITIALGUESSHARRIS_PARALLEL_H
