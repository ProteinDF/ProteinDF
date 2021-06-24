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

#ifndef DFLOCALIZE_PARALLEL_H
#define DFLOCALIZE_PARALLEL_H

#include <vector>

#include "DfLocalize.h"
#include "TlCommunicate.h"

class DfLocalize_Parallel : public DfLocalize {
   public:
    DfLocalize_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfLocalize_Parallel();

   public:
    virtual void exec();
    virtual double localize(TlDenseGeneralMatrix_Lapack* pC);

   protected:
    virtual void initialize();

    // Return true if the task remains, and return false if there are no tasks left.
    // If the molecular orbital in charge is locked during processing,
    // -1 will be placed in the molecular orbitals of TaskItem.
    virtual bool getTaskItem(DfLocalize::TaskItem* pTask, bool isInitialized = false);

   protected:
    void lockMO(const index_type MO);
    void unlockMO(const index_type MO);
    bool isLockedMO(const index_type mo1, const index_type mo2) const;

   protected:
    std::set<DfObject::index_type> lockMOs_;
};

#endif  // DFLOCALIZE_PARALLEL_H
