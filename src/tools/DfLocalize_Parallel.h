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

#include <utility>
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

   protected:
    void getBlockCMatrix(const TlDenseGeneralMatrix_Lapack& C, const index_type MOsPerBlock, const int block1,
                         const int block2, TlDenseGeneralMatrix_Lapack* pBlockC);
    void setBlockCMatrix(const index_type MOsPerBlock, const int block1, const int block2,
                         const TlDenseGeneralMatrix_Lapack& blockC, TlDenseGeneralMatrix_Lapack* pC);

   protected:
    typedef std::pair<int, int> JobItem;  // contains blockId1, blockId2
    void makeJobList(const int dim);
    bool getJobItem(JobItem* pJob, bool isInitialized = false);

    std::vector<JobItem> jobList_;

   protected:
    void initLockBlock(const int block);
    void lockBlock(const int block);
    void unlockBlock(const int block);
    bool isLockedBlock(const int block1, const int block2) const;

   protected:
    std::vector<char> lockBlocks_;
};

#endif  // DFLOCALIZE_PARALLEL_H
