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

#ifndef DFGENERATEGRID_PARALLEL_H
#define DFGENERATEGRID_PARALLEL_H

#include "DfGenerateGrid.h"

class DfGenerateGrid_Parallel : public DfGenerateGrid {
public:
    DfGenerateGrid_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfGenerateGrid_Parallel();

protected:
    virtual void logger(const std::string& str) const;

    virtual void makeTable();

    virtual TlMatrix getOMatrix();

    virtual void generateGrid(const TlMatrix& O);
    void generateGrid_DC(const TlMatrix& O);

    // TODO: implement master-slave model
    //void generateGrid_MS();

    void gatherGridData();
    
protected:
    enum {
        TAG_GENGRID_MSG_TO_ROOT = 1201,
        TAG_GENGRID_SEND_RANGE = 1202,
        TAG_GENGRID_SEND_ATOMLIST = 1203,
        TAG_GENGRID_SEND_ATOM = 1204,
        TAG_GENGRID_SEND_DATA = 1205,

        TAG_GENGRID_GATHER_GRID_DATA = 1206
    };
};

#endif // DFGENERATEGRID_PARALLEL_H
