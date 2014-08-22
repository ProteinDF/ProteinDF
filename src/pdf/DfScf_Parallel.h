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

#ifndef DFSCF_PARALLEL_H
#define DFSCF_PARALLEL_H

#include "DfScf.h"

class DfScf_Parallel : public DfScf {
public:
    DfScf_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfScf_Parallel();

protected:
    virtual void logger(const std::string& str) const;

protected:
    virtual DfGridFreeXC* getDfGridFreeXcObject();
    virtual DfXCFunctional* getDfXCFunctional();

protected:
    virtual void saveParam() const;
    virtual void setScfParam();

protected:
    virtual void diffDensityMatrix();

    virtual DfDensityFittingObject* getDfDensityFittingObject();

    virtual void doXCIntegral();

    virtual void doThreeIndexIntegral();

    virtual DfJMatrix* getDfJMatrixObject();
    
    virtual DfKMatrix* getDfKMatrixObject();

    virtual DfFockMatrix* getDfFockMatrixObject();

    virtual DfTransFmatrix* getDfTransFmatrixObject(bool isExecDiis);

    virtual void doLevelShift();

    virtual DfDiagonal* getDfDiagonalObject();

    virtual DfTransatob* getDfTransatobObject();

    virtual DfDmatrix* getDfDmatrixObject();

    virtual DfTotalEnergy* getDfTotalEnergyObject();

    virtual DfPopulation* getDfPopulationObject();

    virtual DfSummary* getDfSummaryObject();

    virtual bool judge();
    virtual int checkConverge();

    virtual DfConverge* getDfConverge();

    virtual void cleanup();

    virtual bool checkMaxIteration();
};

#endif // DFSCF_PARALLEL_H
