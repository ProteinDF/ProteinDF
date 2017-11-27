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

#ifndef DFCALCGRIDX_PARALLEL_H
#define DFCALCGRIDX_PARALLEL_H

#include <set>
#include <cassert>
#include "DfCalcGridX.h"
#include "TlCommunicate.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlSparseVectorMatrix.h"

class DfCalcGridX_Parallel : public DfCalcGridX {
public:
    DfCalcGridX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCalcGridX_Parallel();

    // for LAPACK --------------------------------------------------------------
public:
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B,
                                       DfFunctional_LDA* pFunctional,
                                       TlSymmetricMatrix* pF_A,
                                       TlSymmetricMatrix* pF_B);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlSymmetricMatrix& P_A,
                                       const TlSymmetricMatrix& P_B,
                                       DfFunctional_GGA* pFunctional,
                                       TlSymmetricMatrix* pF_A,
                                       TlSymmetricMatrix* pF_B);

protected:
    virtual void calcRho_LDA(const TlSymmetricMatrix& P_A);
    virtual void calcRho_LDA(const TlSymmetricMatrix& P_A,
                             const TlSymmetricMatrix& P_B);
    virtual void calcRho_GGA(const TlSymmetricMatrix& P_A);
    virtual void calcRho_GGA(const TlSymmetricMatrix& P_A,
                             const TlSymmetricMatrix& P_B);

    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlSymmetricMatrix* pF_A,
                    TlSymmetricMatrix* pF_B);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlSymmetricMatrix* pF_A,
                    TlSymmetricMatrix* pF_B);

    TlMatrix distributeGridMatrix(const int iteration);
    void gatherAndSaveGridMatrix(const TlMatrix& gridMat);

    // for ScaLAPACK -----------------------------------------------------------
public:
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       DfFunctional_LDA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B,
                                       DfFunctional_LDA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A,
                                       TlDistributeSymmetricMatrix* pF_B);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       DfFunctional_GGA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A);
    double calcXCIntegForFockAndEnergy(const TlDistributeSymmetricMatrix& P_A,
                                       const TlDistributeSymmetricMatrix& P_B,
                                       DfFunctional_GGA* pFunctional,
                                       TlDistributeSymmetricMatrix* pF_A,
                                       TlDistributeSymmetricMatrix* pF_B);

protected:
    TlMatrix getGlobalGridMatrix(const int iteration);
    void allReduceGridMatrix(const TlMatrix& gridMat);

    void calcRho_LDA(const TlDistributeSymmetricMatrix& P_A);
    void calcRho_LDA(const TlDistributeSymmetricMatrix& P_A,
                     const TlDistributeSymmetricMatrix& P_B);
    void calcRho_GGA(const TlDistributeSymmetricMatrix& P_A);
    void calcRho_GGA(const TlDistributeSymmetricMatrix& P_A,
                     const TlDistributeSymmetricMatrix& P_B);

    void calcRho_LDA(const TlDistributeMatrix& P_A,
                     TlMatrix* pGridMatrix);
    void calcRho_LDA(const TlDistributeMatrix& P_A,
                     const TlDistributeMatrix& P_B,
                     TlMatrix* pGridMatrix);
    void calcRho_GGA(const TlDistributeMatrix& P_A,
                     TlMatrix* pGridMatrix);
    void calcRho_GGA(const TlDistributeMatrix& P_A,
                     const TlDistributeMatrix& P_B,
                     TlMatrix* pGridMatrix);


    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlDistributeSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_LDA* pFunctional,
                    TlDistributeSymmetricMatrix* pF_A,
                    TlDistributeSymmetricMatrix* pF_B);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlDistributeSymmetricMatrix* pF_A);
    double buildVxc(DfFunctional_GGA* pFunctional,
                    TlDistributeSymmetricMatrix* pF_A,
                    TlDistributeSymmetricMatrix* pF_B);

    virtual void getWholeDensity(double* pRhoA, double* pRhoB) const;



    virtual void defineCutOffValues(const TlSymmetricMatrix& P);

    virtual void defineCutOffValues(const TlSymmetricMatrix& PA,
                                    const TlSymmetricMatrix& PB);

public:
    virtual TlMatrix energyGradient(const TlSymmetricMatrix& P_A,
                                    DfFunctional_LDA* pFunctional);
    virtual TlMatrix energyGradient(const TlSymmetricMatrix& P_A,
                                    DfFunctional_GGA* pFunctional);


protected:
    void defineCutOffValues(const TlDistributeSymmetricMatrix& P);
    void defineCutOffValues(const TlDistributeSymmetricMatrix& PA,
                            const TlDistributeSymmetricMatrix& PB);

protected:
    // tag for MPI
    enum {
        TAG_REQUEST_JOB = 9001,
        TAG_ASSIGN_JOB = 9002,
        TAG_TERMINATE_SLAVE = 9003,
        TAG_TERMINATE_OK = 9004,

        TAG_CALC_GRID_DISTRIBUTE = 9005,
        TAG_CALCGRID_GATHER = 9006
    };


protected:
    int assignAtomRange_;
    int assignAoRange_;
    int assignJobsPerProc_;
    std::size_t densityMatrixCacheMemSize_;

    int calcMode_;
};

#endif // DFCALCGRIDX_PARALLEL_H
