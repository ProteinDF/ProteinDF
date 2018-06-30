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

#include <cassert>
#include <set>
#include "DfCalcGridX.h"
#include "TlCommunicate.h"

class TlDenseGeneralMatrix_blacs;
class TlDenseSymmetricMatrix_blacs;

class DfCalcGridX_Parallel : public DfCalcGridX {
 public:
  DfCalcGridX_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfCalcGridX_Parallel();

  // for LAPACK --------------------------------------------------------------
 public:
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                     DfFunctional_LDA* pFunctional,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_A);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                     const TlDenseSymmetricMatrix_BLAS_Old& P_B,
                                     DfFunctional_LDA* pFunctional,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_B);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                     DfFunctional_GGA* pFunctional,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_A);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                                     const TlDenseSymmetricMatrix_BLAS_Old& P_B,
                                     DfFunctional_GGA* pFunctional,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                                     TlDenseSymmetricMatrix_BLAS_Old* pF_B);

 protected:
  virtual void calcRho_LDA(const TlDenseSymmetricMatrix_BLAS_Old& P_A);
  virtual void calcRho_LDA(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                           const TlDenseSymmetricMatrix_BLAS_Old& P_B);
  virtual void calcRho_GGA(const TlDenseSymmetricMatrix_BLAS_Old& P_A);
  virtual void calcRho_GGA(const TlDenseSymmetricMatrix_BLAS_Old& P_A,
                           const TlDenseSymmetricMatrix_BLAS_Old& P_B);

  double buildVxc(DfFunctional_LDA* pFunctional,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_A);
  double buildVxc(DfFunctional_LDA* pFunctional,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_B);
  double buildVxc(DfFunctional_GGA* pFunctional,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_A);
  double buildVxc(DfFunctional_GGA* pFunctional,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_A,
                  TlDenseSymmetricMatrix_BLAS_Old* pF_B);

  TlDenseGeneralMatrix_BLAS_old distributeGridMatrix(const int iteration);
  void gatherAndSaveGridMatrix(const TlDenseGeneralMatrix_BLAS_old& gridMat);

  // for ScaLAPACK -----------------------------------------------------------
 public:
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_blacs& P_A,
                                     DfFunctional_LDA* pFunctional,
                                     TlDenseSymmetricMatrix_blacs* pF_A);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_blacs& P_A,
                                     const TlDenseSymmetricMatrix_blacs& P_B,
                                     DfFunctional_LDA* pFunctional,
                                     TlDenseSymmetricMatrix_blacs* pF_A,
                                     TlDenseSymmetricMatrix_blacs* pF_B);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_blacs& P_A,
                                     DfFunctional_GGA* pFunctional,
                                     TlDenseSymmetricMatrix_blacs* pF_A);
  double calcXCIntegForFockAndEnergy(const TlDenseSymmetricMatrix_blacs& P_A,
                                     const TlDenseSymmetricMatrix_blacs& P_B,
                                     DfFunctional_GGA* pFunctional,
                                     TlDenseSymmetricMatrix_blacs* pF_A,
                                     TlDenseSymmetricMatrix_blacs* pF_B);

 protected:
  TlDenseGeneralMatrix_BLAS_old getGlobalGridMatrix(const int iteration);
  void allReduceGridMatrix(const TlDenseGeneralMatrix_BLAS_old& gridMat);

  void calcRho_LDA(const TlDenseSymmetricMatrix_blacs& P_A);
  void calcRho_LDA(const TlDenseSymmetricMatrix_blacs& P_A,
                   const TlDenseSymmetricMatrix_blacs& P_B);
  void calcRho_GGA(const TlDenseSymmetricMatrix_blacs& P_A);
  void calcRho_GGA(const TlDenseSymmetricMatrix_blacs& P_A,
                   const TlDenseSymmetricMatrix_blacs& P_B);

  void calcRho_LDA(const TlDenseGeneralMatrix_blacs& P_A,
                   TlDenseGeneralMatrix_BLAS_old* pGridMatrix);
  void calcRho_LDA(const TlDenseGeneralMatrix_blacs& P_A,
                   const TlDenseGeneralMatrix_blacs& P_B,
                   TlDenseGeneralMatrix_BLAS_old* pGridMatrix);
  void calcRho_GGA(const TlDenseGeneralMatrix_blacs& P_A,
                   TlDenseGeneralMatrix_BLAS_old* pGridMatrix);
  void calcRho_GGA(const TlDenseGeneralMatrix_blacs& P_A,
                   const TlDenseGeneralMatrix_blacs& P_B,
                   TlDenseGeneralMatrix_BLAS_old* pGridMatrix);

  double buildVxc(DfFunctional_LDA* pFunctional,
                  TlDenseSymmetricMatrix_blacs* pF_A);
  double buildVxc(DfFunctional_LDA* pFunctional,
                  TlDenseSymmetricMatrix_blacs* pF_A,
                  TlDenseSymmetricMatrix_blacs* pF_B);
  double buildVxc(DfFunctional_GGA* pFunctional,
                  TlDenseSymmetricMatrix_blacs* pF_A);
  double buildVxc(DfFunctional_GGA* pFunctional,
                  TlDenseSymmetricMatrix_blacs* pF_A,
                  TlDenseSymmetricMatrix_blacs* pF_B);

  virtual void getWholeDensity(double* pRhoA, double* pRhoB) const;

  virtual void defineCutOffValues(const TlDenseSymmetricMatrix_BLAS_Old& P);

  virtual void defineCutOffValues(const TlDenseSymmetricMatrix_BLAS_Old& PA,
                                  const TlDenseSymmetricMatrix_BLAS_Old& PB);

 public:
  virtual TlDenseGeneralMatrix_BLAS_old energyGradient(
      const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_LDA* pFunctional);
  virtual TlDenseGeneralMatrix_BLAS_old energyGradient(
      const TlDenseSymmetricMatrix_BLAS_Old& P_A, DfFunctional_GGA* pFunctional);

 protected:
  void defineCutOffValues(const TlDenseSymmetricMatrix_blacs& P);
  void defineCutOffValues(const TlDenseSymmetricMatrix_blacs& PA,
                          const TlDenseSymmetricMatrix_blacs& PB);

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

#endif  // DFCALCGRIDX_PARALLEL_H
