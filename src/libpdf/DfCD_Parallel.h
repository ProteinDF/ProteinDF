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

#ifndef DFCD_PARALLEL_H
#define DFCD_PARALLEL_H

#include <cstdlib>
#include "DfCD.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"
#include "tl_dense_symmetric_matrix_scalapack.h"

class DfCD_Parallel : public DfCD {
 public:
  DfCD_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfCD_Parallel();

 public:
  virtual void calcCholeskyVectorsForJK();
  virtual void calcCholeskyVectorsForK();
  virtual void calcCholeskyVectorsForGridFree();

  virtual void getJ(TlDenseSymmetricMatrix_Lapack* pJ);
  void getJ_D(TlDenseSymmetricMatrix_Scalapack* pJ);

  virtual void getM(const TlDenseSymmetricMatrix_Lapack& P,
                    TlDenseSymmetricMatrix_Lapack* pM);

  void getM(const TlDenseSymmetricMatrix_Scalapack& P,
            TlDenseSymmetricMatrix_Scalapack* pM);

  // void getJ_distributed(TlDenseSymmetricMatrix_Scalapack *pJ);
  // void getK_distributed(const RUN_TYPE runType,
  //                       TlDenseSymmetricMatrix_Scalapack *pK);

 protected:
  // void makeSuperMatrix_distribute();
  // TlDenseSymmetricMatrix_Scalapack
  // getGMatrix_distribute(const TlOrbitalInfoObject& orbitalInfo,
  //                       const TlSparseSymmetricMatrix& schwarzTable,
  //                       const index_type numOfItilde,
  //                       const PQ2I_Type& PQ2I);
  // void makeL(const TlDenseSymmetricMatrix_Scalapack& G);

  virtual DfTaskCtrl* getDfTaskCtrlObject() const;
  virtual void finalize(TlDenseSymmetricMatrix_Lapack* pMat);
  virtual void finalize(TlSparseSymmetricMatrix* pMat);
  virtual void finalize_I2PQ(PQ_PairArray* pI2PQ);

  virtual void saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath);
  virtual PQ_PairArray getI2PQ(const std::string& filepath);

  // virtual void saveLjk(const TlDenseGeneralMatrix_Lapack& L);
  // virtual TlDenseGeneralMatrix_Lapack getLjk();

  // virtual TlDenseSymmetricMatrix_Lapack getPMatrix();

  virtual void divideCholeskyBasis(const index_type numOfCBs,
                                   index_type* pStart, index_type* pEnd);

  TlDenseSymmetricMatrix_Scalapack getCholeskyVector_distribute(
      const TlDenseVector_Lapack& L_col, const PQ_PairArray& I2PQ);
  TlDenseSymmetricMatrix_Scalapack getCholeskyVectorA_distribute(
      const TlOrbitalInfoObject& orbInfo_p,
      const TlOrbitalInfoObject& orbInfo_q, const TlDenseVector_Lapack& L_col,
      const PQ_PairArray& I2PQ);

  // -------------------------------------------------------------------------
 protected:
  void getJ_cvm(TlDenseSymmetricMatrix_Lapack* pJ);
  void getJ_mmap_DC(TlDenseSymmetricMatrix_Lapack* pJ);

 protected:
  virtual TlDenseSymmetricMatrix_Lapack getPMatrix(const RUN_TYPE runType,
                                                   int itr);

 protected:
  virtual void getK_S_woCD(const RUN_TYPE runType,
                           TlDenseSymmetricMatrix_Lapack* pK);
  virtual void getK_S_woCD_mmap_DC(const RUN_TYPE runType,
                                   TlDenseSymmetricMatrix_Lapack* pK);

  virtual void getK_byLk(const RUN_TYPE runType,
                           TlDenseSymmetricMatrix_Lapack* pK);

 public:
  void getK_D(const RUN_TYPE runType);

 protected:
  void getK_S_woCD_D(const RUN_TYPE runType,
                     TlDenseSymmetricMatrix_Scalapack* pK);

 protected:
  TlDenseVector_Lapack getScreenedDensityMatrixD(const PQ_PairArray& I2PQ);
  void expandJMatrixD(const TlDenseVector_Lapack& vJ, const PQ_PairArray& I2PQ,
                      TlDenseSymmetricMatrix_Scalapack* pJ);

 protected:
  virtual void getM_S(const TlDenseSymmetricMatrix_Lapack& P,
                      TlDenseSymmetricMatrix_Lapack* pM);
  virtual void getM_A(const TlDenseSymmetricMatrix_Lapack& P,
                      TlDenseSymmetricMatrix_Lapack* pM);

  void getM_S(const TlDenseSymmetricMatrix_Scalapack& P,
              TlDenseSymmetricMatrix_Scalapack* pM);
  void getM_A(const TlDenseSymmetricMatrix_Scalapack& P,
              TlDenseSymmetricMatrix_Scalapack* pM);

 protected:
  virtual TlDenseGeneralMatrix_arrays_RowOriented
  calcCholeskyVectorsOnTheFlyS_new(
      const TlOrbitalInfoObject& orbInfo, const std::string& I2PQ_path,
      const double threshold,
      void (DfCD_Parallel::*calcDiagonalsFunc)(const TlOrbitalInfoObject&,
                                               PQ_PairArray*,
                                               std::vector<double>*),
      void (DfCD_Parallel::*getSuperMatrixElements)(
          const TlOrbitalInfoObject&, const index_type,
          const std::vector<index_type>&, const PQ_PairArray&,
          std::vector<double>*));

 protected:
  virtual TlDenseGeneralMatrix_arrays_RowOriented calcCholeskyVectorsOnTheFlyA(
      const TlOrbitalInfoObject& orbInfo_p,
      const TlOrbitalInfoObject& orbInfo_q, const std::string& I2PQ_path);

  // JK --------------------------------------------------------------
  virtual void getSuperMatrixElements(
      const TlOrbitalInfoObject& orbInfo, const index_type G_row,
      const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PQ,
      std::vector<double>* pSuperMatrixElements);

  // FastCDK-full ----------------------------------------------------
  virtual void getSuperMatrixElements_K_full(
      const TlOrbitalInfoObject& orbInfo, const index_type G_row,
      const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PR,
      std::vector<double>* pSuperMatrixElements);

  // FastCDK-halt ----------------------------------------------------
  virtual void getSuperMatrixElements_K_half(
      const TlOrbitalInfoObject& orbInfo, const index_type G_row,
      const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PR,
      std::vector<double>* pSuperMatrixElements);

  // save L matrix ---------------------------------------------------
  void saveL(const TlDenseGeneralMatrix_arrays_RowOriented& L,
             const std::string& path);

  TlDenseGeneralMatrix_arrays_ColOriented transLMatrix(
      const TlDenseGeneralMatrix_arrays_RowOriented& rowVectorMatrix);
  void transLMatrix2mmap(
      const TlDenseGeneralMatrix_arrays_RowOriented& rowVectorMatrix,
      const std::string& savePath);

  // for debug
  TlDenseGeneralMatrix_Lapack mergeL(
      const TlDenseGeneralMatrix_arrays_RowOriented& L);
  TlDenseGeneralMatrix_Lapack mergeL(
      const TlDenseGeneralMatrix_arrays_ColOriented& L);

 private:
  bool isDebugSaveL_;
};

#endif  // DFCD_PARALLEL_H
