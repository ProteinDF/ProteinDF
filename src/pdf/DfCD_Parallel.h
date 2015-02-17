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
#include "TlDistributeSymmetricMatrix.h"
#include "TlRowVectorMatrix.h"
#include "TlColVectorMatrix.h"

class DfCD_Parallel : public DfCD {
public:
    DfCD_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCD_Parallel();

public:
    void calcCholeskyVectorsForJK();
    void calcCholeskyVectorsForGridFree();

    virtual void getJ(TlSymmetricMatrix* pJ);
    void getJ_D(TlDistributeSymmetricMatrix* pJ);

    virtual void getM(const TlSymmetricMatrix& P,
                      TlSymmetricMatrix* pM);

    void getM(const TlDistributeSymmetricMatrix& P,
              TlDistributeSymmetricMatrix* pM);

    // void getJ_distributed(TlDistributeSymmetricMatrix *pJ);
    // void getK_distributed(const RUN_TYPE runType,
    //                       TlDistributeSymmetricMatrix *pK);

protected:
    // void makeSuperMatrix_distribute();
    // TlDistributeSymmetricMatrix 
    // getGMatrix_distribute(const TlOrbitalInfoObject& orbitalInfo, 
    //                       const TlSparseSymmetricMatrix& schwarzTable,
    //                       const index_type numOfItilde,
    //                       const PQ2I_Type& PQ2I);
    // void makeL(const TlDistributeSymmetricMatrix& G);

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;
    virtual void finalize(TlSymmetricMatrix *pMat);
    virtual void finalize(TlSparseSymmetricMatrix *pMat);
    virtual void finalize_I2PQ(PQ_PairArray* pI2PQ);

    virtual void saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath);
    virtual PQ_PairArray getI2PQ(const std::string& filepath);

    // virtual void saveLjk(const TlMatrix& L);
    // virtual TlMatrix getLjk();

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

    TlDistributeSymmetricMatrix 
    getCholeskyVector_distribute(const TlVector& L_col,
                                 const PQ_PairArray& I2PQ);
    TlDistributeSymmetricMatrix 
    getCholeskyVectorA_distribute(const TlOrbitalInfoObject& orbInfo_p,
                                  const TlOrbitalInfoObject& orbInfo_q,
                                  const TlVector& L_col,
                                  const PQ_PairArray& I2PQ);

    // -------------------------------------------------------------------------
protected:
    virtual TlSymmetricMatrix getPMatrix(const RUN_TYPE runType, int itr);

protected:
    virtual void getK_S_woCD(const RUN_TYPE runType,
                             TlSymmetricMatrix* pK);

    virtual void getK_S_fast(const RUN_TYPE runType,
                             TlSymmetricMatrix *pK);

public:
    void getK_D(const RUN_TYPE runType,
                TlDistributeSymmetricMatrix* pK);

protected:
    TlVector getScreenedDensityMatrixD(const PQ_PairArray& I2PQ);
    void expandJMatrixD(const TlVector& vJ, const PQ_PairArray& I2PQ, TlDistributeSymmetricMatrix* pJ);

protected:
    virtual void getM_S(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    virtual void getM_A(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);

    void getM_S(const TlDistributeSymmetricMatrix& P,
                TlDistributeSymmetricMatrix* pM);
    void getM_A(const TlDistributeSymmetricMatrix& P,
                TlDistributeSymmetricMatrix* pM);
    

protected:
    virtual TlRowVectorMatrix 
    calcCholeskyVectorsOnTheFlyS_new(const TlOrbitalInfoObject& orbInfo,
                                     const std::string& I2PQ_path,
                                     const double threshold,
                                     void(DfCD_Parallel::*calcDiagonalsFunc)(
                                         const TlOrbitalInfoObject&,
                                         PQ_PairArray*,
                                         TlVector*),
                                     void(DfCD_Parallel::*getSuperMatrixElements)(
                                         const TlOrbitalInfoObject&,
                                         const index_type,
                                         const std::vector<index_type>&,
                                         const PQ_PairArray&,
                                         std::vector<double>*));

    
protected:
    virtual TlRowVectorMatrix calcCholeskyVectorsOnTheFlyA(const TlOrbitalInfoObject& orbInfo_p,
                                                           const TlOrbitalInfoObject& orbInfo_q,
                                                           const std::string& I2PQ_path);

    // JK --------------------------------------------------------------
    virtual void getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                                        const index_type G_row,
                                        const std::vector<index_type>& G_col_list,
                                        const PQ_PairArray& I2PQ,
                                        std::vector<double>* pSuperMatrixElements);

    // FastCDK-full ----------------------------------------------------
    virtual void getSuperMatrixElements_K_full(const TlOrbitalInfoObject& orbInfo,
                                               const index_type G_row,
                                               const std::vector<index_type>& G_col_list,
                                               const PQ_PairArray& I2PR,
                                               std::vector<double>* pSuperMatrixElements);

    // FastCDK-halt ----------------------------------------------------
    virtual void getSuperMatrixElements_K_half(const TlOrbitalInfoObject& orbInfo,
                                               const index_type G_row,
                                               const std::vector<index_type>& G_col_list,
                                               const PQ_PairArray& I2PR,
                                               std::vector<double>* pSuperMatrixElements);

    // save L matrix ---------------------------------------------------
    void saveL(const TlRowVectorMatrix& L,
               const std::string& path);

    TlColVectorMatrix transLMatrix(const TlRowVectorMatrix& rowVectorMatrix);

    // for debug
    TlMatrix mergeL(const TlRowVectorMatrix& L);
    TlMatrix mergeL(const TlColVectorMatrix& L);

private:
    bool isDebugSaveL_;
};

#endif // DFCD_PARALLEL_H

