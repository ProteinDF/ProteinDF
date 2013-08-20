#ifndef DFCD_PARALLEL_H
#define DFCD_PARALLEL_H

#include <cstdlib>
#include "DfCD.h"
#include "TlDistributeSymmetricMatrix.h"
// #include "TlRowVectorMatrix.h"
#include "TlRowVectorMatrix2.h"
#include "TlColVectorMatrix2.h"

class DfCD_Parallel : public DfCD {
public:
    DfCD_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCD_Parallel();

public:
    // virtual void calcCholeskyVectors();
    void getJ_distributed(TlDistributeSymmetricMatrix *pJ);
    void getK_distributed(const RUN_TYPE runType,
                          TlDistributeSymmetricMatrix *pK);

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

    virtual void saveLjk(const TlMatrix& L);
    virtual TlMatrix getLjk();

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

    TlDistributeSymmetricMatrix 
    getCholeskyVector_distribute(const TlVector& L_col,
                                 const PQ_PairArray& I2PQ);

    // -------------------------------------------------------------------------
public:
    virtual void getJ(TlSymmetricMatrix* pJ);
    virtual void getK(const RUN_TYPE runType,
                      TlSymmetricMatrix* pK);

    void getJ_D(TlDistributeSymmetricMatrix* pJ);
    void getK_D(const RUN_TYPE runType,
                TlDistributeSymmetricMatrix* pK);
    
protected:
    virtual TlRowVectorMatrix2 calcCholeskyVectorsOnTheFlyS(const TlOrbitalInfoObject& orbInfo);

    virtual std::vector<double>
    getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo,
                           const index_type G_row,
                           const std::vector<index_type>& G_col_list,
                           const PQ_PairArray& I2PQ,
                           const TlSparseSymmetricMatrix& schwartzTable);
    void saveLjk(const TlRowVectorMatrix2& L);
    // TlColVectorMatrix2 getColVector(const TlRowVectorMatrix2& L);

    // for debug
    TlMatrix mergeL(const TlRowVectorMatrix2& L);
    TlMatrix mergeL(const TlColVectorMatrix2& L);

private:
    bool isDebugSaveL_;
};

#endif // DFCD_PARALLEL_H

