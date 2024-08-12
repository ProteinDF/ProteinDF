#ifndef DF_CDK_MATRIX_H
#define DF_CDK_MATRIX_H

#include "DfObject.h"
#include "common.h"
#include "datatype.h"

class DfTaskCtrl;

class DfCdkMatrix : public DfObject {
public:
    DfCdkMatrix(TlSerializeData* pPdfParam);
    virtual ~DfCdkMatrix();

public:
    virtual void getK();

protected:
    enum FastCDK_MODE {
        FASTCDK_NONE,
        FASTCDK_PRODUCTIVE,
        FASTCDK_PRODUCTIVE_FULL,
        FASTCDK_DEBUG_FULL_SUPERMATRIX,  // V_pr,qs = <pq|rs>
        FASTCDK_DEBUG_SUPERMATRIX        // V_pr,qs = <pq|rs> + <ps|rq>
    };

protected:
    template <typename SymmetricMatrix, typename Vector,
              typename SparseGeneralMatrix, typename SparseSymmetricMatrix,
              typename TempSymmetricMatrix = SymmetricMatrix>
    void getK_runType();

    template <typename SymmetricMatrix, typename Vector, typename SparseGeneralMatrix, typename SparseSymmetricMatrix, typename TempSymmetricMatrix = SymmetricMatrix>
    void getK_runType_L(DfObject::RUN_TYPE runType);

    // ---------------------------------------------------------------------------
    // calc HFx by Ljk
    // ---------------------------------------------------------------------------
protected:
    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix = SymmetricMatrix>
    void getK_byLjk_useDenseMatrix(const RUN_TYPE runType);

    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename SparseGeneralMatrix>
    void getK_byLjk_useTransMatrix(const RUN_TYPE runType);

    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
    void getK_byLjk_useSparseMatrix(const RUN_TYPE runType);

    // ---------------------------------------------------------------------------
    // calc HFx by Lk
    // ---------------------------------------------------------------------------
protected:
    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector>
    void getK_byLk(const RUN_TYPE runType);

    // ---------------------------------------------------------------------------
    // task control
    // ---------------------------------------------------------------------------
protected:
    DfTaskCtrl* getDfTaskCtrlObject() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    PQ_PairArray getI2PQ(const std::string& filepath);

    template <typename SymmetricMatrix>
    void finalize(SymmetricMatrix* pK);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
protected:
    template <typename SymmetricMatrix, typename Vector>
    SymmetricMatrix getCholeskyVector(const Vector& L_col, const PQ_PairArray& I2PQ);

    template <typename SparseGeneralMatrix>
    SparseGeneralMatrix getTrans_I2PQ_Matrix(const PQ_PairArray& I2PQ);

    template <typename SymmetricMatrix, typename Vector, typename SparseGeneralMatrix>
    SymmetricMatrix convert_I2PQ(const SparseGeneralMatrix& I2PQ_mat, const Vector& L);

    // ---------------------------------------------------------------------------
    template <typename SymmetricMatrix, typename Vector>
    Vector getScreenedDensityMatrix(const RUN_TYPE runType,
                                    const PQ_PairArray& I2PR);

    template <typename SymmetricMatrix, typename Vector>
    void expandKMatrix(const Vector& vK, const PQ_PairArray& I2PR, SymmetricMatrix* pK);

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    bool useMmapMatrix_;
    bool isCvSavedAsMmap_;
    int sparseMatrixLevel_;

    FastCDK_MODE fastCDK_mode_;

    bool useFP32_;
};

#endif  // DF_CDK_MATRIX_H
