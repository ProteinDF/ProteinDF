#ifndef DF_CDK_MATRIX_H
#define DF_CDK_MATRIX_H

#include "DfObject.h"
#include "common.h"
#include "datatype.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"
#include "tl_dense_general_matrix_mmap.h"

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
              typename TempDenseSymmetricMatrix>
    void getK_runType();

    template <typename DenseSymmetricMatrix, typename Vector, typename SparseGeneralMatrix, typename SparseSymmetricMatrix,
              typename TempDenseSymmetricMatrix = DenseSymmetricMatrix>
    void getK_runType_L(DfObject::RUN_TYPE runType);

    // ---------------------------------------------------------------------------
    // calc HFx by Ljk
    // ---------------------------------------------------------------------------
protected:
    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix>
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

// ----------------------------------------------------------------------------
// template
// ----------------------------------------------------------------------------
template <typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix, typename SparseSymmetricMatrix,
          typename TempDenseSymmetricMatrix>
void DfCdkMatrix::getK_runType() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            this->getK_runType_L<SymmetricMatrix, Vector, SparseGeneralMatrix, SparseSymmetricMatrix, TempDenseSymmetricMatrix>(RUN_RKS);
        } break;

        case METHOD_UKS: {
            this->getK_runType_L<SymmetricMatrix, Vector, SparseGeneralMatrix, SparseSymmetricMatrix, TempDenseSymmetricMatrix>(RUN_UKS_ALPHA);
            this->getK_runType_L<SymmetricMatrix, Vector, SparseGeneralMatrix, SparseSymmetricMatrix, TempDenseSymmetricMatrix>(RUN_UKS_BETA);

        } break;

        case METHOD_ROKS: {
            this->getK_runType_L<SymmetricMatrix, Vector, SparseGeneralMatrix, SparseSymmetricMatrix, TempDenseSymmetricMatrix>(RUN_ROKS_ALPHA);
            this->getK_runType_L<SymmetricMatrix, Vector, SparseGeneralMatrix, SparseSymmetricMatrix, TempDenseSymmetricMatrix>(RUN_ROKS_BETA);
        } break;

        default:
            this->log_.critical("program error.");
            CnErr.abort();
            break;
    }
}

template <typename DenseSymmetricMatrix, typename Vector, typename SparseGeneralMatrix, typename SparseSymmetricMatrix, typename TempDenseSymmetricMatrix>
void DfCdkMatrix::getK_runType_L(DfObject::RUN_TYPE runType) {
    switch (this->fastCDK_mode_) {
        case FASTCDK_NONE:
            if (this->isCvSavedAsMmap_) {
                this->log_.info("loading the cholesky vectors on mmap.");
                switch (this->sparseMatrixLevel_) {
                    case 1: {
                        this->log_.info("use sparse matrix to transpose the cholesky vector.");
                        this->getK_byLjk_useTransMatrix<TlDenseGeneralMatrix_mmap, DenseSymmetricMatrix, Vector, SparseGeneralMatrix>(runType);
                    } break;

                    case 2: {
                        this->log_.info("use sparse matrix to transpose the cholesky vector.");
                        this->log_.info("use sparse matrix in temporary K matrix.");
                        this->getK_byLjk_useSparseMatrix<TlDenseGeneralMatrix_mmap, DenseSymmetricMatrix, Vector,
                                                         SparseGeneralMatrix, SparseSymmetricMatrix>(runType);

                    } break;

                    default: {
                        this->log_.info("use dense matrix for building K.");
                        this->getK_byLjk_useDenseMatrix<TlDenseGeneralMatrix_mmap, DenseSymmetricMatrix, Vector, TempDenseSymmetricMatrix>(runType);
                    } break;
                }

            } else {
                this->log_.info("loading the cholesky vectors on arrays.");
                // v1
                this->getK_byLjk_useDenseMatrix<TlDenseGeneralMatrix_arrays_ColOriented, DenseSymmetricMatrix, Vector, TempDenseSymmetricMatrix>(runType);
            }
            break;

        case FASTCDK_DEBUG_FULL_SUPERMATRIX:
        case FASTCDK_DEBUG_SUPERMATRIX:
        case FASTCDK_PRODUCTIVE_FULL:
        case FASTCDK_PRODUCTIVE:
            this->getK_byLk<TlDenseGeneralMatrix_mmap, DenseSymmetricMatrix, Vector>(runType);
            break;

        default:
            this->log_.critical(
                TlUtils::format("%s: %d: program error.", __FILE__, __LINE__));
            CnErr.abort();
    }
}

template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix>
void DfCdkMatrix::getK_byLjk_useDenseMatrix(const RUN_TYPE runType) {
    this->log_.info("calc K by CD method");
    this->log_.info("transformed by dens matrix");

    const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
    this->log_.info(TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    index_type calcNumOfCVs = L.getNumOfCols();

    TempSymmetricMatrix P;
    if ((this->isUpdateMethod_ == true) && (this->m_nIteration > 1)) {
        this->log_.info("update method is applied.");
        P = DfObject::getDiffDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration);

        double threshold = 1.0E-10;
        const double diffP = (*this->pPdfParam_)["DfDiffDensityMatrix"]["max_abs"][this->m_nIteration - 1].getDouble();
        this->log_.info(TlUtils::format("dP max: %e", diffP));
        if (diffP < 1.0) {
            threshold /= diffP;
            threshold = std::sqrt(threshold);
        }
        this->log_.info(TlUtils::format("threshold: %e", threshold));

        std::vector<double> errors;
        {
            TlDenseVector_Lapack errs_vtr = DfObject::loadLjkErrorsVector<TlDenseVector_Lapack>();
            errors = errs_vtr;
        }
        std::vector<double>::iterator it = std::find_if(errors.begin(), errors.end(),
                                                        std::bind(std::less<double>(), std::placeholders::_1, threshold));
        calcNumOfCVs = static_cast<index_type>(std::distance(errors.begin(), it));

    } else {
        P = this->getPInMatrix<SymmetricMatrix>(runType, this->m_nIteration);
    }
    if (runType == RUN_RKS) {
        P *= 0.5;
    }

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

    int numOfThreads = 1;
#ifdef _OPENMP
    numOfThreads = omp_get_num_procs();
#endif  // _OPENMP
    const int taskSize = 100 * numOfThreads;
    this->log_.info(TlUtils::format("task size: %d", taskSize));

    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<std::size_t> tasks;
    bool hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks, true);

    TempSymmetricMatrix K(this->m_nNumOfAOs);
    while (hasTasks == true) {
        const int numOfTasks = tasks.size();
        // #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < numOfTasks; ++i) {
            const TempSymmetricMatrix l = this->getCholeskyVector<SymmetricMatrix, Vector>(L.getColVector(tasks[i]), I2PQ);
            // this->log_.info(TlUtils::format("l: %s", TlUtils::getTypeName(l).c_str()));
            assert(l.getNumOfRows() == this->m_nNumOfAOs);

            // #pragma omp critical
            { K += l * P * l; }
        }
        hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks);
    }

    K *= -1.0;

    SymmetricMatrix Kout = K;
    if ((this->isUpdateMethod_ == true) && (this->m_nIteration > 1)) {
        SymmetricMatrix prevK = this->getHFxMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
        Kout += prevK;
    }

    this->log_.info("finalize");
    // this->finalize(&K);

    delete pTaskCtrl;
    pTaskCtrl = NULL;

    this->saveHFxMatrix(runType, this->m_nIteration, Kout);
}

template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix>
void DfCdkMatrix::getK_byLjk_useTransMatrix(const RUN_TYPE runType) {
    this->log_.info("calc K by CD method");
    this->log_.info("transformed by trans matrix");
    SymmetricMatrix K(this->m_nNumOfAOs);

    const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
    this->log_.info(
        TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const SymmetricMatrix P = this->getSpinDensityMatrix<SymmetricMatrix>(
        runType, this->m_nIteration - 1);

    this->log_.info("start loop");
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const SparseGeneralMatrix I2PQ_mat =
        this->getTrans_I2PQ_Matrix<SparseGeneralMatrix>(I2PQ);

    int numOfThreads = 1;
#ifdef _OPENMP
    numOfThreads = omp_get_num_procs();
#endif  // _OPENMP
    const int taskSize = 100 * numOfThreads;
    this->log_.info(TlUtils::format("task size: %d", taskSize));

    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<std::size_t> tasks;
    bool hasTasks = pTaskCtrl->getQueue(numOfCBs, taskSize, &tasks, true);
    while (hasTasks == true) {
        const int numOfTasks = tasks.size();
        // #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < numOfTasks; ++i) {
            const SymmetricMatrix l = this->convert_I2PQ<SymmetricMatrix, Vector, SparseGeneralMatrix>(I2PQ_mat, L.getColVector(tasks[i]));
            assert(l.getNumOfRows() == this->m_nNumOfAOs);

            // #pragma omp critical
            { K += l * P * l; }
        }

        hasTasks = pTaskCtrl->getQueue(numOfCBs, taskSize, &tasks);
    }

    K *= -1.0;
    this->log_.info("finalize");
    // this->finalize(&K);

    delete pTaskCtrl;
    pTaskCtrl = NULL;

    this->saveHFxMatrix(runType, this->m_nIteration, K);
}

template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
void DfCdkMatrix::getK_byLjk_useSparseMatrix(const RUN_TYPE runType) {
    this->log_.info("calc K by CD method");
    this->log_.info("transformed by sparse matrix");
    SymmetricMatrix K(this->m_nNumOfAOs);

    const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
    this->log_.info(
        TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const SymmetricMatrix P = this->getSpinDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
    const SparseSymmetricMatrix P_SM = P;

    this->log_.info("start loop");
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const SparseGeneralMatrix I2PQ_mat =
        this->getTrans_I2PQ_Matrix<SparseGeneralMatrix>(I2PQ);

    SparseSymmetricMatrix tmpK(this->m_nNumOfAOs);
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<std::size_t> tasks;

    int numOfThreads = 1;
#ifdef _OPENMP
    numOfThreads = omp_get_num_procs();
#endif  // _OPENMP
    const int taskSize = 100 * numOfThreads;
    this->log_.info(TlUtils::format("task size: %d", taskSize));

    bool hasTasks = pTaskCtrl->getQueue(numOfCBs, taskSize, &tasks, true);
    while (hasTasks == true) {
        const int numOfTasks = tasks.size();
        // #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < numOfTasks; ++i) {
            const SparseSymmetricMatrix l =
                this->convert_I2PQ<SymmetricMatrix, Vector,
                                   SparseGeneralMatrix>(
                    I2PQ_mat, L.getColVector(tasks[i]));
            assert(l.getNumOfRows() == this->m_nNumOfAOs);

            SparseSymmetricMatrix X = l * P_SM * l;
            // #pragma omp critical
            { tmpK += X; }
        }

        hasTasks = pTaskCtrl->getQueue(numOfCBs, taskSize, &tasks);
    }

    K = tmpK;
    K *= -1.0;
    this->log_.info("finalize");
    // this->finalize(&K);

    delete pTaskCtrl;
    pTaskCtrl = NULL;

    this->saveHFxMatrix(runType, this->m_nIteration, K);
}

// -----------------------------------------------------------------------------
// getTrans_I2PQ_Matrix
// -----------------------------------------------------------------------------
template <typename SparseGeneralMatrix>
SparseGeneralMatrix DfCdkMatrix::getTrans_I2PQ_Matrix(
    const PQ_PairArray& I2PQ) {
    const TlMatrixObject::index_type numOfItilde = I2PQ.size();
    const TlMatrixObject::index_type numOfAoPairs =
        this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    SparseGeneralMatrix expandL(numOfAoPairs, numOfItilde);
    for (TlMatrixObject::index_type i = 0; i < numOfItilde; ++i) {
        TlMatrixObject::index_type row = I2PQ[i].index1();
        TlMatrixObject::index_type col = I2PQ[i].index2();
        if (row < col) {
            std::swap(row, col);
        }
        const TlMatrixObject::index_type rowcol = row * (row + 1) / 2 + col;
        expandL.set(rowcol, i, 1.0);
    }
    return expandL;
}

// -----------------------------------------------------------------------------
// convert_I2PQ
// -----------------------------------------------------------------------------
template <typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix>
SymmetricMatrix DfCdkMatrix::convert_I2PQ(
    const SparseGeneralMatrix& I2PQ_mat, const Vector& L) {
    Vector expandL = I2PQ_mat * L;
    SymmetricMatrix answer(this->m_nNumOfAOs);
    answer.vtr2mat(expandL);

    return answer;
}

// -----------------------------------------------------------------------------
// convert cholesky vector
// -----------------------------------------------------------------------------
template <typename SymmetricMatrix, typename Vector>
SymmetricMatrix DfCdkMatrix::getCholeskyVector(const Vector& L_col,
                                               const PQ_PairArray& I2PQ) {
    const index_type numOfItilde = L_col.getSize();
    assert(static_cast<std::size_t>(numOfItilde) == I2PQ.size());

    SymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(), I2PQ[i].index2(), L_col.get(i));
    }

    return answer;
}

// -----------------------------------------------------------------------------
// Fast CDK
// -----------------------------------------------------------------------------
template <typename Lk_MatrixType, typename SymmetricMatrix, typename Vector>
void DfCdkMatrix::getK_byLk(const RUN_TYPE runType) {
    this->log_.info("calc K(fast) by CD method (serial).");
    SymmetricMatrix K(this->m_nNumOfAOs);

    const Lk_MatrixType L(DfObject::getLkMatrixPath());
    this->log_.info(
        TlUtils::format("L(J): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const PQ_PairArray I2PR = this->getI2PQ(this->getI2prVtrPath());
    const Vector vP = this->getScreenedDensityMatrix<SymmetricMatrix, Vector>(runType, I2PR);

    const index_type numOfI = I2PR.size();
    Vector vK(numOfI);

    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<std::size_t> tasks;
    bool hasTasks = pTaskCtrl->getQueue(numOfCBs, 100, &tasks, true);

    this->log_.info("start loop");
    while (hasTasks == true) {
        std::vector<std::size_t>::const_iterator itEnd = tasks.end();
        for (std::vector<std::size_t>::const_iterator it = tasks.begin();
             it != itEnd; ++it) {
            const Vector LI = L.getColVector(*it);
            assert(LI.getSize() == vK.getSize());

            Vector tmpLI = LI;
            const double qi = tmpLI.dotInPlace(vP).sum();
            vK += qi * LI;
        }

        hasTasks = pTaskCtrl->getQueue(numOfCBs, 100, &tasks);
    }

    vK *= -1.0;

    this->expandKMatrix(vK, I2PR, &K);
    // this->finalize(&K);

    this->saveHFxMatrix(runType, this->m_nIteration, K);
}

template <typename SymmetricMatrix, typename Vector>
void DfCdkMatrix::expandKMatrix(const Vector& vK, const PQ_PairArray& I2PR,
                                SymmetricMatrix* pK) {
    assert(pK != NULL);
    const index_type numOfI = I2PR.size();
    for (index_type i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PR[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        const double coef = (r != c) ? 0.5 : 1.0;
        pK->add(r, c, coef * vK.get(i));
    }
}

template <typename SymmetricMatrix, typename Vector>
Vector DfCdkMatrix::getScreenedDensityMatrix(const RUN_TYPE runType,
                                             const PQ_PairArray& I2PR) {
    SymmetricMatrix P;
    switch (runType) {
        case RUN_RKS:
            P = 0.5 * this->getPInMatrix<SymmetricMatrix>(RUN_RKS, this->m_nIteration);
            break;

        case RUN_UKS_ALPHA:
        case RUN_UKS_BETA:
            P = this->getPInMatrix<SymmetricMatrix>(runType, this->m_nIteration);
            break;

        case RUN_ROKS_ALPHA: {
            P = 0.5 * this->getPInMatrix<SymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration);
            P += this->getPInMatrix<SymmetricMatrix>(RUN_ROKS_OPEN, this->m_nIteration);
        } break;

        case RUN_ROKS_BETA: {
            P = 0.5 * this->getPInMatrix<SymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration);
        } break;

        default:
            this->log_.critical(TlUtils::format("Program Error: %s:%d", __FILE__, __LINE__));
            CnErr.abort();
    }

    const std::size_t numOfI = I2PR.size();
    Vector answer(numOfI);

    for (std::size_t i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PR[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();
        answer.set(i, P.get(r, c));
    }

    return answer;
}

#endif  // DF_CDK_MATRIX_H
