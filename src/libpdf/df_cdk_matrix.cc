#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include "CnError.h"
#include "DfTaskCtrl.h"
#include "df_cdk_matrix.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"
#include "tl_dense_general_matrix_mmap.h"

// ----------------------------------------------------------------------------
// construct & destruct
// ----------------------------------------------------------------------------
DfCdkMatrix::DfCdkMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
    this->updateLinearAlgebraPackageParam(
        (*(this->pPdfParam_))["linear_algebra_package/K"].getStr());

    this->useMmapMatrix_ = true;
    if (!(*this->pPdfParam_)["CD/use_mmap_matrix"].getStr().empty()) {
        this->useMmapMatrix_ =
            (*this->pPdfParam_)["CD/use_mmap_matrix"].getBoolean();
    }
    this->isCvSavedAsMmap_ = true;
    if (!(*this->pPdfParam_)["CD/is_cv_saved_as_mmap"].getStr().empty()) {
        this->isCvSavedAsMmap_ =
            (*this->pPdfParam_)["CD/is_cv_saved_as_mmap"].getBoolean();
    }
    if (this->useMmapMatrix_ == true) {
        this->isCvSavedAsMmap_ = true;
    }

    this->sparseMatrixLevel_ = 0;
    if (!(*this->pPdfParam_)["CD/sparse_matrix_level"].getStr().empty()) {
        this->sparseMatrixLevel_ =
            (*this->pPdfParam_)["CD/sparse_matrix_level"].getInt();
    }

    if (this->K_engine_ == DfObject::K_ENGINE_CD) {
        this->fastCDK_mode_ = FASTCDK_NONE;
    } else {
        assert(this->K_engine_ == DfObject::K_ENGINE_FASTCDK);
        this->fastCDK_mode_ = FASTCDK_PRODUCTIVE;

        if ((*pPdfParam)["debug/DfCD/FastCDK_mode"].getStr().empty() != true) {
            const std::string fastCDK_mode_str = TlUtils::toUpper(
                (*pPdfParam)["debug/DfCD/FastCDK_mode"].getStr());
            if (fastCDK_mode_str == "PRODUCTIVE_FULL") {
                this->fastCDK_mode_ = FASTCDK_PRODUCTIVE_FULL;
            } else if (fastCDK_mode_str == "FULL_SUPERMATRIX") {
                this->fastCDK_mode_ = FASTCDK_DEBUG_FULL_SUPERMATRIX;
            } else if (fastCDK_mode_str == "SUPERMATRIX") {
                this->fastCDK_mode_ = FASTCDK_DEBUG_SUPERMATRIX;
            }
        }
    }
}

DfCdkMatrix::~DfCdkMatrix() {}

// ----------------------------------------------------------------------------
// public
// ----------------------------------------------------------------------------
void DfCdkMatrix::getK() {
    switch (this->linearAlgebraPackage_) {
        case LAP_LAPACK: {
            this->log_.info("Linear Algebra Package: Lapack");
            this->getK_method<TlDenseSymmetricMatrix_Lapack,
                              TlDenseVector_Lapack, TlDenseGeneralMatrix_Lapack,
                              TlDenseSymmetricMatrix_Lapack>();
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            this->log_.info("Linear Algebra Package: Eigen");
            this->getK_method<TlDenseSymmetricMatrix_Eigen, TlDenseVector_Eigen,
                              TlSparseGeneralMatrix_Eigen,
                              TlSparseSymmetricMatrix_Eigen>();
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            this->log_.info("Linear Algebra Package: ViennaCL");
            this->getK_method<TlDenseSymmetricMatrix_ViennaCL,
                              TlDenseVector_ViennaCL,
                              TlSparseGeneralMatrix_ViennaCL,
                              TlSparseSymmetricMatrix_ViennaCL>();
        } break;
#endif  // HAVE_VIENNACL

        default:
            CnErr.abort(
                TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
    }
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
template <typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
void DfCdkMatrix::getK_method() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            this->getK_L<SymmetricMatrix, Vector, SparseGeneralMatrix,
                         SparseSymmetricMatrix>(RUN_RKS);
        } break;

        case METHOD_UKS: {
            this->getK_L<SymmetricMatrix, Vector, SparseGeneralMatrix,
                         SparseSymmetricMatrix>(RUN_UKS_ALPHA);
            this->getK_L<SymmetricMatrix, Vector, SparseGeneralMatrix,
                         SparseSymmetricMatrix>(RUN_UKS_BETA);

        } break;

        case METHOD_ROKS: {
            this->getK_L<SymmetricMatrix, Vector, SparseGeneralMatrix,
                         SparseSymmetricMatrix>(RUN_ROKS_ALPHA);
            this->getK_L<SymmetricMatrix, Vector, SparseGeneralMatrix,
                         SparseSymmetricMatrix>(RUN_ROKS_BETA);
        } break;

        default:
            this->log_.critical("program error.");
            CnErr.abort();
            break;
    }
}

template <typename SymmetricMatrix, typename Vector,
          typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
void DfCdkMatrix::getK_L(DfObject::RUN_TYPE runType) {
    switch (this->fastCDK_mode_) {
        case FASTCDK_NONE:
            if (this->isCvSavedAsMmap_) {
                this->log_.info("loading the cholesky vectors on mmap.");
                switch (this->sparseMatrixLevel_) {
                    case 1: {
                        this->log_.info(
                            "use sparse matrix to transpose the cholesky "
                            "vector.");
                        this->getK_byLjk_useTransMatrix<
                            TlDenseGeneralMatrix_mmap, SymmetricMatrix, Vector,
                            SparseGeneralMatrix>(runType);
                    } break;

                    case 2: {
                        this->log_.info(
                            "use sparse matrix to transpose the cholesky "
                            "vector.");
                        this->log_.info(
                            "use sparse matrix in temporary K matrix.");
                        this->getK_byLjk_useSparseMatrix<
                            TlDenseGeneralMatrix_mmap, SymmetricMatrix, Vector,
                            SparseGeneralMatrix, SparseSymmetricMatrix>(
                            runType);

                    } break;

                    default: {
                        this->log_.info("use dense matrix for building K.");
                        this->getK_byLjk_useDenseMatrix<
                            TlDenseGeneralMatrix_mmap, SymmetricMatrix, Vector>(
                            runType);
                    } break;
                }

            } else {
                this->log_.info("loading the cholesky vectors on arrays.");
                // v1
                this->getK_byLjk_useDenseMatrix<
                    TlDenseGeneralMatrix_arrays_ColOriented, SymmetricMatrix,
                    Vector>(runType);
            }
            break;

        case FASTCDK_DEBUG_FULL_SUPERMATRIX:
        case FASTCDK_DEBUG_SUPERMATRIX:
        case FASTCDK_PRODUCTIVE_FULL:
        case FASTCDK_PRODUCTIVE:
            this->getK_byLk<TlDenseGeneralMatrix_mmap, SymmetricMatrix, Vector>(
                runType);
            break;

        default:
            this->log_.critical(
                TlUtils::format("%s: %d: program error.", __FILE__, __LINE__));
            CnErr.abort();
    }
}

template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector>
void DfCdkMatrix::getK_byLjk_useDenseMatrix(const RUN_TYPE runType) {
    this->log_.info("calc K by CD method (serial).");
    SymmetricMatrix K(this->m_nNumOfAOs);

    const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
    this->log_.info(
        TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    SymmetricMatrix P = this->getPInMatrix<SymmetricMatrix>(runType, this->m_nIteration);
    if (runType == RUN_RKS) {
        P *= 0.5;
    }

    this->log_.info("start loop");
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

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
            const SymmetricMatrix l =
                this->getCholeskyVector<SymmetricMatrix, Vector>(
                    L.getColVector(tasks[i]), I2PQ);
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
          typename SparseGeneralMatrix>
void DfCdkMatrix::getK_byLjk_useTransMatrix(const RUN_TYPE runType) {
    this->log_.info("calc K by CD method (serial).");
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
            const SymmetricMatrix l =
                this->convert_I2PQ<SymmetricMatrix, Vector,
                                   SparseGeneralMatrix>(
                    I2PQ_mat, L.getColVector(tasks[i]));
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
    this->log_.info("calc K by CD method (serial).");
    SymmetricMatrix K(this->m_nNumOfAOs);

    const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
    this->log_.info(
        TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
    const index_type numOfCBs = L.getNumOfCols();

    const SymmetricMatrix P = this->getSpinDensityMatrix<SymmetricMatrix>(
        runType, this->m_nIteration - 1);
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

// -----------------------------------------------------------------------------
// task control
// -----------------------------------------------------------------------------
DfTaskCtrl* DfCdkMatrix::getDfTaskCtrlObject() const {
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);

    return pDfTaskCtrl;
}

// -----------------------------------------------------------------------------
// I/O
// -----------------------------------------------------------------------------
PQ_PairArray DfCdkMatrix::getI2PQ(const std::string& filepath) {
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        this->log_.critical(
            TlUtils::format("could not found: %s", filepath.c_str()));
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = Index2(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}

template <typename SymmetricMatrix>
void DfCdkMatrix::finalize(SymmetricMatrix* pK) {
    // do nothing on serial code
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

// special version -------------------------------------------------------------
#ifdef HAVE_VIENNACL
template <>
TlDenseSymmetricMatrix_ViennaCL DfCdkMatrix::getCholeskyVector<
    TlDenseSymmetricMatrix_ViennaCL, TlDenseVector_ViennaCL>(
    const TlDenseVector_ViennaCL& L_col, const PQ_PairArray& I2PQ) {
    const index_type numOfItilde = L_col.getSize();
    assert(static_cast<std::size_t>(numOfItilde) == I2PQ.size());

    // this->log_.info("convert matrix is build via Eigen.");
    const TlDenseVector_Eigen L_col_eigen = L_col;
    TlDenseSymmetricMatrix_Eigen answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(), I2PQ[i].index2(), L_col_eigen.get(i));
    }

    return TlDenseSymmetricMatrix_ViennaCL(answer);
}
#endif  // HAVE_VIENNACL

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

// special version -------------------------------------------------------------
#ifdef HAVE_VIENNACL
template <>
TlSparseGeneralMatrix_ViennaCL
DfCdkMatrix::getTrans_I2PQ_Matrix<TlSparseGeneralMatrix_ViennaCL>(
    const PQ_PairArray& I2PQ) {
    this->log_.info("trans matrix is built via Eigen.");
    const TlMatrixObject::index_type numOfItilde = I2PQ.size();
    const TlMatrixObject::index_type numOfAoPairs =
        this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    TlSparseGeneralMatrix_Eigen expandL(numOfAoPairs, numOfItilde);
    for (TlMatrixObject::index_type i = 0; i < numOfItilde; ++i) {
        TlMatrixObject::index_type row = I2PQ[i].index1();
        TlMatrixObject::index_type col = I2PQ[i].index2();
        if (row < col) {
            std::swap(row, col);
        }
        return TlSparseGeneralMatrix_ViennaCL(expandL);
    }
}
#endif  // HAVE_VIENNAACL

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
