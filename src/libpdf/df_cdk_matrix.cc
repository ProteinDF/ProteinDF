#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include "CnError.h"
#include "DfTaskCtrl.h"
#include "df_cdk_matrix.h"

#include "tl_dense_general_matrix_mmap.h"

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen_float.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_sparse_general_matrix_eigen.h"
#include "tl_sparse_general_matrix_eigen_float.h"
#include "tl_sparse_symmetric_matrix_eigen.h"
#include "tl_sparse_symmetric_matrix_eigen_float.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_general_matrix_viennacl_float.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl_float.h"
#include "tl_dense_vector_viennacl.h"
#include "tl_dense_vector_viennacl_float.h"
#include "tl_sparse_general_matrix_viennacl.h"
#include "tl_sparse_general_matrix_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_viennacl.h"
#include "tl_sparse_symmetric_matrix_viennacl_float.h"
#endif  // HAVE_VIENNACL

// ----------------------------------------------------------------------------
// construct & destruct
// ----------------------------------------------------------------------------
DfCdkMatrix::DfCdkMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
    this->updateLinearAlgebraPackageParam((*(this->pPdfParam_))["linear_algebra_package/K"].getStr());

    this->useFP32_ = false;
    if (!(*this->pPdfParam_)["CD/use_fp32"].getStr().empty()) {
        this->useFP32_ = (*this->pPdfParam_)["CD/use_fp32"].getBoolean();
    }

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
            this->getK_runType<TlDenseSymmetricMatrix_Lapack,
                               TlDenseVector_Lapack,
                               TlDenseGeneralMatrix_Lapack,
                               TlDenseSymmetricMatrix_Lapack,
                               TlDenseSymmetricMatrix_Lapack>();
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            if (this->useFP32_ == true) {
                this->log_.info("Linear Algebra Package: Eigen(temporary: FP32)");
                this->getK_runType<TlDenseSymmetricMatrix_Eigen,
                                   TlDenseVector_Eigen,
                                   TlSparseGeneralMatrix_Eigen,
                                   TlSparseSymmetricMatrix_Eigen,
                                   TlDenseSymmetricMatrix_EigenFloat>();
            } else {
                this->log_.info("Linear Algebra Package: Eigen");
                this->getK_runType<TlDenseSymmetricMatrix_Eigen,
                                   TlDenseVector_Eigen,
                                   TlSparseGeneralMatrix_Eigen,
                                   TlDenseSymmetricMatrix_Eigen,
                                   TlDenseSymmetricMatrix_Eigen>();
            }
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            if (this->useFP32_ == true) {
                this->log_.info("Linear Algebra Package: ViennaCL(temporary: FP32)");
                this->getK_runType<TlDenseSymmetricMatrix_ViennaCL,
                                   TlDenseVector_ViennaCL,
                                   TlSparseGeneralMatrix_ViennaCL,
                                   TlDenseSymmetricMatrix_ViennaCL,
                                   TlDenseSymmetricMatrix_ViennaCLFloat>();
            } else {
                this->log_.info("Linear Algebra Package: ViennaCL");
                this->getK_runType<TlDenseSymmetricMatrix_ViennaCL,
                                   TlDenseVector_ViennaCL,
                                   TlSparseGeneralMatrix_ViennaCL,
                                   TlDenseSymmetricMatrix_ViennaCL,
                                   TlDenseSymmetricMatrix_ViennaCL>();
            }
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

// -----------------------------------------------------------------------------
// task control
// -----------------------------------------------------------------------------
DfTaskCtrl* DfCdkMatrix::getDfTaskCtrlObject() const {
    this->log_.info("task control: serial");
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
