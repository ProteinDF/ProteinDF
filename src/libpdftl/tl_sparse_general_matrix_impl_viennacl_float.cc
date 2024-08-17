#include <map>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_EIGEN
#define VIENNACL_WITH_EIGEN 1
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>
#include <viennacl/tools/tools.hpp>

#include "tl_dense_general_matrix_impl_viennacl_float.h"
#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_matrix_utils.h"
#include "tl_sparse_general_matrix_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"
#include "viennacl/matrix.hpp"

TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col)
    : matrix_(row, col) {};

TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()) {
#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType tmp(rhs.getNumOfRows(), rhs.getNumOfCols());
        viennacl::copy(rhs.matrix_, tmp);
        viennacl::copy(tmp, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, float> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd = tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it = tmp[i].begin(); it != itEnd; ++it) {
                this->matrix_(i, it->first) = it->second;
            }
        }
    }
#endif  // HAVE_EIGEN
}

TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs) {
#ifdef HAVE_EIGEN
    {
        Eigen::SparseMatrix<float> tmp(rhs.getNumOfRows(), rhs.getNumOfCols());
        viennacl::copy(rhs.matrix_, tmp);
        Eigen::SparseMatrix<float> tmp2 = tmp.selfadjointView<Eigen::Lower>();
        viennacl::copy(tmp2, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, float> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) = it->second;
            }
        }
    }
#endif  // HAVE_EIGEN
}

TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const MatrixDataType& rhs) {
    this->matrix_ = rhs;
}

TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const TlDenseGeneralMatrix_ImplViennaCLFloat& DM_vcl) {
    TlDenseGeneralMatrix_ImplEigenFloat DM_eigen(DM_vcl);
    TlSparseGeneralMatrix_ImplEigenFloat SM_eigen(DM_eigen);
    viennacl::copy(SM_eigen.matrix_, this->matrix_);
}

#ifdef HAVE_EIGEN
TlSparseGeneralMatrix_ImplViennaCLFloat::TlSparseGeneralMatrix_ImplViennaCLFloat(
    const TlSparseGeneralMatrix_ImplEigenFloat& rhs) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_EIGEN

TlSparseGeneralMatrix_ImplViennaCLFloat::~TlSparseGeneralMatrix_ImplViennaCLFloat() {}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
TlMatrixObject::index_type TlSparseGeneralMatrix_ImplViennaCLFloat::getNumOfRows()
    const {
    return this->matrix_.size1();
}

TlMatrixObject::index_type TlSparseGeneralMatrix_ImplViennaCLFloat::getNumOfCols()
    const {
    return this->matrix_.size2();
}

void TlSparseGeneralMatrix_ImplViennaCLFloat::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
    const bool preserve = true;
    this->matrix_.resize(newRow, newCol, preserve);
}

double TlSparseGeneralMatrix_ImplViennaCLFloat::get(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type col) const {
    return static_cast<double>(this->matrix_(row, col));
}

void TlSparseGeneralMatrix_ImplViennaCLFloat::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                  const double value) {
    this->matrix_(row, col) = static_cast<float>(value);
}

void TlSparseGeneralMatrix_ImplViennaCLFloat::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                  const double value) {
    this->matrix_(row, col) += static_cast<float>(value);
}

void TlSparseGeneralMatrix_ImplViennaCLFloat::mul(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                  const double value) {
    this->matrix_(row, col) *= static_cast<float>(value);
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ImplViennaCLFloat& TlSparseGeneralMatrix_ImplViennaCLFloat::
operator=(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs) {
    // std::cout << "TlSparseGeneralMatrix_ImplViennaCLFloat& "
    //              "TlSparseGeneralMatrix_ImplViennaCLFloat::operator=(const "
    //              "TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)"
    //           << std::endl;
    if (this != &rhs) {
        this->matrix_.clear();
        this->matrix_.resize(rhs.getNumOfRows(), rhs.getNumOfCols());

        std::vector<std::map<unsigned int, float> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) = it->second;
            }
        }
    }

    return *this;
}

TlSparseGeneralMatrix_ImplViennaCLFloat& TlSparseGeneralMatrix_ImplViennaCLFloat::
operator+=(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs) {
#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType rhs1(this->getNumOfRows(), this->getNumOfCols());
        EigenMatrixDataType rhs2(rhs.getNumOfRows(), rhs.getNumOfCols());
        viennacl::copy(this->matrix_, rhs1);
        viennacl::copy(rhs.matrix_, rhs2);
        rhs1 += rhs2;
        viennacl::copy(rhs1, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, float> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) += it->second;
            }
        }
    }
#endif  // HAVE_EIGEN

    return *this;
}

TlSparseGeneralMatrix_ImplViennaCLFloat& TlSparseGeneralMatrix_ImplViennaCLFloat::
operator-=(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs) {
#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType rhs1(this->getNumOfRows(), this->getNumOfCols());
        EigenMatrixDataType rhs2(rhs.getNumOfRows(), rhs.getNumOfCols());
        viennacl::copy(this->matrix_, rhs1);
        viennacl::copy(rhs.matrix_, rhs2);
        rhs1 -= rhs2;
        viennacl::copy(rhs1, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, float> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) -= it->second;
            }
        }
    }
#endif  // HAVE_EIGEN

    return *this;
}

TlSparseGeneralMatrix_ImplViennaCLFloat& TlSparseGeneralMatrix_ImplViennaCLFloat::
operator*=(const float coef) {
#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType rhs1(this->getNumOfRows(), this->getNumOfCols());
        viennacl::copy(this->matrix_, rhs1);
        rhs1 *= coef;
        viennacl::copy(rhs1, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, float> > tmp(this->matrix_.size1());
        viennacl::copy(this->matrix_, tmp);
        const std::size_t size1 = this->matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, float>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, float>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) *= it->second;
            }
        }
    }
#endif  // HAVE_EIGEN
    return *this;
}

// ----------------------------------------------------------------------------
// I/O
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Others
// ----------------------------------------------------------------------------
// SM(G) = SM(G) * SM(G)
TlSparseGeneralMatrix_ImplViennaCLFloat operator*(const TlSparseGeneralMatrix_ImplViennaCLFloat& sm1,
                                                  const TlSparseGeneralMatrix_ImplViennaCLFloat& sm2) {
    assert(sm1.getNumOfCols() == sm2.getNumOfRows());
    TlSparseGeneralMatrix_ImplViennaCLFloat answer;
    answer.matrix_ = viennacl::linalg::prod(sm1.matrix_, sm2.matrix_);

    return answer;
}

// DV = SM(G) * DV
TlDenseVector_ImplViennaCLFloat operator*(const TlSparseGeneralMatrix_ImplViennaCLFloat& mat,
                                          const TlDenseVector_ImplViennaCLFloat& vtr) {
    // std::cout << "SM x DV (1)" << std::endl;
    assert(mat.getNumOfCols() == vtr.getSize());
    // std::cout << "SM x DV (2)" << std::endl;
    TlDenseVector_ImplViennaCLFloat answer;
    answer.vector_ = viennacl::linalg::prod(mat.matrix_, vtr.vector_);

    return answer;
}

// DV = DV * SM(G)
// TlDenseVector_ImplViennaCLFloat operator*(
//     const TlDenseVector_ImplViennaCLFloat& vtr,     const
//     TlSparseGeneralMatrix_ImplViennaCLFloat& mat) {
//   std::cout << "TlDenseVector_ImplViennaCLFloat operator*(const
//   TlDenseVector_ImplViennaCLFloat& vtr, const TlSparseGeneralMatrix_ImplViennaCLFloat&
//   mat)" << std::endl; assert(mat.getNumOfRows() == vtr.getSize());
//   TlDenseVector_ImplViennaCLFloat answer;
//   answer.vector_ = viennacl::linalg::prod(viennacl::trans(mat.matrix_),
//   vtr.vector_);

//   return answer;
// }
