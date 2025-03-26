#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <map>

#include "tl_dense_vector_impl_viennacl_float.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl_float.h"

#ifdef HAVE_EIGEN
#include "tl_dense_symmetric_matrix_impl_eigen_float.h"
#include "tl_sparse_symmetric_matrix_impl_eigen_float.h"
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ImplViennaCLFloat::TlSparseSymmetricMatrix_ImplViennaCLFloat(
    const TlMatrixObject::index_type dim)
    : TlSparseGeneralMatrix_ImplViennaCLFloat(dim, dim) {}

TlSparseSymmetricMatrix_ImplViennaCLFloat::TlSparseSymmetricMatrix_ImplViennaCLFloat(
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs)
    : TlSparseGeneralMatrix_ImplViennaCLFloat(rhs.getNumOfRows(),
                                              rhs.getNumOfCols()) {
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

TlSparseSymmetricMatrix_ImplViennaCLFloat::TlSparseSymmetricMatrix_ImplViennaCLFloat(
    const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs)
    : TlSparseGeneralMatrix_ImplViennaCLFloat(rhs.getNumOfRows(),
                                              rhs.getNumOfCols()) {
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

TlSparseSymmetricMatrix_ImplViennaCLFloat::TlSparseSymmetricMatrix_ImplViennaCLFloat(
    const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs) {
    TlSparseSymmetricMatrix_ImplEigenFloat eigenSM = TlDenseSymmetricMatrix_ImplEigenFloat(rhs);
    viennacl::copy(eigenSM.matrix_, this->matrix_);
}

#ifdef HAVE_EIGEN
TlSparseSymmetricMatrix_ImplViennaCLFloat::TlSparseSymmetricMatrix_ImplViennaCLFloat(
    const TlSparseSymmetricMatrix_ImplEigenFloat& rhs) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ImplViennaCLFloat::~TlSparseSymmetricMatrix_ImplViennaCLFloat() {}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
double TlSparseSymmetricMatrix_ImplViennaCLFloat::get(const TlMatrixObject::index_type row,
                                                      const TlMatrixObject::index_type col) const {
    float answer;
    // access elements in top-right angular matrix
    if (row >= col) {
        answer = this->matrix_(row, col);
    } else {
        answer = this->matrix_(col, row);
    }

    return static_cast<double>(answer);
}

void TlSparseSymmetricMatrix_ImplViennaCLFloat::set(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                    const double value) {
    const float v = static_cast<float>(value);
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) = v;
    } else {
        this->matrix_(col, row) = v;
    }
}

void TlSparseSymmetricMatrix_ImplViennaCLFloat::add(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                    const double value) {
    const float v = static_cast<float>(value);
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) += v;
    } else {
        this->matrix_(col, row) += v;
    }
}

void TlSparseSymmetricMatrix_ImplViennaCLFloat::mul(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                                                    const double value) {
    const float v = static_cast<float>(value);
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) *= v;
    } else {
        this->matrix_(col, row) *= v;
    }
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseSymmetricMatrix_ImplViennaCLFloat& TlSparseSymmetricMatrix_ImplViennaCLFloat::
operator=(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs) {
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

TlSparseSymmetricMatrix_ImplViennaCLFloat& TlSparseSymmetricMatrix_ImplViennaCLFloat::
operator+=(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs) {
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

TlSparseSymmetricMatrix_ImplViennaCLFloat& TlSparseSymmetricMatrix_ImplViennaCLFloat::
operator-=(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs) {
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

TlSparseSymmetricMatrix_ImplViennaCLFloat& TlSparseSymmetricMatrix_ImplViennaCLFloat::
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
// protected
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// others
// ----------------------------------------------------------------------------
// SM(G) = SM(G) * SM(S)
TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
    const TlSparseGeneralMatrix_ImplViennaCLFloat& sm1,
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm2) {
    return sm1 * TlSparseGeneralMatrix_ImplViennaCLFloat(sm2);
}
// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm1,
    const TlSparseGeneralMatrix_ImplViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ImplViennaCLFloat(sm1) * sm2;
}
// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm1,
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm2) {
    return TlSparseGeneralMatrix_ImplViennaCLFloat(sm1) *
           TlSparseGeneralMatrix_ImplViennaCLFloat(sm2);
}

TlDenseVector_ImplViennaCLFloat operator*(
    const TlSparseSymmetricMatrix_ImplViennaCLFloat& mat,
    const TlDenseVector_ImplViennaCLFloat& vtr) {
    assert(mat.getNumOfCols() == vtr.getSize());
    return (TlSparseGeneralMatrix_ImplViennaCLFloat(mat) * vtr);
}
