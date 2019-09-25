#include <map>

#include "tl_dense_vector_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ImplViennaCL::TlSparseSymmetricMatrix_ImplViennaCL(
    const TlMatrixObject::index_type dim)
    : TlSparseGeneralMatrix_ImplViennaCL(dim, dim) {}

TlSparseSymmetricMatrix_ImplViennaCL::TlSparseSymmetricMatrix_ImplViennaCL(
    const TlSparseSymmetricMatrix_ImplViennaCL& rhs)
    : TlSparseGeneralMatrix_ImplViennaCL(rhs.getNumOfRows(),
                                         rhs.getNumOfCols()) {
    std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
    viennacl::copy(rhs.matrix_, tmp);
    const std::size_t size1 = rhs.matrix_.size1();
    for (std::size_t i = 0; i < size1; ++i) {
        typename std::map<unsigned int, double>::const_iterator itEnd =
            tmp[i].end();
        for (typename std::map<unsigned int, double>::const_iterator it =
                 tmp[i].begin();
             it != itEnd; ++it) {
            this->matrix_(i, it->first) = it->second;
        }
    }
}

TlSparseSymmetricMatrix_ImplViennaCL::TlSparseSymmetricMatrix_ImplViennaCL(
    const TlSparseGeneralMatrix_ImplViennaCL& rhs)
    : TlSparseGeneralMatrix_ImplViennaCL(rhs.getNumOfRows(),
                                         rhs.getNumOfCols()) {
    std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
    viennacl::copy(rhs.matrix_, tmp);
    const std::size_t size1 = rhs.matrix_.size1();
    for (std::size_t i = 0; i < size1; ++i) {
        typename std::map<unsigned int, double>::const_iterator itEnd =
            tmp[i].end();
        for (typename std::map<unsigned int, double>::const_iterator it =
                 tmp[i].begin();
             it != itEnd; ++it) {
            this->matrix_(i, it->first) = it->second;
        }
    }
}

TlSparseSymmetricMatrix_ImplViennaCL::TlSparseSymmetricMatrix_ImplViennaCL(
    const TlDenseSymmetricMatrix_ImplViennaCL& rhs) {
    TlSparseSymmetricMatrix_ImplEigen eigenSM =
        TlDenseSymmetricMatrix_ImplEigen(TlDenseSymmetricMatrix_ImplEigen(rhs));
    viennacl::copy(eigenSM.matrix_, this->matrix_);
}

#ifdef HAVE_EIGEN
TlSparseSymmetricMatrix_ImplViennaCL::TlSparseSymmetricMatrix_ImplViennaCL(
    const TlSparseSymmetricMatrix_ImplEigen& rhs) {
    viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_EIGEN

TlSparseSymmetricMatrix_ImplViennaCL::~TlSparseSymmetricMatrix_ImplViennaCL() {}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
double TlSparseSymmetricMatrix_ImplViennaCL::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
    double answer;
    // access elements in top-right angular matrix
    if (row >= col) {
        answer = this->matrix_(row, col);
    } else {
        answer = this->matrix_(col, row);
    }

    return answer;
}

void TlSparseSymmetricMatrix_ImplViennaCL::set(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) = value;
    } else {
        this->matrix_(col, row) = value;
    }
}

void TlSparseSymmetricMatrix_ImplViennaCL::add(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) += value;
    } else {
        this->matrix_(col, row) += value;
    }
}

void TlSparseSymmetricMatrix_ImplViennaCL::mul(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
    // set elements in top-right angular matrix
    if (row >= col) {
        this->matrix_(row, col) *= value;
    } else {
        this->matrix_(col, row) *= value;
    }
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseSymmetricMatrix_ImplViennaCL& TlSparseSymmetricMatrix_ImplViennaCL::
operator=(const TlSparseSymmetricMatrix_ImplViennaCL& rhs) {
    if (this != &rhs) {
        this->matrix_.clear();
        this->matrix_.resize(rhs.getNumOfRows(), rhs.getNumOfCols());

        std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, double>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, double>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) = it->second;
            }
        }
    }

    return *this;
}

TlSparseSymmetricMatrix_ImplViennaCL& TlSparseSymmetricMatrix_ImplViennaCL::
operator+=(const TlSparseSymmetricMatrix_ImplViennaCL& rhs) {
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
        std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, double>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, double>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) += it->second;
            }
        }
    }
#endif  // HAVE_EIGEN

    return *this;
}

TlSparseSymmetricMatrix_ImplViennaCL& TlSparseSymmetricMatrix_ImplViennaCL::
operator-=(const TlSparseSymmetricMatrix_ImplViennaCL& rhs) {
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
        std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
        viennacl::copy(rhs.matrix_, tmp);
        const std::size_t size1 = rhs.matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, double>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, double>::const_iterator it =
                     tmp[i].begin();
                 it != itEnd; ++it) {
                this->matrix_(i, it->first) -= it->second;
            }
        }
    }
#endif  // HAVE_EIGEN

    return *this;
}

TlSparseSymmetricMatrix_ImplViennaCL& TlSparseSymmetricMatrix_ImplViennaCL::
operator*=(const double coef) {
#ifdef HAVE_EIGEN
    {
        EigenMatrixDataType rhs1(this->getNumOfRows(), this->getNumOfCols());
        viennacl::copy(this->matrix_, rhs1);
        rhs1 *= coef;
        viennacl::copy(rhs1, this->matrix_);
    }
#else
    {
        std::vector<std::map<unsigned int, double> > tmp(this->matrix_.size1());
        viennacl::copy(this->matrix_, tmp);
        const std::size_t size1 = this->matrix_.size1();
        for (std::size_t i = 0; i < size1; ++i) {
            typename std::map<unsigned int, double>::const_iterator itEnd =
                tmp[i].end();
            for (typename std::map<unsigned int, double>::const_iterator it =
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
TlSparseGeneralMatrix_ImplViennaCL operator*(
    const TlSparseGeneralMatrix_ImplViennaCL& sm1,
    const TlSparseSymmetricMatrix_ImplViennaCL& sm2) {
    return sm1 * TlSparseGeneralMatrix_ImplViennaCL(sm2);
}
// SM(G) = SM(S) * SM(G)
TlSparseGeneralMatrix_ImplViennaCL operator*(
    const TlSparseSymmetricMatrix_ImplViennaCL& sm1,
    const TlSparseGeneralMatrix_ImplViennaCL& sm2) {
    return TlSparseGeneralMatrix_ImplViennaCL(sm1) * sm2;
}
// SM(G) = SM(S) * SM(S)
TlSparseGeneralMatrix_ImplViennaCL operator*(
    const TlSparseSymmetricMatrix_ImplViennaCL& sm1,
    const TlSparseSymmetricMatrix_ImplViennaCL& sm2) {
    return TlSparseGeneralMatrix_ImplViennaCL(sm1) *
           TlSparseGeneralMatrix_ImplViennaCL(sm2);
}

TlDenseVector_ImplViennaCL operator*(
    const TlSparseSymmetricMatrix_ImplViennaCL& mat,
    const TlDenseVector_ImplViennaCL& vtr) {
    assert(mat.getNumOfCols() == vtr.getSize());
    return (TlSparseGeneralMatrix_ImplViennaCL(mat) * vtr);
}
