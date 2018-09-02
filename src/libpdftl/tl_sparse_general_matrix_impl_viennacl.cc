#include <map>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_EIGEN
#define VIENNACL_WITH_EIGEN 1
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "tl_dense_general_matrix_impl_eigen.h"
#include "tl_sparse_general_matrix_impl_eigen.h"
#endif  // HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>
#include <viennacl/tools/tools.hpp>
#include "viennacl/matrix.hpp"

#include "tl_dense_general_matrix_impl_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_matrix_utils.h"
#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_symmetric_matrix_impl_viennacl.h"

TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col)
    : matrix_(row, col){};

TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const TlSparseGeneralMatrix_ImplViennaCL& rhs)
    : matrix_(rhs.getNumOfRows(), rhs.getNumOfCols()){
      std::cout << "TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(const TlSparseGeneralMatrix_ImplViennaCL& rhs)" << std::endl;
#ifdef HAVE_EIGEN
          {
            EigenMatrixDataType tmp(rhs.getNumOfRows(), rhs.getNumOfCols());
  viennacl::copy(rhs.matrix_, tmp);
  viennacl::copy(tmp, this->matrix_);
}
#else
          {std::vector<std::map<unsigned int, double> > tmp(
              rhs.matrix_.size1());
viennacl::copy(rhs.matrix_, tmp);
const std::size_t size1 = rhs.matrix_.size1();
for (std::size_t i = 0; i < size1; ++i) {
  typename std::map<unsigned int, double>::const_iterator itEnd = tmp[i].end();
  for (typename std::map<unsigned int, double>::const_iterator it =
           tmp[i].begin();
       it != itEnd; ++it) {
    this->matrix_(i, it->first) = it->second;
  }
}
}
#endif  // HAVE_EIGEN
}

TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const TlSparseSymmetricMatrix_ImplViennaCL& rhs){
#ifdef HAVE_EIGEN
    {Eigen::SparseMatrix<double> tmp(rhs.getNumOfRows(), rhs.getNumOfCols());
viennacl::copy(rhs.matrix_, tmp);
Eigen::SparseMatrix<double> tmp2 = tmp.selfadjointView<Eigen::Lower>();
viennacl::copy(tmp2, this->matrix_);
}
#else
    {std::vector<std::map<unsigned int, double> > tmp(rhs.matrix_.size1());
viennacl::copy(rhs.matrix_, tmp);
const std::size_t size1 = rhs.matrix_.size1();
for (std::size_t i = 0; i < size1; ++i) {
  typename std::map<unsigned int, double>::const_iterator itEnd = tmp[i].end();
  for (typename std::map<unsigned int, double>::const_iterator it =
           tmp[i].begin();
       it != itEnd; ++it) {
    this->matrix_(i, it->first) = it->second;
  }
}
}
#endif  // HAVE_EIGEN
}

TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const MatrixDataType& rhs) {
  this->matrix_ = rhs;
}

TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const TlDenseGeneralMatrix_ImplViennaCL& DM_vcl) {
  TlDenseGeneralMatrix_ImplEigen DM_eigen(DM_vcl);
  TlSparseGeneralMatrix_ImplEigen SM_eigen(DM_eigen);
  viennacl::copy(SM_eigen.matrix_, this->matrix_);
}

#ifdef HAVE_EIGEN
TlSparseGeneralMatrix_ImplViennaCL::TlSparseGeneralMatrix_ImplViennaCL(
    const TlSparseGeneralMatrix_ImplEigen& rhs) {
  viennacl::copy(rhs.matrix_, this->matrix_);
}
#endif  // HAVE_EIGEN

TlSparseGeneralMatrix_ImplViennaCL::~TlSparseGeneralMatrix_ImplViennaCL() {}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
TlMatrixObject::index_type TlSparseGeneralMatrix_ImplViennaCL::getNumOfRows()
    const {
  return this->matrix_.size1();
}

TlMatrixObject::index_type TlSparseGeneralMatrix_ImplViennaCL::getNumOfCols()
    const {
  return this->matrix_.size2();
}

void TlSparseGeneralMatrix_ImplViennaCL::resize(
    const TlMatrixObject::index_type newRow,
    const TlMatrixObject::index_type newCol) {
  const bool preserve = true;
  this->matrix_.resize(newRow, newCol, preserve);
}

double TlSparseGeneralMatrix_ImplViennaCL::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->matrix_(row, col);
}

void TlSparseGeneralMatrix_ImplViennaCL::set(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  this->matrix_(row, col) = value;
}

void TlSparseGeneralMatrix_ImplViennaCL::add(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  this->matrix_(row, col) += value;
}

void TlSparseGeneralMatrix_ImplViennaCL::mul(
    const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
    const double value) {
  this->matrix_(row, col) *= value;
}

// ----------------------------------------------------------------------------
// operators
// ----------------------------------------------------------------------------
TlSparseGeneralMatrix_ImplViennaCL& TlSparseGeneralMatrix_ImplViennaCL::
operator=(const TlSparseGeneralMatrix_ImplViennaCL& rhs) {
  std::cout << "TlSparseGeneralMatrix_ImplViennaCL& TlSparseGeneralMatrix_ImplViennaCL::operator=(const TlSparseGeneralMatrix_ImplViennaCL& rhs)" << std::endl;
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

TlSparseGeneralMatrix_ImplViennaCL& TlSparseGeneralMatrix_ImplViennaCL::
operator+=(const TlSparseGeneralMatrix_ImplViennaCL& rhs) {
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

TlSparseGeneralMatrix_ImplViennaCL& TlSparseGeneralMatrix_ImplViennaCL::
operator-=(const TlSparseGeneralMatrix_ImplViennaCL& rhs) {
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

TlSparseGeneralMatrix_ImplViennaCL& TlSparseGeneralMatrix_ImplViennaCL::
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
// I/O
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Others
// ----------------------------------------------------------------------------
// SM(G) = SM(G) * SM(G)
TlSparseGeneralMatrix_ImplViennaCL operator*(
    const TlSparseGeneralMatrix_ImplViennaCL& sm1,
    const TlSparseGeneralMatrix_ImplViennaCL& sm2) {
  assert(sm1.getNumOfCols() == sm2.getNumOfRows());
  TlSparseGeneralMatrix_ImplViennaCL answer;
  answer.matrix_ = viennacl::linalg::prod(sm1.matrix_, sm2.matrix_);

  return answer;
}

// DV = SM(G) * DV
TlDenseVector_ImplViennaCL operator*(
    const TlSparseGeneralMatrix_ImplViennaCL& mat,
    const TlDenseVector_ImplViennaCL& vtr) {
  //std::cout << "SM x DV (1)" << std::endl;
  assert(mat.getNumOfCols() == vtr.getSize());
  //std::cout << "SM x DV (2)" << std::endl;
  TlDenseVector_ImplViennaCL answer;
  answer.vector_ = viennacl::linalg::prod(mat.matrix_, vtr.vector_);

  return answer;
}

// DV = DV * SM(G)
// TlDenseVector_ImplViennaCL operator*(
//     const TlDenseVector_ImplViennaCL& vtr,     const TlSparseGeneralMatrix_ImplViennaCL& mat) {
//   std::cout << "TlDenseVector_ImplViennaCL operator*(const TlDenseVector_ImplViennaCL& vtr, const TlSparseGeneralMatrix_ImplViennaCL& mat)" << std::endl;
//   assert(mat.getNumOfRows() == vtr.getSize());
//   TlDenseVector_ImplViennaCL answer;
//   answer.vector_ = viennacl::linalg::prod(viennacl::trans(mat.matrix_), vtr.vector_);

//   return answer;
// }
