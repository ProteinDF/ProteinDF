#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_impl_eigen.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#include "tl_sparse_symmetric_matrix_eigen.h"
#include "tl_sparse_symmetric_matrix_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_symmetric_matrix_impl_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlMatrixObject::index_type dim) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(dim);
}

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlDenseGeneralMatrix_Eigen& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<const TlDenseGeneralMatrix_ImplEigen*>(rhs.pImpl_)));
}

TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlSparseSymmetricMatrix_Eigen& sm) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<TlSparseSymmetricMatrix_ImplEigen*>(sm.pImpl_)));
}

#ifdef HAVE_VIENNACL
TlDenseSymmetricMatrix_Eigen::TlDenseSymmetricMatrix_Eigen(
    const TlDenseSymmetricMatrix_ViennaCL& rhs) {
  this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplViennaCL*>(rhs.pImpl_)));
}
#endif  // HAVE_VIENNACL

TlDenseSymmetricMatrix_Eigen::~TlDenseSymmetricMatrix_Eigen() {
  delete this->pImpl_;
  this->pImpl_ = NULL;
}

void TlDenseSymmetricMatrix_Eigen::vtr2mat(const std::vector<double>& vtr) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)->vtr2mat(vtr);
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  if (this != &rhs) {
    delete this->pImpl_;
    this->pImpl_ = new TlDenseSymmetricMatrix_ImplEigen(
        *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_)));
  }

  return *this;
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator+(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer += rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator-(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer -= rhs;
  return answer;
}

const TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::operator*(
    const TlDenseSymmetricMatrix_Eigen& rhs) const {
  TlDenseSymmetricMatrix_Eigen answer = *this;
  answer *= rhs;
  return answer;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator+=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) +=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));

  return *this;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator-=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) -=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));

  return *this;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator*=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) *= coef;

  return *this;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator/=(
    const double coef) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) /= coef;

  return *this;
}

TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::operator*=(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)) *=
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));

  return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
const TlDenseSymmetricMatrix_Eigen& TlDenseSymmetricMatrix_Eigen::dotInPlace(
    const TlDenseSymmetricMatrix_Eigen& rhs) {
  dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
      ->dotInPlace(
          *dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(rhs.pImpl_));

  return *this;
}

bool TlDenseSymmetricMatrix_Eigen::eig(
    TlDenseVector_Eigen* pEigVal, TlDenseGeneralMatrix_Eigen* pEigVec) const {
  TlDenseVector_ImplEigen* pImpl_eigval =
      dynamic_cast<TlDenseVector_ImplEigen*>(pEigVal->pImpl_);
  TlDenseGeneralMatrix_ImplEigen* pImpl_eigvec =
      dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(pEigVec->pImpl_);

  const bool answer =
      dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
          ->eig(pImpl_eigval, pImpl_eigvec);
  return answer;
}

TlDenseSymmetricMatrix_Eigen TlDenseSymmetricMatrix_Eigen::inverse() const {
  TlDenseSymmetricMatrix_ImplEigen tmp =
      dynamic_cast<const TlDenseSymmetricMatrix_ImplEigen*>(this->pImpl_)
          ->inverse();
  TlDenseSymmetricMatrix_Eigen answer(tmp);
  return answer;
}

// ---------------------------------------------------------------------------
// friend functions
// ---------------------------------------------------------------------------
TlDenseGeneralMatrix_Eigen operator*(const TlDenseGeneralMatrix_Eigen& mat1,
                                     const TlDenseSymmetricMatrix_Eigen& mat2) {
  const TlDenseGeneralMatrix_ImplEigen tmp =
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(mat1.pImpl_)) *
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(mat2.pImpl_));
  TlDenseGeneralMatrix_Eigen answer(tmp);
  return answer;
}

TlDenseGeneralMatrix_Eigen operator*(const TlDenseSymmetricMatrix_Eigen& mat1,
                                     const TlDenseGeneralMatrix_Eigen& mat2) {
  const TlDenseGeneralMatrix_ImplEigen tmp =
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(mat1.pImpl_)) *
      *(dynamic_cast<TlDenseGeneralMatrix_ImplEigen*>(mat2.pImpl_));
  TlDenseGeneralMatrix_Eigen answer(tmp);
  return answer;
}

TlDenseVector_Eigen operator*(const TlDenseSymmetricMatrix_Eigen& dms1,
                              const TlDenseVector_Eigen& dv) {
  return TlDenseVector_Eigen(
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(dms1.pImpl_)) *
      *(dynamic_cast<TlDenseVector_ImplEigen*>(dv.pImpl_)));
}

TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& dv,
                              const TlDenseSymmetricMatrix_Eigen& dms1) {
  return TlDenseVector_Eigen(
      *(dynamic_cast<TlDenseVector_ImplEigen*>(dv.pImpl_)) *
      *(dynamic_cast<TlDenseSymmetricMatrix_ImplEigen*>(dms1.pImpl_)));
}

TlDenseSymmetricMatrix_Eigen operator*(const double coef,
                                       const TlDenseSymmetricMatrix_Eigen& DM) {
  TlDenseSymmetricMatrix_Eigen answer = DM;
  answer *= coef;
  return answer;
}

TlDenseSymmetricMatrix_Eigen operator*(const TlDenseSymmetricMatrix_Eigen& DM,
                                       const double coef) {
  TlDenseSymmetricMatrix_Eigen answer = DM;
  answer *= coef;
  return answer;
}
