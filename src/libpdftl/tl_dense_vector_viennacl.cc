#include "tl_dense_vector_viennacl.h"
#include "tl_dense_vector_impl_viennacl.h"

#ifdef HAVE_EIGEN
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_impl_eigen.h"
#endif  // HAVE_EIGEN

TlDenseVector_ViennaCL::TlDenseVector_ViennaCL(
    const TlDenseVectorObject::index_type size) {
    this->pImpl_ = new TlDenseVector_ImplViennaCL(size);
}

TlDenseVector_ViennaCL::TlDenseVector_ViennaCL(
    const TlDenseVector_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCL(
        *dynamic_cast<const TlDenseVector_ImplViennaCL*>(rhs.pImpl_));
}

TlDenseVector_ViennaCL::TlDenseVector_ViennaCL(const std::vector<double>& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCL(rhs);
}

#ifdef HAVE_EIGEN
TlDenseVector_ViennaCL::TlDenseVector_ViennaCL(const TlDenseVector_Eigen& rhs) {
    this->pImpl_ = new TlDenseVector_ImplViennaCL(
        *(dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_)));
}
#endif  // HAVE_EIGEN

TlDenseVector_ViennaCL::operator std::vector<double>() const {
    return std::vector<double>(
        *(dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)));
}

TlDenseVector_ViennaCL::~TlDenseVector_ViennaCL() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::operator=(
    const TlDenseVector_ViennaCL& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseVector_ImplViennaCL(
            *dynamic_cast<const TlDenseVector_ImplViennaCL*>(rhs.pImpl_));
    }

    return *this;
}

TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::operator+=(
    const TlDenseVector_ViennaCL& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)
        ->
        operator+=(*dynamic_cast<TlDenseVector_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::operator-=(
    const TlDenseVector_ViennaCL& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)
        ->
        operator-=(*dynamic_cast<TlDenseVector_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::operator*=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)->operator*=(rhs);

    return *this;
}

TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::operator/=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)->operator/=(rhs);

    return *this;
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::dotInPlace(
    const TlDenseVector_ViennaCL& rhs) {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)
        ->dotInPlace(*dynamic_cast<TlDenseVector_ImplViennaCL*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_ViennaCL& TlDenseVector_ViennaCL::reverse() {
    dynamic_cast<TlDenseVector_ImplViennaCL*>(this->pImpl_)->reverse();

    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_ViennaCL operator+(const TlDenseVector_ViennaCL& rhs1,
                                 const TlDenseVector_ViennaCL& rhs2) {
    TlDenseVector_ViennaCL answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_ViennaCL operator-(const TlDenseVector_ViennaCL& rhs1,
                                 const TlDenseVector_ViennaCL& rhs2) {
    TlDenseVector_ViennaCL answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_ViennaCL operator*(const TlDenseVector_ViennaCL& rhs1,
                                 const double rhs2) {
    TlDenseVector_ViennaCL answer = rhs1;
    answer *= rhs2;

    return answer;
}

TlDenseVector_ViennaCL operator*(const double rhs1,
                                 const TlDenseVector_ViennaCL& rhs2) {
    return rhs2 * rhs1;
}
