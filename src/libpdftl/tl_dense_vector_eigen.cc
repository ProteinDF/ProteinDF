#include "tl_dense_vector_eigen.h"

#include "tl_dense_vector_impl_eigen.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_Eigen::TlDenseVector_Eigen(TlDenseVectorObject::index_type size) {
    this->pImpl_ = new TlDenseVector_ImplEigen(size);
}

TlDenseVector_Eigen::TlDenseVector_Eigen(const TlDenseVector_Eigen& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigen(
        *dynamic_cast<const TlDenseVector_ImplEigen*>(rhs.pImpl_));
}

TlDenseVector_Eigen::TlDenseVector_Eigen(const std::vector<double>& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigen(rhs);
}

TlDenseVector_Eigen::operator std::vector<double>() const {
    return std::vector<double>(*(dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)));
}

// TlDenseVector_Eigen::TlDenseVector_Eigen(const double* p, size_type size) {
//   this->pImpl_ = new TlDenseVector_ImplEigen(p, size);
// }

TlDenseVector_Eigen::TlDenseVector_Eigen(const TlDenseVector_ImplEigen& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigen(rhs);
}

#ifdef HAVE_VIENNACL
TlDenseVector_Eigen::TlDenseVector_Eigen(const TlDenseVector_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigen(
        *dynamic_cast<const TlDenseVector_ImplViennaCL*>(rhs.pImpl_));
}
#endif  // HAVE_VIENNACL

TlDenseVector_Eigen::~TlDenseVector_Eigen() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_Eigen& TlDenseVector_Eigen::operator=(
    const TlDenseVector_Eigen& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseVector_ImplEigen(
            *dynamic_cast<const TlDenseVector_ImplEigen*>(rhs.pImpl_));
    }

    return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator+=(
    const TlDenseVector_Eigen& rhs) {
    dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
        ->
        operator+=(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator-=(
    const TlDenseVector_Eigen& rhs) {
    dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
        ->
        operator-=(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator*=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)->operator*=(rhs);

    return *this;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::operator/=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)->operator/=(rhs);

    return *this;
}

double TlDenseVector_Eigen::operator*(const TlDenseVector_Eigen& rhs) const {
    return dynamic_cast<const TlDenseVector_ImplEigen*>(this->pImpl_)
        ->
        operator*(*dynamic_cast<const TlDenseVector_ImplEigen*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_Eigen::dot(const TlDenseVector_Eigen& rhs) const {
    double answer =
        dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
            ->dot(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

    return answer;
}

TlDenseVector_Eigen& TlDenseVector_Eigen::dotInPlace(
    const TlDenseVector_Eigen& rhs) {
    dynamic_cast<TlDenseVector_ImplEigen*>(this->pImpl_)
        ->dotInPlace(*dynamic_cast<TlDenseVector_ImplEigen*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_Eigen operator+(const TlDenseVector_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
    TlDenseVector_Eigen answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_Eigen operator-(const TlDenseVector_Eigen& rhs1,
                              const TlDenseVector_Eigen& rhs2) {
    TlDenseVector_Eigen answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                              const double rhs2) {
    TlDenseVector_Eigen answer = rhs1;
    answer *= rhs2;

    return answer;
}

TlDenseVector_Eigen operator*(const double rhs1,
                              const TlDenseVector_Eigen& rhs2) {
    return rhs2 * rhs1;
}
