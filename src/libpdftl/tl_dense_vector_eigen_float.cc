#include "tl_dense_vector_eigen_float.h"

#include "tl_dense_vector_impl_eigen_float.h"

#ifdef HAVE_VIENNACL
#include "tl_dense_vector_impl_viennacl.h"
#include "tl_dense_vector_viennacl.h"
#endif  // HAVE_VIENNACL

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(TlDenseVectorObject::index_type size) {
    this->pImpl_ = new TlDenseVector_ImplEigenFloat(size);
}

TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(const TlDenseVector_EigenFloat& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigenFloat(*dynamic_cast<const TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));
}

TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(const std::vector<double>& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigenFloat(rhs);
}

TlDenseVector_EigenFloat::operator std::vector<double>() const {
    return std::vector<double>(*(dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)));
}

// TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(const double* p, size_type size) {
//   this->pImpl_ = new TlDenseVector_ImplEigenFloat(p, size);
// }

TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(const TlDenseVector_ImplEigenFloat& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigenFloat(rhs);
}

#ifdef HAVE_VIENNACL
TlDenseVector_EigenFloat::TlDenseVector_EigenFloat(const TlDenseVector_ViennaCL& rhs) {
    this->pImpl_ = new TlDenseVector_ImplEigenFloat(
        *dynamic_cast<const TlDenseVector_ImplViennaCL*>(rhs.pImpl_));
}
#endif  // HAVE_VIENNACL

TlDenseVector_EigenFloat::~TlDenseVector_EigenFloat() {
    delete this->pImpl_;
    this->pImpl_ = NULL;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::operator=(
    const TlDenseVector_EigenFloat& rhs) {
    if (this != &rhs) {
        delete this->pImpl_;
        this->pImpl_ = NULL;

        this->pImpl_ = new TlDenseVector_ImplEigenFloat(*dynamic_cast<const TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));
    }

    return *this;
}

TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::operator+=(
    const TlDenseVector_EigenFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->operator+=(*dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::operator-=(
    const TlDenseVector_EigenFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->operator-=(*dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::operator*=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->operator*=(rhs);

    return *this;
}

TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::operator/=(const double rhs) {
    dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->operator/=(rhs);

    return *this;
}

double TlDenseVector_EigenFloat::operator*(const TlDenseVector_EigenFloat& rhs) const {
    return dynamic_cast<const TlDenseVector_ImplEigenFloat*>(this->pImpl_)->operator*(*dynamic_cast<const TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVector_EigenFloat::dot(const TlDenseVector_EigenFloat& rhs) const {
    double answer = dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->dot(*dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));

    return answer;
}

TlDenseVector_EigenFloat& TlDenseVector_EigenFloat::dotInPlace(const TlDenseVector_EigenFloat& rhs) {
    dynamic_cast<TlDenseVector_ImplEigenFloat*>(this->pImpl_)->dotInPlace(*dynamic_cast<TlDenseVector_ImplEigenFloat*>(rhs.pImpl_));

    return *this;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
TlDenseVector_EigenFloat operator+(const TlDenseVector_EigenFloat& rhs1,
                                   const TlDenseVector_EigenFloat& rhs2) {
    TlDenseVector_EigenFloat answer = rhs1;
    answer += rhs2;

    return answer;
}

TlDenseVector_EigenFloat operator-(const TlDenseVector_EigenFloat& rhs1,
                                   const TlDenseVector_EigenFloat& rhs2) {
    TlDenseVector_EigenFloat answer = rhs1;
    answer -= rhs2;

    return answer;
}

TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& rhs1,
                                   const double rhs2) {
    TlDenseVector_EigenFloat answer = rhs1;
    answer *= rhs2;

    return answer;
}

TlDenseVector_EigenFloat operator*(const double rhs1,
                                   const TlDenseVector_EigenFloat& rhs2) {
    return rhs2 * rhs1;
}
