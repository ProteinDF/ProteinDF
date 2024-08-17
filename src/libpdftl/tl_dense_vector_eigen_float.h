#ifndef TL_DENSE_VECTOR_EIGEN_FLOAT_H
#define TL_DENSE_VECTOR_EIGEN_FLOAT_H

#include <vector>

#include "tl_dense_vector_impl_eigen_float.h"
#include "tl_dense_vector_object.h"
class TlDenseGeneralMatrix_EigenFloat;
class TlDenseSymmetricMatrix_EigenFloat;

class TlSparseGeneralMatrix_EigenFloat;
class TlSparseSymmetricMatrix_EigenFloat;

class TlDenseVector_ViennaCLFloat;

class TlDenseVector_EigenFloat : public TlDenseVectorObject {
public:
    explicit TlDenseVector_EigenFloat(TlDenseVectorObject::index_type size = 0);
    TlDenseVector_EigenFloat(const TlDenseVector_EigenFloat& rhs);
    TlDenseVector_EigenFloat(const std::vector<double>& rhs);
    //   TlDenseVector_EigenFloat(const double* p, const TlDenseVectorObject::size_type length);
    TlDenseVector_EigenFloat(const TlDenseVector_ImplEigenFloat& rhs);

#ifdef HAVE_VIENNACL
    TlDenseVector_EigenFloat(const TlDenseVector_ViennaCLFloat& rhs);
#endif  // HAVE_VIENNACL

    operator std::vector<double>() const;

    virtual ~TlDenseVector_EigenFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_EigenFloat& operator=(const TlDenseVector_EigenFloat& rhs);

    TlDenseVector_EigenFloat& operator+=(const TlDenseVector_EigenFloat& rhs);
    TlDenseVector_EigenFloat& operator-=(const TlDenseVector_EigenFloat& rhs);
    TlDenseVector_EigenFloat& operator*=(const double rhs);
    TlDenseVector_EigenFloat& operator/=(const double rhs);

    double operator*(const TlDenseVector_EigenFloat& rhs) const;

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    double dot(const TlDenseVector_EigenFloat& rhs) const;
    TlDenseVector_EigenFloat& dotInPlace(const TlDenseVector_EigenFloat& rhs);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_EigenFloat;
    friend class TlDenseSymmetricMatrix_EigenFloat;
    friend class TlDenseVector_ViennaCLFloat;

    friend TlDenseVector_EigenFloat operator+(const TlDenseVector_EigenFloat& rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);
    friend TlDenseVector_EigenFloat operator-(const TlDenseVector_EigenFloat& rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);
    friend TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& rhs1,
                                              const double rhs2);
    friend TlDenseVector_EigenFloat operator*(const double rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);

    // DM(G) * DV
    friend TlDenseVector_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);
    // DV * DM(G)
    friend TlDenseVector_EigenFloat operator*(
        const TlDenseVector_EigenFloat& rhs1,
        const TlDenseGeneralMatrix_EigenFloat& rhs2);
    // DM(S) * DV
    friend TlDenseVector_EigenFloat operator*(
        const TlDenseSymmetricMatrix_EigenFloat& rhs1,
        const TlDenseVector_EigenFloat& rhs2);
    // DV * DM(S)
    friend TlDenseVector_EigenFloat operator*(
        const TlDenseVector_EigenFloat& rhs1,
        const TlDenseSymmetricMatrix_EigenFloat& rhs2);

    friend TlDenseVector_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& mat,
                                              const TlDenseVector_EigenFloat& vtr);
    friend TlDenseVector_EigenFloat operator*(
        const TlSparseSymmetricMatrix_EigenFloat& mat,
        const TlDenseVector_EigenFloat& vtr);
};

#endif  // TL_DENSE_VECTOR_EIGEN_FLOAT_H
