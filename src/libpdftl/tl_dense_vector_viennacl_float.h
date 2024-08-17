#ifndef TL_DENSE_VECTOR_VIENNACL_FLOAT_H
#define TL_DENSE_VECTOR_VIENNACL_FLOAT_H

#include <vector>

#include "tl_dense_vector_object.h"

class TlDenseGeneralMatrix_ViennaCLFloat;
class TlDenseSymmetricMatrix_ViennaCLFloat;
class TlDenseVector_EigenFloat;
class TlSparseGeneralMatrix_ViennaCLFloat;
class TlSparseSymmetricMatrix_ViennaCLFloat;

class TlDenseVector_ViennaCLFloat : public TlDenseVectorObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseVector_ViennaCLFloat(const TlDenseVectorObject::index_type size = 0);
    TlDenseVector_ViennaCLFloat(const TlDenseVector_ViennaCLFloat& rhs);
    TlDenseVector_ViennaCLFloat(const std::vector<double>& rhs);

#ifdef HAVE_EIGEN
    TlDenseVector_ViennaCLFloat(const TlDenseVector_EigenFloat& rhs);
#endif  // HAVE_EIGEN

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ViennaCLFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_ViennaCLFloat& operator=(const TlDenseVector_ViennaCLFloat& rhs);

    TlDenseVector_ViennaCLFloat& operator+=(const TlDenseVector_ViennaCLFloat& rhs);
    TlDenseVector_ViennaCLFloat& operator-=(const TlDenseVector_ViennaCLFloat& rhs);
    TlDenseVector_ViennaCLFloat& operator*=(const float rhs);
    TlDenseVector_ViennaCLFloat& operator/=(const float rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    TlDenseVector_ViennaCLFloat& dotInPlace(const TlDenseVector_ViennaCLFloat& rhs);

    TlDenseVector_ViennaCLFloat& reverse();

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ViennaCLFloat;
    friend class TlDenseVector_EigenFloat;

    friend TlDenseVector_ViennaCLFloat operator+(const TlDenseVector_ViennaCLFloat& rhs1,
                                                 const TlDenseVector_ViennaCLFloat& rhs2);
    friend TlDenseVector_ViennaCLFloat operator-(const TlDenseVector_ViennaCLFloat& rhs1,
                                                 const TlDenseVector_ViennaCLFloat& rhs2);
    friend TlDenseVector_ViennaCLFloat operator*(const TlDenseVector_ViennaCLFloat& rhs1,
                                                 const float rhs2);
    friend TlDenseVector_ViennaCLFloat operator*(const double rhs1,
                                                 const TlDenseVector_ViennaCLFloat& rhs2);

    friend TlDenseVector_ViennaCLFloat operator*(const TlDenseGeneralMatrix_ViennaCLFloat& rhs1,
                                                 const TlDenseVector_ViennaCLFloat& rhs2);
    friend TlDenseVector_ViennaCLFloat operator*(const TlDenseVector_ViennaCLFloat& rhs1,
                                                 const TlDenseGeneralMatrix_ViennaCLFloat& rhs2);

    friend TlDenseVector_ViennaCLFloat operator*(const TlSparseGeneralMatrix_ViennaCLFloat& mat,
                                                 const TlDenseVector_ViennaCLFloat& vtr);
    friend TlDenseVector_ViennaCLFloat operator*(const TlSparseSymmetricMatrix_ViennaCLFloat& mat,
                                                 const TlDenseVector_ViennaCLFloat& vtr);
};

#endif  // TL_DENSE_VECTOR_VIENNACL_FLOAT_H
