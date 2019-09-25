#ifndef TL_DENSE_VECTOR_VIENNACL_H
#define TL_DENSE_VECTOR_VIENNACL_H

#include <vector>

#include "tl_dense_vector_object.h"

class TlDenseGeneralMatrix_ViennaCL;
class TlDenseSymmetricMatrix_ViennaCL;
class TlDenseVector_Eigen;
class TlSparseGeneralMatrix_ViennaCL;
class TlSparseSymmetricMatrix_ViennaCL;

class TlDenseVector_ViennaCL : public TlDenseVectorObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlDenseVector_ViennaCL(
        const TlDenseVectorObject::index_type size = 0);
    TlDenseVector_ViennaCL(const TlDenseVector_ViennaCL& rhs);
    TlDenseVector_ViennaCL(const std::vector<double>& rhs);

#ifdef HAVE_EIGEN
    TlDenseVector_ViennaCL(const TlDenseVector_Eigen& rhs);
#endif  // HAVE_EIGEN

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ViennaCL();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlDenseVector_ViennaCL& operator=(const TlDenseVector_ViennaCL& rhs);

    TlDenseVector_ViennaCL& operator+=(const TlDenseVector_ViennaCL& rhs);
    TlDenseVector_ViennaCL& operator-=(const TlDenseVector_ViennaCL& rhs);
    TlDenseVector_ViennaCL& operator*=(const double rhs);
    TlDenseVector_ViennaCL& operator/=(const double rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    TlDenseVector_ViennaCL& dotInPlace(const TlDenseVector_ViennaCL& rhs);

    TlDenseVector_ViennaCL& reverse();

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ViennaCL;
    friend class TlDenseVector_Eigen;

    friend TlDenseVector_ViennaCL operator+(const TlDenseVector_ViennaCL& rhs1,
                                            const TlDenseVector_ViennaCL& rhs2);
    friend TlDenseVector_ViennaCL operator-(const TlDenseVector_ViennaCL& rhs1,
                                            const TlDenseVector_ViennaCL& rhs2);
    friend TlDenseVector_ViennaCL operator*(const TlDenseVector_ViennaCL& rhs1,
                                            const double rhs2);
    friend TlDenseVector_ViennaCL operator*(const double rhs1,
                                            const TlDenseVector_ViennaCL& rhs2);

    friend TlDenseVector_ViennaCL operator*(
        const TlDenseGeneralMatrix_ViennaCL& rhs1,
        const TlDenseVector_ViennaCL& rhs2);
    friend TlDenseVector_ViennaCL operator*(
        const TlDenseVector_ViennaCL& rhs1,
        const TlDenseGeneralMatrix_ViennaCL& rhs2);

    friend TlDenseVector_ViennaCL operator*(
        const TlSparseGeneralMatrix_ViennaCL& mat,
        const TlDenseVector_ViennaCL& vtr);
    friend TlDenseVector_ViennaCL operator*(
        const TlSparseSymmetricMatrix_ViennaCL& mat,
        const TlDenseVector_ViennaCL& vtr);
};

#endif  // TL_DENSE_VECTOR_VIENNACL_H
