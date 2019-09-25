#ifndef TL_DENSE_VECTOR_SCALAPACK_H
#define TL_DENSE_VECTOR_SCALAPACK_H

#include <vector>
#include "tl_dense_vector_object.h"

class TlDenseGeneralMatrix_Scalapack;
class TlDenseSymmetricMatrix_Scalapack;
class TlDenseVector_Lapack;

class TlDenseVector_Scalapack : public TlDenseVectorObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlDenseVector_Scalapack(
        const TlDenseVectorObject::index_type size = 0);
    TlDenseVector_Scalapack(const TlDenseVector_Scalapack& rhs);
    TlDenseVector_Scalapack(const TlDenseVector_Lapack& rhs);
    // TlDenseVector_Scalapack(const std::vector<double>& rhs);
    // TlDenseVector_Scalapack(const double* p,
    //                      const TlDenseVectorObject::size_type length);
    virtual ~TlDenseVector_Scalapack();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    std::vector<double> getVector() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlDenseVector_Scalapack& operator=(const TlDenseVector_Scalapack& rhs);

    TlDenseVector_Scalapack& operator+=(const TlDenseVector_Scalapack& rhs);
    TlDenseVector_Scalapack& operator-=(const TlDenseVector_Scalapack& rhs);
    TlDenseVector_Scalapack& operator*=(const double rhs);
    TlDenseVector_Scalapack& operator/=(const double rhs);

    double operator*(const TlDenseVector_Scalapack& rhs) const;

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
   public:
    double dot(const TlDenseVector_Scalapack& rhs) const;
    TlDenseVector_Scalapack& dotInPlace(const TlDenseVector_Scalapack& rhs);

   public:
    // double* data();
    // const double* data() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
   public:
    virtual bool load(const std::string& filePath);
    virtual bool save(const std::string& filePath) const;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_Scalapack;
    friend class TlDenseSymmetricMatrix_Scalapack;

    friend TlDenseVector_Scalapack operator+(
        const TlDenseVector_Scalapack& rhs1,
        const TlDenseVector_Scalapack& rhs2);
    friend TlDenseVector_Scalapack operator-(
        const TlDenseVector_Scalapack& rhs1,
        const TlDenseVector_Scalapack& rhs2);
    friend TlDenseVector_Scalapack operator*(
        const TlDenseVector_Scalapack& rhs1, const double rhs2);
    friend TlDenseVector_Scalapack operator*(
        const double rhs1, const TlDenseVector_Scalapack& rhs2);

    friend TlDenseVector_Scalapack operator*(
        const TlDenseGeneralMatrix_Scalapack& rhs1,
        const TlDenseVector_Scalapack& rhs2);
    friend TlDenseVector_Scalapack operator*(
        const TlDenseVector_Scalapack& rhs1,
        const TlDenseGeneralMatrix_Scalapack& rhs2);

    friend TlDenseVector_Scalapack operator*(
        const TlDenseSymmetricMatrix_Scalapack& rhs1,
        const TlDenseVector_Scalapack& rhs2);
    friend TlDenseVector_Scalapack operator*(
        const TlDenseVector_Scalapack& rhs1,
        const TlDenseSymmetricMatrix_Scalapack& rhs2);
};

#endif  // TL_DENSE_VECTOR_SCALAPACK_H
