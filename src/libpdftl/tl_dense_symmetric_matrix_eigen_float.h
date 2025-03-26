#ifndef TL_DENSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H
#define TL_DENSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_symmetric_matrix_object.h"

class TlCommunicate;
class TlDenseGeneralMatrix_EigenFloat;
class TlDenseSymmetricMatrix_Eigen;
class TlDenseVector_EigenFloat;
class TlSparseSymmetricMatrix_EigenFloat;

class TlDenseSymmetricMatrix_ViennaCLFloat;

class TlDenseSymmetricMatrix_EigenFloat : public TlDenseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_EigenFloat(const TlMatrixObject::index_type dim = 1,
                                               double const* const pBuf = NULL);
    TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs);
    TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs);
    TlDenseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& sm);

#ifdef HAVE_VIENNACL
    TlDenseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlDenseSymmetricMatrix_EigenFloat();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseSymmetricMatrix_EigenFloat& operator=(const TlDenseSymmetricMatrix_EigenFloat& rhs);

    const TlDenseSymmetricMatrix_EigenFloat operator+(const TlDenseSymmetricMatrix_EigenFloat& rhs) const;
    const TlDenseSymmetricMatrix_EigenFloat operator-(const TlDenseSymmetricMatrix_EigenFloat& rhs) const;
    const TlDenseSymmetricMatrix_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& rhs) const;

    TlDenseSymmetricMatrix_EigenFloat& operator+=(const TlDenseSymmetricMatrix_EigenFloat& rhs);
    TlDenseSymmetricMatrix_EigenFloat& operator-=(const TlDenseSymmetricMatrix_EigenFloat& rhs);
    TlDenseSymmetricMatrix_EigenFloat& operator*=(const double coef);
    TlDenseSymmetricMatrix_EigenFloat& operator/=(const double coef);
    TlDenseSymmetricMatrix_EigenFloat& operator*=(const TlDenseSymmetricMatrix_EigenFloat& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    // double sum() const;
    const TlDenseSymmetricMatrix_EigenFloat& dotInPlace(const TlDenseSymmetricMatrix_EigenFloat& rhs);

    bool eig(TlDenseVector_EigenFloat* pEigVal, TlDenseGeneralMatrix_EigenFloat* pEigVec) const;
    bool diagonal(TlDenseVector_EigenFloat* pEigVal, TlDenseGeneralMatrix_EigenFloat* pEigVec) const {
        return this->eig(pEigVal, pEigVec);
    }

    TlDenseSymmetricMatrix_EigenFloat inverse() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    TlMatrixObject::size_type getNumOfElements() const;
    float* data();
    const float* data() const;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlCommunicate;
    friend class TlDenseGeneralMatrix_EigenFloat;
    friend class TlDenseSymmetricMatrix_Eigen;
    friend class TlSparseSymmetricMatrix_EigenFloat;
    friend class TlDenseSymmetricMatrix_ViennaCLFloat;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& mat1,
                                                     const TlDenseSymmetricMatrix_EigenFloat& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& mat1,
                                                     const TlDenseGeneralMatrix_EigenFloat& mat2);

    // DM(S) * DV
    friend TlDenseVector_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);
    // DV * DM(S)
    friend TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& rhs1,
                                              const TlDenseSymmetricMatrix_EigenFloat& rhs2);
};

TlDenseSymmetricMatrix_EigenFloat operator*(const double coef,
                                            const TlDenseSymmetricMatrix_EigenFloat& DM);

TlDenseSymmetricMatrix_EigenFloat operator*(const TlDenseSymmetricMatrix_EigenFloat& DM,
                                            const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H
