#ifndef TL_DENSE_GENERAL_MATRIX_EIGEN_FLOAT_H
#define TL_DENSE_GENERAL_MATRIX_EIGEN_FLOAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_eigen_float.h"
#include "tl_dense_general_matrix_object.h"

class TlDenseGeneralMatrix_Eigen;
class TlDenseSymmetricMatrix_EigenFloat;
class TlDenseVector_EigenFloat;
class TlSparseGeneralMatrix_EigenFloat;
class TlSparseSymmetricMatrix_EigenFloat;

class TlCommunicate;
class TlDenseGeneralMatrix_ViennaCLFloat;

class TlDenseGeneralMatrix_EigenFloat : public TlDenseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_EigenFloat(const TlMatrixObject::index_type row = 1,
                                             const TlMatrixObject::index_type col = 1,
                                             double const* const pBuf = NULL);
    TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs);
    TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_Eigen& rhs);
    TlDenseGeneralMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs);
    TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& sm);
#ifdef HAVE_VIENNACL
    TlDenseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
#endif  // HAVE_VIENNACL

    operator std::vector<double>() const;

    virtual ~TlDenseGeneralMatrix_EigenFloat();

    // ---------------------------------------------------------------------------
    // transformation
    // ---------------------------------------------------------------------------
    void vtr2mat(const std::vector<double>& vtr);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseGeneralMatrix_EigenFloat& operator=(const TlDenseGeneralMatrix_EigenFloat& rhs);

    const TlDenseGeneralMatrix_EigenFloat operator+(
        const TlDenseGeneralMatrix_EigenFloat& rhs) const;
    const TlDenseGeneralMatrix_EigenFloat operator-(
        const TlDenseGeneralMatrix_EigenFloat& rhs) const;
    const TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseGeneralMatrix_EigenFloat& rhs) const;

    TlDenseGeneralMatrix_EigenFloat& operator+=(const TlDenseGeneralMatrix_EigenFloat& rhs);
    TlDenseGeneralMatrix_EigenFloat& operator-=(const TlDenseGeneralMatrix_EigenFloat& rhs);
    TlDenseGeneralMatrix_EigenFloat& operator*=(const double coef);
    TlDenseGeneralMatrix_EigenFloat& operator/=(const double coef);
    TlDenseGeneralMatrix_EigenFloat& operator*=(const TlDenseGeneralMatrix_EigenFloat& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    double sum() const;
    double getRMS() const;

    TlDenseGeneralMatrix_EigenFloat dot(const TlDenseGeneralMatrix_EigenFloat& rhs) const;
    const TlDenseGeneralMatrix_EigenFloat& dotInPlace(
        const TlDenseGeneralMatrix_EigenFloat& rhs);

    TlDenseGeneralMatrix_EigenFloat transpose() const;
    TlDenseGeneralMatrix_EigenFloat inverse() const;

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

    friend class TlDenseGeneralMatrix_Eigen;
    friend class TlDenseSymmetricMatrix_EigenFloat;
    friend class TlDenseGeneralMatrix_ViennaCLFloat;
    friend class TlSparseGeneralMatrix_EigenFloat;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseGeneralMatrix_EigenFloat& mat1,
        const TlDenseSymmetricMatrix_EigenFloat& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseSymmetricMatrix_EigenFloat& mat1,
        const TlDenseGeneralMatrix_EigenFloat& mat2);

    // DM(G) * DV
    friend TlDenseVector_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& rhs1,
                                              const TlDenseVector_EigenFloat& rhs2);
    // DV * DM(G)
    friend TlDenseVector_EigenFloat operator*(const TlDenseVector_EigenFloat& rhs1,
                                              const TlDenseGeneralMatrix_EigenFloat& rhs2);

    // DM(G) * SM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseGeneralMatrix_EigenFloat& mat1,
        const TlSparseGeneralMatrix_EigenFloat& mat2);
    // SM(G) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlSparseGeneralMatrix_EigenFloat& mat1,
        const TlDenseGeneralMatrix_EigenFloat& mat2);

    // DM(G) * SM(S)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseGeneralMatrix_EigenFloat& mat1,
        const TlSparseSymmetricMatrix_EigenFloat& mat2);
    // SM(S) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlSparseSymmetricMatrix_EigenFloat& mat1,
        const TlDenseGeneralMatrix_EigenFloat& mat2);
};

TlDenseGeneralMatrix_EigenFloat operator*(const double coef,
                                          const TlDenseGeneralMatrix_EigenFloat& DM);
TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& DM,
                                          const double coef);

#endif  // TL_DENSE_GENERAL_MATRIX_EIGEN_FLOAT_H
