#ifndef TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H
#define TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_symmetric_matrix_object.h"

class TlCommunicate;
class TlDenseGeneralMatrix_Eigen;
class TlDenseSymmetricMatrix_EigenFloat;
class TlSparseSymmetricMatrix_Eigen;
class TlDenseVector_Eigen;

#ifdef HAVE_VIENNACL
class TlDenseSymmetricMatrix_ViennaCL;
#endif  // HAVE_VIENNACL

class TlDenseSymmetricMatrix_Eigen : public TlDenseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_Eigen(const TlMatrixObject::index_type dim = 1,
                                          double const* const pBuf = NULL);
    TlDenseSymmetricMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_Eigen(const TlDenseSymmetricMatrix_EigenFloat& rhs);
    TlDenseSymmetricMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_Eigen(const TlSparseSymmetricMatrix_Eigen& sm);

#ifdef HAVE_VIENNACL
    TlDenseSymmetricMatrix_Eigen(const TlDenseSymmetricMatrix_ViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlDenseSymmetricMatrix_Eigen();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseSymmetricMatrix_Eigen& operator=(
        const TlDenseSymmetricMatrix_Eigen& rhs);

    const TlDenseSymmetricMatrix_Eigen operator+(
        const TlDenseSymmetricMatrix_Eigen& rhs) const;
    const TlDenseSymmetricMatrix_Eigen operator-(
        const TlDenseSymmetricMatrix_Eigen& rhs) const;
    const TlDenseSymmetricMatrix_Eigen operator*(
        const TlDenseSymmetricMatrix_Eigen& rhs) const;

    TlDenseSymmetricMatrix_Eigen& operator+=(
        const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_Eigen& operator-=(
        const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseSymmetricMatrix_Eigen& operator*=(const double coef);
    TlDenseSymmetricMatrix_Eigen& operator/=(const double coef);
    TlDenseSymmetricMatrix_Eigen& operator*=(
        const TlDenseSymmetricMatrix_Eigen& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    // double sum() const;
    const TlDenseSymmetricMatrix_Eigen& dotInPlace(
        const TlDenseSymmetricMatrix_Eigen& rhs);

    bool eig(TlDenseVector_Eigen* pEigVal,
             TlDenseGeneralMatrix_Eigen* pEigVec) const;
    bool diagonal(TlDenseVector_Eigen* pEigVal,
                  TlDenseGeneralMatrix_Eigen* pEigVec) const {
        return this->eig(pEigVal, pEigVec);
    }

    TlDenseSymmetricMatrix_Eigen inverse() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    TlMatrixObject::size_type getNumOfElements() const;
    double* data();
    const double* data() const;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlCommunicate;
    friend class TlDenseGeneralMatrix_Eigen;
    friend class TlDenseSymmetricMatrix_EigenFloat;
    friend class TlSparseSymmetricMatrix_Eigen;
    friend class TlDenseSymmetricMatrix_ViennaCL;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& mat1,
        const TlDenseSymmetricMatrix_Eigen& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseSymmetricMatrix_Eigen& mat1,
        const TlDenseGeneralMatrix_Eigen& mat2);

    // DM(S) * DV
    friend TlDenseVector_Eigen operator*(
        const TlDenseSymmetricMatrix_Eigen& rhs1,
        const TlDenseVector_Eigen& rhs2);
    // DV * DM(S)
    friend TlDenseVector_Eigen operator*(
        const TlDenseVector_Eigen& rhs1,
        const TlDenseSymmetricMatrix_Eigen& rhs2);
};

TlDenseSymmetricMatrix_Eigen operator*(const double coef,
                                       const TlDenseSymmetricMatrix_Eigen& DM);

TlDenseSymmetricMatrix_Eigen operator*(const TlDenseSymmetricMatrix_Eigen& DM,
                                       const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H
