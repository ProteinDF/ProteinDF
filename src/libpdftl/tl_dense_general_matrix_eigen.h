#ifndef TL_DENSE_GENERAL_MATRIX_EIGEN_H
#define TL_DENSE_GENERAL_MATRIX_EIGEN_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_object.h"

class TlCommunicate;
class TlDenseGeneralMatrix_EigenFloat;
class TlDenseSymmetricMatrix_Eigen;
class TlDenseVector_Eigen;
class TlDenseGeneralMatrix_ImplEigen;
class TlDenseGeneralMatrix_ViennaCL;
class TlSparseGeneralMatrix_Eigen;
class TlSparseSymmetricMatrix_Eigen;

class TlDenseGeneralMatrix_Eigen : public TlDenseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_Eigen(const TlMatrixObject::index_type row = 1,
                                        const TlMatrixObject::index_type col = 1,
                                        double const* const pBuf = NULL);
    TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs);
    TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_EigenFloat& rhs);
    TlDenseGeneralMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs);
    TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_ImplEigen& rhs);
    TlDenseGeneralMatrix_Eigen(const TlSparseGeneralMatrix_Eigen& sm);
#ifdef HAVE_VIENNACL
    TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_ViennaCL& rhs);
#endif  // HAVE_VIENNACL

    operator std::vector<double>() const;

    virtual ~TlDenseGeneralMatrix_Eigen();

    // ---------------------------------------------------------------------------
    // transformation
    // ---------------------------------------------------------------------------
    void vtr2mat(const std::vector<double>& vtr);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseGeneralMatrix_Eigen& operator=(const TlDenseGeneralMatrix_Eigen& rhs);

    const TlDenseGeneralMatrix_Eigen operator+(
        const TlDenseGeneralMatrix_Eigen& rhs) const;
    const TlDenseGeneralMatrix_Eigen operator-(
        const TlDenseGeneralMatrix_Eigen& rhs) const;
    const TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& rhs) const;

    TlDenseGeneralMatrix_Eigen& operator+=(const TlDenseGeneralMatrix_Eigen& rhs);
    TlDenseGeneralMatrix_Eigen& operator-=(const TlDenseGeneralMatrix_Eigen& rhs);
    TlDenseGeneralMatrix_Eigen& operator*=(const double coef);
    TlDenseGeneralMatrix_Eigen& operator/=(const double coef);
    TlDenseGeneralMatrix_Eigen& operator*=(const TlDenseGeneralMatrix_Eigen& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    double sum() const;
    double getRMS() const;

    TlDenseGeneralMatrix_Eigen dot(const TlDenseGeneralMatrix_Eigen& rhs) const;
    const TlDenseGeneralMatrix_Eigen& dotInPlace(
        const TlDenseGeneralMatrix_Eigen& rhs);

    TlDenseGeneralMatrix_Eigen transpose() const;
    TlDenseGeneralMatrix_Eigen inverse() const;

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

    friend class TlDenseGeneralMatrix_EigenFloat;
    friend class TlDenseSymmetricMatrix_Eigen;
    friend class TlDenseGeneralMatrix_ViennaCL;
    friend class TlSparseGeneralMatrix_Eigen;

    // DM(G) * DM(S)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& mat1,
        const TlDenseSymmetricMatrix_Eigen& mat2);
    // DM(S) * DM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseSymmetricMatrix_Eigen& mat1,
        const TlDenseGeneralMatrix_Eigen& mat2);

    // DM(G) * DV
    friend TlDenseVector_Eigen operator*(const TlDenseGeneralMatrix_Eigen& rhs1,
                                         const TlDenseVector_Eigen& rhs2);
    // DV * DM(G)
    friend TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                                         const TlDenseGeneralMatrix_Eigen& rhs2);

    // DM(G) * SM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& mat1,
        const TlSparseGeneralMatrix_Eigen& mat2);
    // SM(G) * DM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlSparseGeneralMatrix_Eigen& mat1,
        const TlDenseGeneralMatrix_Eigen& mat2);

    // DM(G) * SM(S)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& mat1,
        const TlSparseSymmetricMatrix_Eigen& mat2);
    // SM(S) * DM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlSparseSymmetricMatrix_Eigen& mat1,
        const TlDenseGeneralMatrix_Eigen& mat2);
};

TlDenseGeneralMatrix_Eigen operator*(const double coef,
                                     const TlDenseGeneralMatrix_Eigen& DM);
TlDenseGeneralMatrix_Eigen operator*(const TlDenseGeneralMatrix_Eigen& DM,
                                     const double coef);

#endif  // TL_DENSE_GENERAL_MATRIX_EIGEN_H
