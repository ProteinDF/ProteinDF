#ifndef TL_DENSE_GENERAL_MATRIX_VIENNACL_H
#define TL_DENSE_GENERAL_MATRIX_VIENNACL_H

#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_ViennaCL;
class TlDenseGeneralMatrix_ImplViennaCL;
class TlDenseVector_ViennaCL;
class TlSparseGeneralMatrix_ViennaCL;

class TlDenseGeneralMatrix_Eigen;

class TlDenseGeneralMatrix_ViennaCL : public TlDenseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlDenseGeneralMatrix_ViennaCL(
        const TlMatrixObject::index_type row = 1,
        const TlMatrixObject::index_type col = 1,
        double const * const pBuf = NULL);
    TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL(const TlDenseSymmetricMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_ImplViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL(const TlSparseGeneralMatrix_ViennaCL& rhs);

#ifdef HAVE_EIGEN
    TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_Eigen& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlDenseGeneralMatrix_ViennaCL();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlDenseGeneralMatrix_ViennaCL& operator=(
        const TlDenseGeneralMatrix_ViennaCL& rhs);

    const TlDenseGeneralMatrix_ViennaCL operator+(
        const TlDenseGeneralMatrix_ViennaCL& rhs) const;
    const TlDenseGeneralMatrix_ViennaCL operator-(
        const TlDenseGeneralMatrix_ViennaCL& rhs) const;
    const TlDenseGeneralMatrix_ViennaCL operator*(
        const TlDenseGeneralMatrix_ViennaCL& rhs) const;

    TlDenseGeneralMatrix_ViennaCL& operator+=(
        const TlDenseGeneralMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL& operator-=(
        const TlDenseGeneralMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL& operator*=(const double coef);
    TlDenseGeneralMatrix_ViennaCL& operator/=(const double coef);
    TlDenseGeneralMatrix_ViennaCL& operator*=(
        const TlDenseGeneralMatrix_ViennaCL& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
   public:
    virtual double sum() const;
    // virtual double trace() const;
    virtual double getRMS() const;
    // virtual double getMaxAbsoluteElement(
    //  TlMatrixObject::index_type* outRow,
    //  TlMatrixObject::index_type* outCol) const;

    const TlDenseGeneralMatrix_ViennaCL& dotInPlace(
        const TlDenseGeneralMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCL transpose() const;
    TlDenseGeneralMatrix_ViennaCL inverse() const;

    TlDenseGeneralMatrix_ViennaCL& reverseColumns();

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
   protected:

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ViennaCL;
    friend class TlSparseGeneralMatrix_ViennaCL;
    friend class TlDenseGeneralMatrix_Eigen;

    // DM(G) = DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ViennaCL operator*(
        const TlDenseGeneralMatrix_ViennaCL& mat1,
        const TlDenseSymmetricMatrix_ViennaCL& mat2);
    // DM(G) = DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ViennaCL operator*(
        const TlDenseSymmetricMatrix_ViennaCL& mat1,
        const TlDenseGeneralMatrix_ViennaCL& mat2);

    // DV = DM(G) * DV
    friend TlDenseVector_ViennaCL operator*(
        const TlDenseGeneralMatrix_ViennaCL& rhs1,
        const TlDenseVector_ViennaCL& rhs2);
    // DV = DV * DM(G)
    friend TlDenseVector_ViennaCL operator*(
        const TlDenseVector_ViennaCL& rhs1,
        const TlDenseGeneralMatrix_ViennaCL& rhs2);
};

TlDenseGeneralMatrix_ViennaCL operator*(
    const double coef, const TlDenseGeneralMatrix_ViennaCL& dm);

#endif  // TL_DENSE_GENERAL_MATRIX_VIENNACL_H
