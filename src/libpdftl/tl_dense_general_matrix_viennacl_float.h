#ifndef TL_DENSE_GENERAL_MATRIX_VIENNACL_FLOAT_H
#define TL_DENSE_GENERAL_MATRIX_VIENNACL_FLOAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_object.h"

class TlDenseGeneralMatrix_ViennaCL;
class TlDenseSymmetricMatrix_ViennaCLFloat;
class TlDenseGeneralMatrix_ImplViennaCLFloat;
class TlDenseVector_ViennaCLFloat;
class TlSparseGeneralMatrix_ViennaCLFloat;

#ifdef HAVE_EIGEN
class TlDenseGeneralMatrix_EigenFloat;
#endif  // HAVE_EIGEN

class TlDenseGeneralMatrix_ViennaCLFloat : public TlDenseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_ViennaCLFloat(const TlMatrixObject::index_type row = 1,
                                                const TlMatrixObject::index_type col = 1,
                                                double const* const pBuf = NULL);
    TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCL& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_ViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlDenseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_EigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlDenseGeneralMatrix_ViennaCLFloat();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseGeneralMatrix_ViennaCLFloat& operator=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);

    const TlDenseGeneralMatrix_ViennaCLFloat operator+(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const;
    const TlDenseGeneralMatrix_ViennaCLFloat operator-(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const;
    const TlDenseGeneralMatrix_ViennaCLFloat operator*(const TlDenseGeneralMatrix_ViennaCLFloat& rhs) const;

    TlDenseGeneralMatrix_ViennaCLFloat& operator+=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat& operator-=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat& operator*=(const double coef);
    TlDenseGeneralMatrix_ViennaCLFloat& operator/=(const double coef);
    TlDenseGeneralMatrix_ViennaCLFloat& operator*=(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);

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

    const TlDenseGeneralMatrix_ViennaCLFloat& dotInPlace(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ViennaCLFloat transpose() const;
    TlDenseGeneralMatrix_ViennaCLFloat inverse() const;

    TlDenseGeneralMatrix_ViennaCLFloat& reverseColumns();

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
protected:
    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_ViennaCL;
    friend class TlDenseSymmetricMatrix_ViennaCLFloat;
    friend class TlSparseGeneralMatrix_ViennaCLFloat;
    friend class TlDenseGeneralMatrix_EigenFloat;

    // DM(G) = DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ViennaCLFloat& mat1,
        const TlDenseSymmetricMatrix_ViennaCLFloat& mat2);
    // DM(G) = DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ViennaCLFloat operator*(
        const TlDenseSymmetricMatrix_ViennaCLFloat& mat1,
        const TlDenseGeneralMatrix_ViennaCLFloat& mat2);

    // DV = DM(G) * DV
    friend TlDenseVector_ViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ViennaCLFloat& rhs1,
        const TlDenseVector_ViennaCLFloat& rhs2);
    // DV = DV * DM(G)
    friend TlDenseVector_ViennaCLFloat operator*(
        const TlDenseVector_ViennaCLFloat& rhs1,
        const TlDenseGeneralMatrix_ViennaCLFloat& rhs2);
};

TlDenseGeneralMatrix_ViennaCLFloat operator*(
    const double coef, const TlDenseGeneralMatrix_ViennaCLFloat& dm);

#endif  // TL_DENSE_GENERAL_MATRIX_VIENNACL_FLOAT_H
