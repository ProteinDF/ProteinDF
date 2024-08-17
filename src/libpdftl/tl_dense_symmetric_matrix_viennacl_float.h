#ifndef TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H
#define TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_ViennaCLFloat;
class TlDenseVector_ViennaCLFloat;
class TlDenseSymmetricMatrix_EigenFloat;
class TlSparseSymmetricMatrix_ViennaCLFloat;

class TlDenseSymmetricMatrix_ViennaCLFloat : public TlDenseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    //
    // ---------------------------------------------------------------------------
public:
    enum EIG_METHOD { EIG_QR,
                      EIG_POWERITERATION };

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_ViennaCLFloat(const TlMatrixObject::index_type dim = 1);
    TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ViennaCLFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
#ifdef HAVE_EIGEN
    TlDenseSymmetricMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlDenseSymmetricMatrix_ViennaCLFloat();

    virtual void vtr2mat(const std::vector<double>& vtr);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseSymmetricMatrix_ViennaCLFloat& operator=(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ViennaCLFloat& operator=(
        const TlDenseSymmetricMatrix_EigenFloat& rhs);

    const TlDenseSymmetricMatrix_ViennaCLFloat operator+(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const;
    const TlDenseSymmetricMatrix_ViennaCLFloat operator-(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const;
    const TlDenseSymmetricMatrix_ViennaCLFloat operator*(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs) const;

    TlDenseSymmetricMatrix_ViennaCLFloat& operator+=(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ViennaCLFloat& operator-=(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);
    TlDenseSymmetricMatrix_ViennaCLFloat& operator*=(const float coef);
    TlDenseSymmetricMatrix_ViennaCLFloat& operator/=(const float coef);
    TlDenseSymmetricMatrix_ViennaCLFloat& operator*=(
        const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    const TlDenseSymmetricMatrix_ViennaCLFloat& dotInPlace(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);

    bool eig(TlDenseVector_ViennaCLFloat* pEigVal, TlDenseGeneralMatrix_ViennaCLFloat* pEigVec,
             EIG_METHOD eigMethod = EIG_QR) const;
    TlDenseSymmetricMatrix_ViennaCLFloat inverse() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
    virtual bool load(const std::string& filePath);
    virtual bool save(const std::string& filePath) const;

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
protected:
    bool eig_powerIteration(TlDenseVector_ViennaCLFloat* pEigVal,
                            TlDenseGeneralMatrix_ViennaCLFloat* pEigVec) const;

    bool eig_QR(TlDenseVector_ViennaCLFloat* pEigVal,
                TlDenseGeneralMatrix_ViennaCLFloat* pEigVec) const;

    // ---------------------------------------------------------------------------
    // Friends
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_ViennaCLFloat;
    friend class TlDenseSymmetricMatrix_EigenFloat;
    friend class TlSparseSymmetricMatrix_ViennaCLFloat;

    friend TlDenseGeneralMatrix_ViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ViennaCLFloat& mat1,
        const TlDenseSymmetricMatrix_ViennaCLFloat& mat2);
    friend TlDenseGeneralMatrix_ViennaCLFloat operator*(
        const TlDenseSymmetricMatrix_ViennaCLFloat& mat1,
        const TlDenseGeneralMatrix_ViennaCLFloat& mat2);
};

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_ViennaCLFloat operator*(
    const float coef, const TlDenseSymmetricMatrix_ViennaCLFloat& dm);
TlDenseSymmetricMatrix_ViennaCLFloat operator*(
    const TlDenseSymmetricMatrix_ViennaCLFloat& dm, const float coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H
