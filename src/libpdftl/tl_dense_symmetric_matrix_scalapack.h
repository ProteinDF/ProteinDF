#ifndef TL_DENSE_SYMMETRIC_MATRIX_SCALAPACK_H
#define TL_DENSE_SYMMETRIC_MATRIX_SCALAPACK_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrixObject;
class TlDenseGeneralMatrix_Scalapack;
class TlDenseVector_Scalapack;
class TlDenseVector_Lapack;
class TlSparseMatrix;

class TlDenseSymmetricMatrix_Scalapack : public TlDenseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlDenseSymmetricMatrix_Scalapack(
        const TlMatrixObject::index_type dim = 1);
    TlDenseSymmetricMatrix_Scalapack(
        const TlDenseSymmetricMatrix_Scalapack& rhs);
    TlDenseSymmetricMatrix_Scalapack(const TlDenseGeneralMatrix_Scalapack& rhs);
    //   TlDenseSymmetricMatrix_Scalapack(const TlDenseVector_Scalapack& v,
    //                                    const TlMatrixObject::index_type dim);
    virtual ~TlDenseSymmetricMatrix_Scalapack();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
   public:
    virtual TlMatrixObject::size_type getNumOfElements() const;
    virtual double getLocal(TlMatrixObject::index_type row,
                            TlMatrixObject::index_type col) const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    // operator TlDenseGeneralMatrix_Scalapack() const;

    TlDenseSymmetricMatrix_Scalapack& operator=(
        const TlDenseSymmetricMatrix_Scalapack& rhs);

    const TlDenseSymmetricMatrix_Scalapack operator+(
        const TlDenseSymmetricMatrix_Scalapack& rhs) const;
    const TlDenseSymmetricMatrix_Scalapack operator-(
        const TlDenseSymmetricMatrix_Scalapack& rhs) const;
    const TlDenseGeneralMatrix_Scalapack operator*(
        const TlDenseSymmetricMatrix_Scalapack& rhs) const;

    TlDenseSymmetricMatrix_Scalapack& operator+=(
        const TlDenseSymmetricMatrix_Scalapack& rhs);
    TlDenseSymmetricMatrix_Scalapack& operator-=(
        const TlDenseSymmetricMatrix_Scalapack& rhs);
    TlDenseSymmetricMatrix_Scalapack& operator*=(const double coef);
    TlDenseSymmetricMatrix_Scalapack& operator/=(const double coef);
    TlDenseSymmetricMatrix_Scalapack& operator*=(
        const TlDenseSymmetricMatrix_Scalapack& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
    const TlDenseSymmetricMatrix_Scalapack& dotInPlace(
        const TlDenseSymmetricMatrix_Scalapack& rhs);

    bool eig(TlDenseVector_Lapack* pEigVal,
             TlDenseGeneralMatrix_Scalapack* pEigVec) const;

    TlDenseSymmetricMatrix_Scalapack inverse() const;

    bool getSparseMatrix(TlSparseMatrix* pMatrix,
                         bool isFinalize = false) const;
    void mergeSparseMatrix(const TlSparseMatrix& M);

    std::vector<TlMatrixObject::index_type> getRowIndexTable() const;
    std::vector<TlMatrixObject::index_type> getColIndexTable() const;
    void getLocalMatrix(TlDenseGeneralMatrixObject* pOutMatrix) const;

   public:
    // double* data();
    // const double* data() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
   public:
    virtual bool load(const std::string& filePath);
    virtual bool save(const std::string& filePath) const;

    virtual void dump(TlDenseVector_Scalapack* v) const;
    virtual void restore(const TlDenseVector_Scalapack& v);

    // ---------------------------------------------------------------------------
    // friend
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_Scalapack;

    // matrix x matrix
    friend TlDenseGeneralMatrix_Scalapack operator*(
        const TlDenseSymmetricMatrix_Scalapack& rhs1,
        const TlDenseGeneralMatrix_Scalapack& rhs2);
    friend TlDenseGeneralMatrix_Scalapack operator*(
        const TlDenseGeneralMatrix_Scalapack& rhs1,
        const TlDenseSymmetricMatrix_Scalapack& rhs2);

    // matrix x vector
    friend TlDenseVector_Scalapack operator*(
        const TlDenseSymmetricMatrix_Scalapack& rhs1,
        const TlDenseVector_Scalapack& rhs2);
    friend TlDenseVector_Scalapack operator*(
        const TlDenseVector_Scalapack& rhs1,
        const TlDenseSymmetricMatrix_Scalapack& rhs2);
};

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Scalapack operator*(
    const double coef, const TlDenseSymmetricMatrix_Scalapack& matrix);
TlDenseSymmetricMatrix_Scalapack operator*(
    const TlDenseSymmetricMatrix_Scalapack& matrix, const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_SCALAPACK_H
