#ifndef TL_DENSE_GENERAL_MATRIX_LAPACK_H
#define TL_DENSE_GENERAL_MATRIX_LAPACK_H

#include "tl_dense_general_matrix_impl_lapack.h"
#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_Lapack;
class TlDenseVector_Lapack;

class TlDenseGeneralMatrix_Lapack : public TlDenseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_Lapack(const TlMatrixObject::index_type row = 1,
                                         const TlMatrixObject::index_type col = 1, double const* const pBuf = NULL);
    TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_Lapack& rhs);
    TlDenseGeneralMatrix_Lapack(const TlDenseSymmetricMatrix_Lapack& rhs);
    TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_ImplLapack& rhs);

    virtual ~TlDenseGeneralMatrix_Lapack();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    virtual TlMatrixObject::size_type getNumOfElements() const;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseGeneralMatrix_Lapack& operator=(const TlDenseGeneralMatrix_Lapack& rhs);

    const TlDenseGeneralMatrix_Lapack operator+(const TlDenseGeneralMatrix_Lapack& rhs) const;
    const TlDenseGeneralMatrix_Lapack operator-(const TlDenseGeneralMatrix_Lapack& rhs) const;
    const TlDenseGeneralMatrix_Lapack operator*(const double coef) const;
    const TlDenseGeneralMatrix_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs) const;

    TlDenseGeneralMatrix_Lapack& operator+=(const TlDenseGeneralMatrix_Lapack& rhs);
    TlDenseGeneralMatrix_Lapack& operator-=(const TlDenseGeneralMatrix_Lapack& rhs);
    TlDenseGeneralMatrix_Lapack& operator*=(const double coef);
    TlDenseGeneralMatrix_Lapack& operator/=(const double coef);
    TlDenseGeneralMatrix_Lapack& operator*=(const TlDenseGeneralMatrix_Lapack& rhs);

public:
    virtual TlMatrixObject::index_type getRowVector(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type length, double* pBuf) const;
    virtual TlMatrixObject::index_type getColVector(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type length, double* pBuf) const;
    virtual std::vector<double> getRowVector(const TlMatrixObject::index_type row) const;
    virtual std::vector<double> getColVector(const TlMatrixObject::index_type col) const;

    virtual TlMatrixObject::index_type setRowVector(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type length, const double* pBuf);
    virtual void setRowVector(const TlMatrixObject::index_type row, const std::vector<double>& v);
    virtual TlMatrixObject::index_type setColVector(const TlMatrixObject::index_type col,
                                                    const TlMatrixObject::index_type length, const double* pBuf);
    virtual void setColVector(const TlMatrixObject::index_type row, const std::vector<double>& v);

    // public:
    // virtual TlMatrixObject::index_type getRowVector(const TlMatrixObject::index_type row,
    //                                                 const TlMatrixObject::index_type length, double* pBuf) const;
    // virtual TlMatrixObject::index_type getColVector(const TlMatrixObject::index_type row,
    //                                                 const TlMatrixObject::index_type length, double* pBuf) const;
    // virtual TlMatrixObject::index_type setRowVector(const TlMatrixObject::index_type row,
    //                                                 const TlMatrixObject::index_type length, const double* pBuf);
    // virtual TlMatrixObject::index_type setColVector(const TlMatrixObject::index_type col,
    //                                                 const TlMatrixObject::index_type length, const double* pBuf);
    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    double sum() const;
    double getRMS() const;

    TlDenseGeneralMatrix_Lapack dot(const TlDenseGeneralMatrix_Lapack& rhs) const;
    const TlDenseGeneralMatrix_Lapack& dotInPlace(const TlDenseGeneralMatrix_Lapack& rhs);

    TlDenseGeneralMatrix_Lapack transpose() const;
    TlDenseGeneralMatrix_Lapack inverse() const;

    // solve Ax = B
    TlDenseGeneralMatrix_Lapack getLeastSquaresSolution(const TlDenseGeneralMatrix_Lapack& B) const;

public:
    double* data();
    const double* data() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
public:
    virtual void dump(TlDenseVector_Lapack* v) const;
    virtual void restore(const TlDenseVector_Lapack& v);

private:
    using TlDenseGeneralMatrixObject::dump;
    using TlDenseGeneralMatrixObject::restore;

    // ---------------------------------------------------------------------------
    // friends
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_Lapack;

    // matrix(sym) x matrix(gen)
    friend TlDenseGeneralMatrix_Lapack operator*(const TlDenseSymmetricMatrix_Lapack& rhs1,
                                                 const TlDenseGeneralMatrix_Lapack& rhs2);
    // matrix(gen) x matrix(sym)
    friend TlDenseGeneralMatrix_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1,
                                                 const TlDenseSymmetricMatrix_Lapack& rhs2);

    friend TlDenseVector_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1, const TlDenseVector_Lapack& rhs2);
    friend TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1, const TlDenseGeneralMatrix_Lapack& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_LAPACK_H
