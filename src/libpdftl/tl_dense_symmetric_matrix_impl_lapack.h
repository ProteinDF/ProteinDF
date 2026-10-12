#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H

#include "tl_dense_general_matrix_impl_lapack.h"
class TlDenseVector_ImplLapack;

class TlDenseSymmetricMatrix_ImplLapack
    : public TlDenseGeneralMatrix_ImplLapack {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    // column-major
    // a11, a21, a31, a41, ..., am1, a12, a22, a32, ...
    explicit TlDenseSymmetricMatrix_ImplLapack(
        const TlMatrixObject::index_type dim = 0,
        double const* const pBuf = NULL);
    TlDenseSymmetricMatrix_ImplLapack(
        const TlDenseSymmetricMatrix_ImplLapack& rhs);
    TlDenseSymmetricMatrix_ImplLapack(
        const TlDenseGeneralMatrix_ImplLapack& rhs);
    virtual ~TlDenseSymmetricMatrix_ImplLapack();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
   public:
    virtual void resize(TlMatrixObject::index_type row,
                        TlMatrixObject::index_type col);

    virtual TlMatrixObject::index_type getRowVector(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type length, double* pBuf) const;
    virtual TlMatrixObject::index_type getColVector(const TlMatrixObject::index_type col,
                                                    const TlMatrixObject::index_type length, double* pBuf) const;
    virtual std::vector<double> getRowVector(const TlMatrixObject::index_type row) const;
    virtual std::vector<double> getColVector(const TlMatrixObject::index_type col) const;

    virtual TlMatrixObject::index_type setRowVector(const TlMatrixObject::index_type row,
                                                    const TlMatrixObject::index_type length, const double* pBuf);
    virtual TlMatrixObject::index_type setColVector(const TlMatrixObject::index_type col,
                                                    const TlMatrixObject::index_type length, const double* pBuf);
    virtual void setRowVector(const TlMatrixObject::index_type row, const std::vector<double>& v);
    virtual void setColVector(const TlMatrixObject::index_type col, const std::vector<double>& v);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
       TlDenseSymmetricMatrix_ImplLapack& operator=(
           const TlDenseSymmetricMatrix_ImplLapack& rhs);
       TlDenseSymmetricMatrix_ImplLapack& operator*=(const double coef);
       // const TlDenseGeneralMatrix_ImplLapack operator*=(
       //   const TlDenseSymmetricMatrix_ImplLapack& rhs);

       // ---------------------------------------------------------------------------
       // operations
       // ---------------------------------------------------------------------------
   public:
    // virtual double sum() const;
    // virtual double getRMS() const;
    // virtual double getMaxAbsoluteElement(
    //     TlMatrixObject::index_type* outRow,
    //     TlMatrixObject::index_type* outCol) const;

    // const TlDenseGeneralMatrix_ImplLapack& dotInPlace(
    //     const TlDenseGeneralMatrix_ImplLapack& rhs);
    TlDenseSymmetricMatrix_ImplLapack transpose() const;
    virtual void transposeInPlace();
    TlDenseSymmetricMatrix_ImplLapack inverse() const;

    bool eig(TlDenseVector_ImplLapack* pEigVal,
             TlDenseGeneralMatrix_ImplLapack* pEigVec) const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
   public:
    // virtual void dump(double* buf, const std::size_t size) const;
    // virtual void restore(const double* buf, const std::size_t size);

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
   protected:
    virtual TlMatrixObject::size_type getNumOfElements() const;

    virtual TlMatrixObject::size_type index(
        TlMatrixObject::index_type row, TlMatrixObject::index_type col) const;

    virtual void vtr2mat(double const* const pBuf);

    // ---------------------------------------------------------------------------
    // private
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // friends
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_Lapack;

    friend TlDenseGeneralMatrix_ImplLapack operator*(
        const TlDenseSymmetricMatrix_ImplLapack& rhs1,
        const TlDenseGeneralMatrix_ImplLapack& rhs2);
    friend TlDenseGeneralMatrix_ImplLapack operator*(
        const TlDenseGeneralMatrix_ImplLapack& rhs1,
        const TlDenseSymmetricMatrix_ImplLapack& rhs2);

    friend TlDenseVector_ImplLapack operator*(
        const TlDenseSymmetricMatrix_ImplLapack& mat,
        const TlDenseVector_ImplLapack& vec);
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H
