#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_EIGEN_FLOAT_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_EIGEN_FLOAT_H

#include <Eigen/Core>

#include "tl_dense_matrix_impl_object.h"
#include "tl_dense_vector_impl_eigen_float.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlDenseGeneralMatrix_ImplEigenFloat;
class TlDenseSymmetricMatrix_ImplEigenFloat;
class TlSparseGeneralMatrix_ImplEigenFloat;
class TlSparseSymmetricMatrix_ImplEigenFloat;
class TlDenseGeneralMatrix_ImplViennaCL;

class TlDenseGeneralMatrix_ImplEigenFloat : public TlDenseMatrix_ImplObject {
public:
    typedef Eigen::MatrixXf MatrixDataType;
    typedef Eigen::Map<MatrixDataType> MapType;
    typedef Eigen::Map<const MatrixDataType> MapTypeConst;

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_ImplEigenFloat(
        const TlMatrixObject::index_type row = 0,
        const TlMatrixObject::index_type col = 0,
        double const* const pBuf = NULL);
    TlDenseGeneralMatrix_ImplEigenFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat(const TlDenseSymmetricMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat(const MatrixDataType& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& sm);

#ifdef HAVE_VIENNACL
    TlDenseGeneralMatrix_ImplEigenFloat(
        const TlDenseGeneralMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlDenseGeneralMatrix_ImplEigenFloat();

    operator std::vector<double>() const;

    // ---------------------------------------------------------------------------
    // transformation
    // ---------------------------------------------------------------------------
    void vtr2mat(const std::vector<double>& vtr);

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    virtual TlMatrixObject::index_type getNumOfRows() const;
    virtual TlMatrixObject::index_type getNumOfCols() const;
    virtual void resize(const TlMatrixObject::index_type row,
                        const TlMatrixObject::index_type col);

    virtual double get(const TlMatrixObject::index_type row,
                       const TlMatrixObject::index_type col) const;

    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseGeneralMatrix_ImplEigenFloat& operator=(
        const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat& operator+=(
        const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat& operator-=(
        const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
    TlDenseGeneralMatrix_ImplEigenFloat& operator*=(const double coef);
    TlDenseGeneralMatrix_ImplEigenFloat& operator/=(const double coef);
    TlDenseGeneralMatrix_ImplEigenFloat& operator*=(
        const TlDenseGeneralMatrix_ImplEigenFloat& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    virtual double sum() const;
    virtual double getRMS() const;
    virtual double getMaxAbsoluteElement(
        TlMatrixObject::index_type* outRow,
        TlMatrixObject::index_type* outCol) const;
    virtual void transposeInPlace();

    TlDenseGeneralMatrix_ImplEigenFloat dot(const TlDenseGeneralMatrix_ImplEigenFloat& rhs) const;
    const TlDenseGeneralMatrix_ImplEigenFloat& dotInPlace(
        const TlDenseGeneralMatrix_ImplEigenFloat& rhs);

    TlDenseGeneralMatrix_ImplEigenFloat transpose() const;
    TlDenseGeneralMatrix_ImplEigenFloat inverse() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    virtual void vtr2mat(double const* const pBuf);

protected:
    virtual TlMatrixObject::size_type getNumOfElements() const;
    float* data();
    const float* data() const;

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    mutable MatrixDataType matrix_;

#if __cplusplus >= 201103L
    mutable std::mutex matrix_mutex_;
#endif  // __cplusplus

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_EigenFloat;

    friend class TlDenseSymmetricMatrix_ImplEigenFloat;
    friend class TlSparseGeneralMatrix_ImplEigenFloat;

    friend class TlDenseGeneralMatrix_ImplViennaCL;

    // DM(G) = DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat2);
    // DM(G) = DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

    // DV = DM(G) * DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vec);
    // DV = DV * DM(G)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& vec,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat);

    // v.s. sparse matrix
    // DM(G) = DM(G) * SM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat2);
    // DM(G) = SM(G) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

    // DM(G) = DM(G) * SM(S)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat1,
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat2);
    // DM(G) = SM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat1,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat2);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const double coef, const TlDenseGeneralMatrix_ImplEigenFloat& DM);
TlDenseGeneralMatrix_ImplEigenFloat operator*(
    const TlDenseGeneralMatrix_ImplEigenFloat& DM, const double coef);

#endif  // TL_DENSE_GENERAL_MATRIX_IMPLE_EIGEN_FLOAT_H
