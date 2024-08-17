#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include <viennacl/matrix.hpp>

#include "tl_dense_matrix_impl_object.h"

#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif  // HAVE_EIGEN

class TlDenseSymmetricMatrix_ImplViennaCLFloat;
class TlDenseVector_ImplViennaCLFloat;
class TlDenseGeneralMatrix_ImplEigenFloat;
class TlSparseGeneralMatrix_ImplViennaCLFloat;

class TlDenseGeneralMatrix_ImplViennaCLFloat : public TlDenseMatrix_ImplObject {
public:
    typedef viennacl::vector<float> VectorDataType;
    typedef viennacl::matrix<float> MatrixDataType;

#ifdef HAVE_EIGEN
    typedef Eigen::MatrixXf EigenMatrixDataType;
    typedef Eigen::Map<MatrixDataType> EigenMapType;
    typedef Eigen::Map<const MatrixDataType> EigenMapTypeConst;
#endif  // HAVE_EIGEN

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseGeneralMatrix_ImplViennaCLFloat(const TlMatrixObject::index_type row = 0,
                                                    const TlMatrixObject::index_type col = 0,
                                                    double const* const pBuf = NULL);
    TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ImplViennaCLFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
#ifdef HAVE_EIGEN
    TlDenseGeneralMatrix_ImplViennaCLFloat(const TlDenseGeneralMatrix_ImplEigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlDenseGeneralMatrix_ImplViennaCLFloat();

    operator std::vector<double>() const;

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
    TlDenseGeneralMatrix_ImplViennaCLFloat& operator=(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);

    TlDenseGeneralMatrix_ImplViennaCLFloat& operator+=(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ImplViennaCLFloat& operator-=(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ImplViennaCLFloat& operator*=(const float coef);
    TlDenseGeneralMatrix_ImplViennaCLFloat& operator/=(const float coef);
    TlDenseGeneralMatrix_ImplViennaCLFloat& operator*=(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    virtual double sum() const;
    // virtual double trace() const;
    virtual double getRMS() const;
    virtual double getMaxAbsoluteElement(TlMatrixObject::index_type* outRow, TlMatrixObject::index_type* outCol) const;
    virtual void transposeInPlace();

    TlDenseGeneralMatrix_ImplViennaCLFloat& dotInPlace(const TlDenseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlDenseGeneralMatrix_ImplViennaCLFloat transpose() const;
    TlDenseGeneralMatrix_ImplViennaCLFloat inverse() const;

    TlDenseGeneralMatrix_ImplViennaCLFloat& reverseColumns();

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
protected:
    void vtr2mat(const double* pBuf);

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    MatrixDataType matrix_;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ImplViennaCLFloat;
    friend class TlDenseGeneralMatrix_ImplEigenFloat;

    // DM(G) = DM(G) * DM(S)
    friend TlDenseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat1,
        const TlDenseSymmetricMatrix_ImplViennaCLFloat& mat2);
    // DM(G) = DM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlDenseSymmetricMatrix_ImplViennaCLFloat& mat1,
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat2);

    // DV = DM(G) * DV
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vec);
    // DV = DV * DM(G)
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseVector_ImplViennaCLFloat& vec,
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat);
};

// DM(G) = float * DM(G)
TlDenseGeneralMatrix_ImplViennaCLFloat operator*(
    const float coef, const TlDenseGeneralMatrix_ImplViennaCLFloat& dm);

#endif  // TL_DENSE_GENERAL_MATRIX_IMPL_VIENNACL_FLOAT_H
