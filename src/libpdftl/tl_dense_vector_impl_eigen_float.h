#ifndef TL_DENSE_VECTOR_IMPL_EIGEN_FLOAT_H
#define TL_DENSE_VECTOR_IMPL_EIGEN_FLOAT_H

#include <Eigen/Core>
#include <vector>

#include "tl_dense_vector_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlDenseGeneralMatrix_ImplEigenFloat;
class TlDenseSymmetricMatrix_ImplEigenFloat;
class TlSparseGeneralMatrix_ImplEigenFloat;
class TlSparseSymmetricMatrix_ImplEigenFloat;
class TlDenseVector_ImplViennaCL;

class TlDenseVector_ImplEigenFloat : public TlDenseVector_ImplObject {
public:
    typedef Eigen::VectorXf VectorDataType;  // Matrix<float, Dynamic, 1>
    typedef Eigen::Map<VectorDataType> MapType;
    typedef Eigen::Map<const VectorDataType> MapTypeConst;

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseVector_ImplEigenFloat(const TlDenseVectorObject::index_type size = 0);
    TlDenseVector_ImplEigenFloat(const TlDenseVector_ImplEigenFloat& rhs);
    TlDenseVector_ImplEigenFloat(const std::vector<double>& rhs);
    TlDenseVector_ImplEigenFloat(const double* p,
                            const TlDenseVectorObject::index_type size);
    TlDenseVector_ImplEigenFloat(const VectorDataType& rhs);

#ifdef HAVE_VIENNACL
    TlDenseVector_ImplEigenFloat(const TlDenseVector_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ImplEigenFloat();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    virtual TlDenseVectorObject::size_type getSize() const;
    virtual void resize(const TlDenseVectorObject::index_type newSize);

    virtual double get(const TlDenseVectorObject::index_type i) const;
    virtual void set(const TlDenseVectorObject::index_type i,
                     const double value);
    virtual void add(const TlDenseVectorObject::index_type i,
                     const double value);
    virtual void mul(const TlDenseVectorObject::index_type i,
                     const double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_ImplEigenFloat& operator=(const TlDenseVector_ImplEigenFloat& rhs);

    TlDenseVector_ImplEigenFloat& operator+=(const TlDenseVector_ImplEigenFloat& rhs);
    TlDenseVector_ImplEigenFloat& operator-=(const TlDenseVector_ImplEigenFloat& rhs);
    TlDenseVector_ImplEigenFloat& operator*=(const double rhs);
    TlDenseVector_ImplEigenFloat& operator/=(const double rhs);

    double operator*(const TlDenseVector_ImplEigenFloat& rhs) const;

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    virtual double sum() const;
    virtual void sortByGreater();

    double dot(const TlDenseVector_ImplEigenFloat& rhs) const;
    TlDenseVector_ImplEigenFloat& dotInPlace(const TlDenseVector_ImplEigenFloat& rhs);

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    VectorDataType vector_;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseVector_ImplViennaCL;

    friend TlDenseVector_ImplEigenFloat operator+(
        const TlDenseVector_ImplEigenFloat& rhs1,
        const TlDenseVector_ImplEigenFloat& rhs2);
    friend TlDenseVector_ImplEigenFloat operator-(
        const TlDenseVector_ImplEigenFloat& rhs1,
        const TlDenseVector_ImplEigenFloat& rhs2);

    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& rhs1, const double rhs2);
    friend TlDenseVector_ImplEigenFloat operator*(
        const double rhs1, const TlDenseVector_ImplEigenFloat& rhs2);

    // DM(G) * DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseGeneralMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vec);
    // DV * DM(G)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& vec,
        const TlDenseGeneralMatrix_ImplEigenFloat& mat);
    // DM(S) * DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vec);
    // DV * DM(S)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& vec,
        const TlDenseSymmetricMatrix_ImplEigenFloat& mat);

    // SM(G) x DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlSparseGeneralMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vtr);
    // DV x SM(G)
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlDenseVector_ImplEigenFloat& vtr,
        const TlSparseGeneralMatrix_ImplEigenFloat& mat);
    // SM(S) x DV
    friend TlDenseVector_ImplEigenFloat operator*(
        const TlSparseSymmetricMatrix_ImplEigenFloat& mat,
        const TlDenseVector_ImplEigenFloat& vtr);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_DENSE_VECTOR_IMPL_EIGEN_FLOAT_H
