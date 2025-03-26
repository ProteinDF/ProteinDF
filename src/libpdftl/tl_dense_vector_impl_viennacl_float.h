#ifndef TL_DENSE_VECTOR_IMPL_VIENNACL_FLOAT_H
#define TL_DENSE_VECTOR_IMPL_VIENNACL_FLOAT_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include <vector>

#include "tl_dense_vector_impl_object.h"
#include "viennacl/vector.hpp"

class TlDenseGeneralMatrix_ImplViennaCLFloat;
class TlSparseGeneralMatrix_ImplViennaCLFloat;
class TlSparseSymmetricMatrix_ImplViennaCLFloat;
class TlDenseVector_ImplEigenFloat;

class TlDenseVector_ImplViennaCLFloat : public TlDenseVector_ImplObject {
public:
    typedef viennacl::vector<float> VectorDataType;

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_ImplViennaCLFloat(TlDenseVectorObject::index_type dim = 0);
    TlDenseVector_ImplViennaCLFloat(const TlDenseVector_ImplViennaCLFloat& rhs);
    TlDenseVector_ImplViennaCLFloat(const std::vector<double>& rhs);
#ifdef HAVE_EIGEN
    TlDenseVector_ImplViennaCLFloat(const TlDenseVector_ImplEigenFloat& rhs);
#endif  // HAVE_EIGEN

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ImplViennaCLFloat();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    virtual TlDenseVectorObject::size_type getSize() const;
    virtual void resize(TlDenseVectorObject::index_type newSize);

public:
    virtual double get(const TlDenseVectorObject::index_type i) const;
    virtual void set(const TlDenseVectorObject::index_type i, const double value);
    virtual void add(const TlDenseVectorObject::index_type i, const double value);
    virtual void mul(const TlDenseVectorObject::index_type i, const double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_ImplViennaCLFloat& operator=(
        const TlDenseVector_ImplViennaCLFloat& rhs);

    TlDenseVector_ImplViennaCLFloat& operator+=(
        const TlDenseVector_ImplViennaCLFloat& rhs);
    TlDenseVector_ImplViennaCLFloat& operator-=(
        const TlDenseVector_ImplViennaCLFloat& rhs);
    TlDenseVector_ImplViennaCLFloat& operator*=(const float rhs);
    TlDenseVector_ImplViennaCLFloat& operator/=(const float rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    virtual double sum() const;
    virtual void sortByGreater();
    TlDenseVector_ImplViennaCLFloat& dotInPlace(const TlDenseVector_ImplViennaCLFloat& rhs);

    TlDenseVector_ImplViennaCLFloat& reverse();

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
protected:
    VectorDataType vector_;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ImplViennaCLFloat;
    friend class TlDenseVector_ImplEigenFloat;

    friend TlDenseVector_ImplViennaCLFloat operator+(
        const TlDenseVector_ImplViennaCLFloat& rhs1,
        const TlDenseVector_ImplViennaCLFloat& rhs2);
    friend TlDenseVector_ImplViennaCLFloat operator-(
        const TlDenseVector_ImplViennaCLFloat& rhs1,
        const TlDenseVector_ImplViennaCLFloat& rhs2);

    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vec);
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseVector_ImplViennaCLFloat& rhs1, const float rhs2);
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const float rhs1, const TlDenseVector_ImplViennaCLFloat& rhs2);

    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vec);
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseVector_ImplViennaCLFloat& vec,
        const TlDenseGeneralMatrix_ImplViennaCLFloat& mat);

    // SM(G) x DV
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vtr);
    // DV x SM(G)
    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlDenseVector_ImplViennaCLFloat& vtr,
        const TlSparseGeneralMatrix_ImplViennaCLFloat& mat);

    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vtr);
};

#endif  // TL_DENSE_VECTOR_IMPL_VIENNACL_FLOAT_H
