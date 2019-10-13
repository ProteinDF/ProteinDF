#ifndef TL_DENSE_VECTOR_IMPL_VIENNACL_H
#define TL_DENSE_VECTOR_IMPL_VIENNACL_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include <vector>

#include "tl_dense_vector_impl_object.h"
#include "viennacl/vector.hpp"

class TlDenseGeneralMatrix_ImplViennaCL;
class TlSparseGeneralMatrix_ImplViennaCL;
class TlSparseSymmetricMatrix_ImplViennaCL;
class TlDenseVector_ImplEigen;

class TlDenseVector_ImplViennaCL : public TlDenseVector_ImplObject {
   public:
    typedef viennacl::vector<double> VectorDataType;

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    TlDenseVector_ImplViennaCL(TlDenseVectorObject::index_type dim = 0);
    TlDenseVector_ImplViennaCL(const TlDenseVector_ImplViennaCL& rhs);
    TlDenseVector_ImplViennaCL(const std::vector<double>& rhs);
#ifdef HAVE_EIGEN
    TlDenseVector_ImplViennaCL(const TlDenseVector_ImplEigen& rhs);
#endif  // HAVE_EIGEN

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ImplViennaCL();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
   public:
    virtual TlDenseVectorObject::size_type getSize() const;
    virtual void resize(TlDenseVectorObject::index_type newSize);

   public:
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
    TlDenseVector_ImplViennaCL& operator=(
        const TlDenseVector_ImplViennaCL& rhs);

    TlDenseVector_ImplViennaCL& operator+=(
        const TlDenseVector_ImplViennaCL& rhs);
    TlDenseVector_ImplViennaCL& operator-=(
        const TlDenseVector_ImplViennaCL& rhs);
    TlDenseVector_ImplViennaCL& operator*=(const double rhs);
    TlDenseVector_ImplViennaCL& operator/=(const double rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
   public:
    virtual double sum() const;
    virtual void sortByGreater();
    TlDenseVector_ImplViennaCL& dotInPlace(
        const TlDenseVector_ImplViennaCL& rhs);

    TlDenseVector_ImplViennaCL& reverse();

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
   protected:
    VectorDataType vector_;

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ImplViennaCL;
    friend class TlDenseVector_ImplEigen;

    friend TlDenseVector_ImplViennaCL operator+(
        const TlDenseVector_ImplViennaCL& rhs1,
        const TlDenseVector_ImplViennaCL& rhs2);
    friend TlDenseVector_ImplViennaCL operator-(
        const TlDenseVector_ImplViennaCL& rhs1,
        const TlDenseVector_ImplViennaCL& rhs2);

    friend TlDenseVector_ImplViennaCL operator*(
        const TlDenseGeneralMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vec);
    friend TlDenseVector_ImplViennaCL operator*(
        const TlDenseVector_ImplViennaCL& rhs1, const double rhs2);
    friend TlDenseVector_ImplViennaCL operator*(
        const double rhs1, const TlDenseVector_ImplViennaCL& rhs2);

    friend TlDenseVector_ImplViennaCL operator*(
        const TlDenseGeneralMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vec);
    friend TlDenseVector_ImplViennaCL operator*(
        const TlDenseVector_ImplViennaCL& vec,
        const TlDenseGeneralMatrix_ImplViennaCL& mat);

    // SM(G) x DV
    friend TlDenseVector_ImplViennaCL operator*(
        const TlSparseGeneralMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vtr);
    // DV x SM(G)
    friend TlDenseVector_ImplViennaCL operator*(
        const TlDenseVector_ImplViennaCL& vtr,
        const TlSparseGeneralMatrix_ImplViennaCL& mat);

    friend TlDenseVector_ImplViennaCL operator*(
        const TlSparseSymmetricMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vtr);
};

#endif  // TL_DENSE_VECTOR_IMPL_VIENNACL_H
