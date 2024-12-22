#ifndef TL_DENSE_VECTOR_IMPL_SCALAPACK_H
#define TL_DENSE_VECTOR_IMPL_SCALAPACK_H

#include <valarray>
#include <vector>

#include "tl_dense_scalapack_object.h"
#include "tl_dense_vector_impl_object.h"

class TlDenseGeneralMatrix_ImplScalapack;
class TlDenseVector_ImplLapack;
class TlDenseVector_Lapack;

class TlDenseVector_ImplScalapack : public TlDenseVector_ImplObject,
                                    public TlDenseScalapackObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseVector_ImplScalapack(TlDenseVectorObject::index_type size = 1);
    TlDenseVector_ImplScalapack(const TlDenseVector_ImplScalapack& rhs);
    TlDenseVector_ImplScalapack(const TlDenseVector_ImplLapack& rhs);
    TlDenseVector_ImplScalapack(const std::vector<double>& rhs);

    operator std::vector<double>() const;

    virtual ~TlDenseVector_ImplScalapack();

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

    std::vector<double> getVector() const;

private:
    using TlDenseScalapackObject::add;
    using TlDenseScalapackObject::get;
    using TlDenseScalapackObject::mul;
    using TlDenseScalapackObject::resize;
    using TlDenseScalapackObject::set;

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlDenseVector_ImplScalapack& operator=(
        const TlDenseVector_ImplScalapack& rhs);

    TlDenseVector_ImplScalapack& operator+=(
        const TlDenseVector_ImplScalapack& rhs);
    TlDenseVector_ImplScalapack& operator-=(
        const TlDenseVector_ImplScalapack& rhs);
    TlDenseVector_ImplScalapack& operator*=(const double rhs);
    TlDenseVector_ImplScalapack& operator/=(const double rhs);

    double operator*(const TlDenseVector_ImplScalapack& rhs) const;

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    virtual void sortByGreater();

    double dot(const TlDenseVector_ImplScalapack& rhs) const;
    TlDenseVector_ImplScalapack& dotInPlace(
        const TlDenseVector_ImplScalapack& rhs);

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
public:
    bool load(const std::string& sFilePath);
    bool save(const std::string& sFilePath) const;

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
protected:
    bool load(std::ifstream& ifs);
    //   bool save(std::ofstream& ofs) const;
    std::vector<TlVectorObject::VectorElement> getVectorElementsInLocal() const;

    void saveElements(
        TlDenseVector_Lapack* pVector,
        const std::vector<TlVectorObject::VectorElement>& elements) const;

protected:
    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // friends
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_ImplScalapack;
    friend class TlDenseSymmetricMatrix_ImplScalapack;

    friend TlDenseVector_ImplScalapack operator*(
        const TlDenseGeneralMatrix_ImplScalapack& A,
        const TlDenseVector_ImplScalapack& X);
    friend TlDenseVector_ImplScalapack operator*(
        const TlDenseVector_ImplScalapack& X,
        const TlDenseGeneralMatrix_ImplScalapack& A);
};

#endif  // TL_DENSE_VECTOR_IMPLE_SCALAPACK_H
