#ifndef TLSPARSEVECTOR_H
#define TLSPARSEVECTOR_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <map>
#include "TlVectorObject.h"
#include "TlStlUtils.h"
#include "TlUtils.h"

/// 疎ベクトルクラス
//class TlSparseVector : public TlVectorObject {
class TlSparseVector {
public:
    typedef TlVectorObject::size_type size_type;

public:
    explicit TlSparseVector(size_type size = 1);
    virtual ~TlSparseVector();

public:
    virtual size_type getSize() const;
    virtual void resize(size_type newSize);

public:
    virtual double get(const size_type index) const;
    virtual void set(const size_type index, const double value);
    virtual void add(const size_type index, const double value);

    double operator[](const size_type index) const;
    double& operator[](const size_type index);

    TlSparseVector& operator+=(const TlSparseVector& rhs);
    TlSparseVector& operator-=(const TlSparseVector& rhs);

public:
    typedef size_type KeyType;
    typedef std::map<KeyType, double> SparseVectorData;

    typedef SparseVectorData::const_iterator const_iterator;
    typedef SparseVectorData::iterator iterator;

    // operator
public:
    /// 最初のイテレータを返す
    const_iterator begin() const {
        return this->data_.begin();
    }

    /// 最初のイテレータを返す
    iterator begin() {
        return this->data_.begin();
    }

    /// 最後のイテレータを返す
    const_iterator end() const {
        return this->data_.end();
    }

    /// 最後のイテレータを返す
    iterator end() {
        return this->data_.end();
    }

protected:
    size_type size_;
    mutable SparseVectorData data_;
};

#endif // TLSPARSEVECTOR_H
