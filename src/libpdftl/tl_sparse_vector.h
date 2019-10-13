// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef TLSPARSEVECTOR_H
#define TLSPARSEVECTOR_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <map>
#include "TlStlUtils.h"
#include "TlUtils.h"
#include "tl_vector_abstract.h"

/// 疎ベクトルクラス
// class TlSparseVector : public TlVectorAbstract {
class TlSparseVector {
   public:
    typedef TlVectorAbstract::size_type size_type;

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
    const_iterator begin() const { return this->data_.begin(); }

    /// 最初のイテレータを返す
    iterator begin() { return this->data_.begin(); }

    /// 最後のイテレータを返す
    const_iterator end() const { return this->data_.end(); }

    /// 最後のイテレータを返す
    iterator end() { return this->data_.end(); }

   protected:
    size_type size_;
    mutable SparseVectorData data_;
};

#endif  // TLSPARSEVECTOR_H
