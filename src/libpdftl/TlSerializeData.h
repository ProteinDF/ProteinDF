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

#ifndef TLSERIALIZEDATA_H
#define TLSERIALIZEDATA_H

#include <map>
#include <string>
#include <vector>

template <typename IteratorType, typename ValueType>
class TlSerializeDataVectorIterator {
 public:
  TlSerializeDataVectorIterator(const IteratorType& it) : it_(it) {}

  TlSerializeDataVectorIterator(const TlSerializeDataVectorIterator& rhs)
      : it_(rhs.it_) {}

  ~TlSerializeDataVectorIterator() {}

  TlSerializeDataVectorIterator& operator=(
      const TlSerializeDataVectorIterator& rhs) {
    if (&rhs != this) {
      this->it_ = rhs.it_;
    }
  }

 public:
  ValueType& operator*() { return **(this->it_); }

  ValueType* operator->() { return *(this->it_); }

  TlSerializeDataVectorIterator& operator++() {
    ++(this->it_);
    return *this;
  }

  TlSerializeDataVectorIterator& operator++(int n) {
    this->it_++(n);
    return *this;
  }

  TlSerializeDataVectorIterator& operator--() {
    --(this->it_);
    return *this;
  }

  TlSerializeDataVectorIterator& operator--(int n) {
    this->it_--(n);
    return *this;
  }

  bool operator==(const TlSerializeDataVectorIterator& rhs) const {
    return (this->it_ == rhs.it_);
  }

  bool operator!=(const TlSerializeDataVectorIterator& rhs) const {
    return (this->it_ != rhs.it_);
  }

 private:
  IteratorType it_;
};

template <typename IteratorType, typename ValueType>
class TlSerializeDataVectorConstIterator {
 public:
  TlSerializeDataVectorConstIterator(const IteratorType& it) : it_(it) {}

  TlSerializeDataVectorConstIterator(
      const TlSerializeDataVectorConstIterator& rhs)
      : it_(rhs.it_) {}

  ~TlSerializeDataVectorConstIterator() {}

  TlSerializeDataVectorConstIterator& operator=(
      const TlSerializeDataVectorConstIterator& rhs) {
    if (&rhs != this) {
      this->it_ = rhs.it_;
    }
  }

 public:
  ValueType& operator*() const { return **(this->it_); }

  ValueType* operator->() const { return *(this->it_); }

  TlSerializeDataVectorConstIterator& operator++() {
    ++(this->it_);
    return *this;
  }

  TlSerializeDataVectorConstIterator& operator++(int n) {
    this->it_++(n);
    return *this;
  }

  TlSerializeDataVectorConstIterator& operator--() {
    --(this->it_);
    return *this;
  }

  TlSerializeDataVectorConstIterator& operator--(int n) {
    this->it_--(n);
    return *this;
  }

  bool operator==(const TlSerializeDataVectorConstIterator& rhs) const {
    return (this->it_ == rhs.it_);
  }

  bool operator!=(const TlSerializeDataVectorConstIterator& rhs) const {
    return (this->it_ != rhs.it_);
  }

 private:
  IteratorType it_;
};

template <typename IteratorType, typename KeyType, typename ValueType>
class TlSerializeDataMapIterator {
 public:
  typedef std::pair<KeyType, ValueType> PairType;

 public:
  TlSerializeDataMapIterator(const IteratorType& it) : it_(it) {}

  TlSerializeDataMapIterator(const TlSerializeDataMapIterator& rhs)
      : it_(rhs.it_) {}

  ~TlSerializeDataMapIterator() {}

  TlSerializeDataMapIterator& operator=(const TlSerializeDataMapIterator& rhs) {
    if (&rhs != this) {
      this->it_ = rhs.it_;
    }
    return *this;
  }

 public:
  PairType& operator*() {
    this->pair_.first = *(this->it_->first);
    this->pair_.second = *(this->it_->second);
    return this->pair_;
  }

  PairType* operator->() {
    this->pair_.first = *(this->it_->first);
    this->pair_.second = *(this->it_->second);
    return &(this->pair_);
  }

  TlSerializeDataMapIterator& operator++() {
    ++(this->it_);
    return *this;
  }

  TlSerializeDataMapIterator& operator++(int n) {
    this->it_++(n);
    return *this;
  }

  TlSerializeDataMapIterator& operator--() {
    --(this->it_);
    return *this;
  }

  TlSerializeDataMapIterator& operator--(int n) {
    this->it_--(n);
    return *this;
  }

  bool operator==(const TlSerializeDataMapIterator& rhs) const {
    return (this->it_ == rhs.it_);
  }

  bool operator!=(const TlSerializeDataMapIterator& rhs) const {
    return (this->it_ != rhs.it_);
  }

 private:
  IteratorType it_;
  PairType pair_;
};

template <typename IteratorType, typename KeyType, typename ValueType>
class TlSerializeDataMapConstIterator {
 public:
  typedef std::pair<KeyType, ValueType> PairType;

 public:
  TlSerializeDataMapConstIterator(const IteratorType& it) : it_(it) {}

  TlSerializeDataMapConstIterator(const TlSerializeDataMapConstIterator& rhs)
      : it_(rhs.it_) {}

  ~TlSerializeDataMapConstIterator() {}

  TlSerializeDataMapConstIterator& operator=(
      const TlSerializeDataMapConstIterator& rhs) {
    if (&rhs != this) {
      this->it_ = rhs.it_;
    }
  }

 public:
  const PairType& operator*() const {
    this->pair_.first = *(this->it_->first);
    this->pair_.second = *(this->it_->second);
    return this->pair_;
  }

  const PairType* operator->() const {
    this->pair_.first = *(this->it_->first);
    this->pair_.second = *(this->it_->second);
    return &(this->pair_);
  }

  TlSerializeDataMapConstIterator& operator++() {
    ++(this->it_);
    return *this;
  }

  TlSerializeDataMapConstIterator& operator++(int n) {
    this->it_++(n);
    return *this;
  }

  TlSerializeDataMapConstIterator& operator--() {
    --(this->it_);
    return *this;
  }

  TlSerializeDataMapConstIterator& operator--(int n) {
    this->it_--(n);
    return *this;
  }

  bool operator==(const TlSerializeDataMapConstIterator& rhs) const {
    return (this->it_ == rhs.it_);
  }

  bool operator!=(const TlSerializeDataMapConstIterator& rhs) const {
    return (this->it_ != rhs.it_);
  }

 private:
  IteratorType it_;
  mutable PairType pair_;
};

class TlSerializeData {
 protected:
  typedef std::vector<TlSerializeData*> ArrayContainerType;
  typedef std::map<TlSerializeData*, TlSerializeData*> MapContainerType;

 public:
  typedef TlSerializeDataVectorIterator<ArrayContainerType::iterator,
                                        TlSerializeData>
      ArrayIterator;
  typedef TlSerializeDataVectorConstIterator<ArrayContainerType::const_iterator,
                                             TlSerializeData>
      ArrayConstIterator;
  typedef TlSerializeDataMapIterator<MapContainerType::iterator,
                                     TlSerializeData, TlSerializeData>
      MapIterator;
  typedef TlSerializeDataMapConstIterator<MapContainerType::const_iterator,
                                          TlSerializeData, TlSerializeData>
      MapConstIterator;

  enum DataType {
    BOOLEAN,
    STRING,
    INT,
    UINT,
    LONG,
    ULONG,
    DOUBLE,
    ARRAY,
    MAP,
    NONE
  };

 public:
  TlSerializeData(DataType dataType = NONE);
  TlSerializeData(bool value);
  TlSerializeData(char value);
  TlSerializeData(unsigned char value);
  TlSerializeData(int value);
  TlSerializeData(unsigned int value);
  TlSerializeData(long value);
  TlSerializeData(unsigned long value);
  TlSerializeData(double value);
  TlSerializeData(const char* pStr);
  TlSerializeData(const char* pStr, const std::size_t size);
  TlSerializeData(const std::string& str);
  TlSerializeData(const TlSerializeData& rhs);
  virtual ~TlSerializeData();

  TlSerializeData& operator=(const TlSerializeData& rhs);

  bool operator==(const TlSerializeData& rhs) const;
  bool operator!=(const TlSerializeData& rhs) const {
    return !(this->operator==(rhs));
  }

 public:
  // 格納されている型を返す
  DataType getType() const;

  std::size_t getSize() const;

  // list container関連
  void resize(std::size_t size);
  void pushBack(const TlSerializeData& value);
  const TlSerializeData& getAt(std::size_t index) const;
  TlSerializeData& getAt(std::size_t index);
  void setAt(std::size_t index, const TlSerializeData& value);

  ArrayIterator beginArray();
  ArrayIterator endArray();
  ArrayConstIterator beginArray() const;
  ArrayConstIterator endArray() const;

  // map container関連
  void add(const TlSerializeData& key, const TlSerializeData& value);
  TlSerializeData& operator[](const TlSerializeData& key);
  const TlSerializeData& operator[](const TlSerializeData& key) const;
  bool hasKey(const TlSerializeData& key) const;
  void erase(const TlSerializeData& key);

  MapIterator beginMap();
  MapIterator endMap();
  MapConstIterator beginMap() const;
  MapConstIterator endMap() const;

  void set(const bool value);
  void set(const char value);
  void set(const unsigned char value);
  void set(const int value);
  void set(const unsigned int value);
  void set(const long value);
  void set(const unsigned long value);
  void set(const double value);
  void set(const char* pStr);
  void set(const char* pStr, const std::size_t size);
  void set(const std::string& value);

  bool getBoolean() const;
  int getInt() const;
  unsigned int getUInt() const;
  long getLong() const;
  unsigned long getULong() const;
  double getDouble() const;
  std::string getStr() const;

  void merge(const TlSerializeData& rhs);

  /// 内容をデバッグ用の文字列として返す
  std::string str() const;

 protected:
  void clearChildren();
  void copyChildren(const TlSerializeData& rhs);

  MapContainerType::iterator find(const TlSerializeData& rhs);
  MapContainerType::const_iterator find(const TlSerializeData& rhs) const;

  static const TlSerializeData& getNullObject();

 protected:
  union Scalar {
    Scalar(int v) { this->int_ = v; }

    Scalar(unsigned int v) { this->uint_ = v; }

    Scalar(long v) { this->long_ = v; }

    Scalar(unsigned long v) { this->ulong_ = v; }

    Scalar(double v) { this->double_ = v; }

    int int_;
    unsigned int uint_;
    long long_;
    unsigned ulong_;
    double double_;
  };

  DataType type_;
  Scalar scalar_;
  std::string str_;
  // std::string* pStr_;
  std::vector<TlSerializeData*> array_;
  std::map<TlSerializeData*, TlSerializeData*> map_;

  static TlSerializeData* pNullObject_;
};

#endif  // TLSERIALIZEDATA_H
