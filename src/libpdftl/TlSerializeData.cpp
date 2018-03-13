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

#include <cassert>
#include <cmath>
#include <limits>

#include "TlSerializeData.h"
#include "TlUtils.h"

TlSerializeData* TlSerializeData::pNullObject_ = NULL;

const TlSerializeData& TlSerializeData::getNullObject() {
  if (TlSerializeData::pNullObject_ == NULL) {
    TlSerializeData::pNullObject_ = new TlSerializeData;
  }

  return *(TlSerializeData::pNullObject_);
}

TlSerializeData::TlSerializeData(DataType dataType)
    : type_(dataType), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();
}

TlSerializeData::TlSerializeData(const bool value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const char value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const unsigned char value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const int value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const unsigned int value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const long value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const unsigned long value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const double value)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(value);
}

TlSerializeData::TlSerializeData(const char* pStr)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(std::string(pStr));
}

TlSerializeData::TlSerializeData(const char* pStr, const std::size_t size)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(std::string(pStr, size));
}

TlSerializeData::TlSerializeData(const std::string& str)
    : type_(NONE), scalar_(0), str_("") {
  this->array_.clear();
  this->map_.clear();

  this->set(str);
}

TlSerializeData::TlSerializeData(const TlSerializeData& rhs)
    : type_(rhs.type_), scalar_(rhs.scalar_), str_(rhs.str_) {
  this->array_.clear();
  this->map_.clear();
  this->copyChildren(rhs);
}

TlSerializeData::~TlSerializeData() { this->clearChildren(); }

TlSerializeData& TlSerializeData::operator=(const TlSerializeData& rhs) {
  if (this != &rhs) {
    this->clearChildren();

    this->type_ = rhs.type_;
    this->scalar_ = rhs.scalar_;
    this->str_ = rhs.str_;

    this->copyChildren(rhs);
  }

  return *this;
}

bool TlSerializeData::operator==(const TlSerializeData& rhs) const {
  bool answer = false;
  if (this->getType() == rhs.getType()) {
    switch (this->getType()) {
      case BOOLEAN:
      case INT:
        answer = (this->scalar_.int_ == rhs.scalar_.int_);
        break;

      case UINT:
        answer = (this->scalar_.uint_ == rhs.scalar_.uint_);
        break;

      case LONG:
        answer = (this->scalar_.long_ == rhs.scalar_.long_);
        break;

      case ULONG:
        answer = (this->scalar_.ulong_ == rhs.scalar_.ulong_);
        break;

      case DOUBLE:
        answer = (std::fabs(this->scalar_.double_ - rhs.scalar_.double_) <
                  std::numeric_limits<double>::epsilon());
        break;

      case STRING:
        answer = (this->str_ == rhs.str_);
        break;

      case ARRAY:
        if (this->getSize() == rhs.getSize()) {
          answer = true;
          const std::size_t maxSize = this->getSize();
          for (std::size_t i = 0; i < maxSize; ++i) {
            if (this->getAt(i) != rhs.getAt(i)) {
              answer = false;
              break;
            }
          }
        }
        break;

      case MAP:
        if (this->getSize() == rhs.getSize()) {
        }
        break;

      default:
        answer = true;
    }
  }

  return answer;
}

void TlSerializeData::clearChildren() {
  for (ArrayContainerType::iterator p = this->array_.begin();
       p != this->array_.end(); ++p) {
    delete *p;
    *p = NULL;
  }
  this->array_.clear();

  for (MapContainerType::iterator p = this->map_.begin(); p != this->map_.end();
       ++p) {
    delete p->first;
    // p->first = NULL;

    delete p->second;
    p->second = NULL;
  }
  this->map_.clear();
}

void TlSerializeData::copyChildren(const TlSerializeData& rhs) {
  assert(this->array_.size() == 0);
  if (rhs.array_.empty() != true) {
    this->array_.reserve(rhs.array_.size());
    for (ArrayContainerType::const_iterator p = rhs.array_.begin();
         p != rhs.array_.end(); ++p) {
      TlSerializeData* pNew = new TlSerializeData(*(*p));
      this->array_.push_back(pNew);
    }
  }

  assert(this->map_.size() == 0);
  if (rhs.map_.empty() != true) {
    for (MapContainerType::const_iterator p = rhs.map_.begin();
         p != rhs.map_.end(); ++p) {
      TlSerializeData* pKey = new TlSerializeData(*(p->first));
      TlSerializeData* pValue = new TlSerializeData(*(p->second));
      this->map_.insert(
          std::pair<TlSerializeData*, TlSerializeData*>(pKey, pValue));
    }
  }
}

TlSerializeData::DataType TlSerializeData::getType() const {
  return this->type_;
}

std::size_t TlSerializeData::getSize() const {
  std::size_t ans = 1;
  if (this->type_ == ARRAY) {
    ans = this->array_.size();
  } else if (this->type_ == MAP) {
    ans = this->map_.size();
  }

  return ans;
}

void TlSerializeData::set(const bool value) {
  this->type_ = BOOLEAN;
  this->scalar_.int_ = (value == true) ? 1 : 0;
}

void TlSerializeData::set(const char value) {
  this->type_ = INT;
  this->scalar_.int_ = value;
}

void TlSerializeData::set(const unsigned char value) {
  this->type_ = UINT;
  this->scalar_.uint_ = value;
}

void TlSerializeData::set(const int value) {
  this->type_ = INT;
  this->scalar_.int_ = value;
}

void TlSerializeData::set(const unsigned int value) {
  this->type_ = UINT;
  this->scalar_.uint_ = value;
}

void TlSerializeData::set(const long value) {
  this->type_ = LONG;
  this->scalar_.long_ = value;
}

void TlSerializeData::set(const unsigned long value) {
  this->type_ = ULONG;
  this->scalar_.ulong_ = value;
}

void TlSerializeData::set(const double value) {
  this->type_ = DOUBLE;
  this->scalar_.double_ = value;
}

void TlSerializeData::set(const char* pStr) { this->set(std::string(pStr)); }

void TlSerializeData::set(const char* pStr, const std::size_t size) {
  this->set(std::string(pStr, size));
}

void TlSerializeData::set(const std::string& value) {
  this->type_ = STRING;
  this->str_ = value;
}

// =====================================================================
void TlSerializeData::resize(const std::size_t newSize) {
  this->type_ = ARRAY;

  const std::size_t oldSize = this->array_.size();
  if (newSize < oldSize) {
    // shrink
    for (std::size_t i = newSize; i < oldSize; ++i) {
      delete this->array_[i];
      this->array_[i] = NULL;
    }
    this->array_.resize(newSize);
  } else if (newSize > oldSize) {
    // expand
    this->array_.resize(newSize);
    for (std::size_t i = oldSize; i < newSize; ++i) {
      TlSerializeData* p = new TlSerializeData;
      this->array_[i] = p;
    }
  }
}

void TlSerializeData::pushBack(const TlSerializeData& value) {
  this->type_ = ARRAY;

  TlSerializeData* pNew = new TlSerializeData(value);
  this->array_.push_back(pNew);
}

const TlSerializeData& TlSerializeData::getAt(const std::size_t index) const {
  if ((this->type_ == ARRAY) && (index < this->array_.size())) {
    return *(this->array_[index]);
  } else {
    return TlSerializeData::getNullObject();
  }
}

TlSerializeData& TlSerializeData::getAt(const std::size_t index) {
  assert(this->type_ == ARRAY);
  if ((index + 1) > this->array_.size()) {
    this->resize(index + 1);
  }

  return *(this->array_[index]);
}

void TlSerializeData::setAt(const std::size_t index,
                            const TlSerializeData& value) {
  this->type_ = ARRAY;
  if ((index + 1) > this->array_.size()) {
    this->resize(index + 1);
  }
  *(this->array_[index]) = value;
}

TlSerializeData::ArrayIterator TlSerializeData::beginArray() {
  return ArrayIterator(this->array_.begin());
}

TlSerializeData::ArrayIterator TlSerializeData::endArray() {
  return ArrayIterator(this->array_.end());
}

TlSerializeData::ArrayConstIterator TlSerializeData::beginArray() const {
  return ArrayConstIterator(this->array_.begin());
}

TlSerializeData::ArrayConstIterator TlSerializeData::endArray() const {
  return ArrayConstIterator(this->array_.end());
}

// =====================================================================

void TlSerializeData::add(const TlSerializeData& key,
                          const TlSerializeData& value) {
  this->type_ = MAP;
  TlSerializeData* pKey = new TlSerializeData(key);
  TlSerializeData* pValue = new TlSerializeData(value);
  this->map_[pKey] = pValue;
}

TlSerializeData& TlSerializeData::operator[](const TlSerializeData& key) {
  this->type_ = MAP;
  MapContainerType::iterator p = this->find(key);
  if (p != this->map_.end()) {
    return *(p->second);
  } else {
    TlSerializeData* pKey = new TlSerializeData(key);
    TlSerializeData* pValue = new TlSerializeData;
    this->map_[pKey] = pValue;
    return *(pValue);
  }
}

const TlSerializeData& TlSerializeData::operator[](
    const TlSerializeData& key) const {
  if (this->getType() == MAP) {
    MapContainerType::const_iterator p = this->find(key);
    if (p != this->map_.end()) {
      return *(p->second);
    }
  }

  return this->getNullObject();
}

bool TlSerializeData::hasKey(const TlSerializeData& key) const {
  bool answer = false;
  if (this->getType() == MAP) {
    MapContainerType::const_iterator p = this->find(key);
    if (p != this->map_.end()) {
      answer = true;
    }
  }
  return answer;
}

TlSerializeData::MapContainerType::iterator TlSerializeData::find(
    const TlSerializeData& key) {
  MapContainerType::iterator answer = this->map_.end();
  if (this->map_.size() > 0) {
    MapContainerType::iterator pEnd = this->map_.end();
    for (MapContainerType::iterator p = this->map_.begin(); p != pEnd; ++p) {
      const TlSerializeData tmp = *(p->first);
      if (tmp == key) {
        answer = p;
        break;
      }
    }
  }

  return answer;
}

TlSerializeData::MapContainerType::const_iterator TlSerializeData::find(
    const TlSerializeData& key) const {
  MapContainerType::const_iterator answer = this->map_.end();

  MapContainerType::const_iterator pEnd = this->map_.end();
  for (MapContainerType::const_iterator p = this->map_.begin(); p != pEnd;
       ++p) {
    if (*(p->first) == key) {
      answer = p;
      break;
    }
  }

  return answer;
}

void TlSerializeData::erase(const TlSerializeData& key) {
  MapContainerType::iterator p = this->map_.begin();
  MapContainerType::iterator pEnd = this->map_.end();
  while (p != pEnd) {
    if (*(p->first) == key) {
      delete p->second;
      p->second = NULL;

      delete p->first;
      this->map_.erase(p++);
      break;
    } else {
      ++p;
    }
  }
}

TlSerializeData::MapIterator TlSerializeData::beginMap() {
  return MapIterator(this->map_.begin());
}

TlSerializeData::MapIterator TlSerializeData::endMap() {
  return MapIterator(this->map_.end());
}

TlSerializeData::MapConstIterator TlSerializeData::beginMap() const {
  return MapConstIterator(this->map_.begin());
}

TlSerializeData::MapConstIterator TlSerializeData::endMap() const {
  return MapConstIterator(this->map_.end());
}

// =====================================================================

bool TlSerializeData::getBoolean() const {
  bool answer = false;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0);
      break;

    case STRING: {
      const std::string check = TlUtils::toUpper(this->str_);
      if ((check == "TRUE") || (check == "YES") || (check == "ON") ||
          (check == "1")) {
        answer = true;
      }
    } break;

    case INT:
      answer = (this->scalar_.int_ != 0);
      break;

    case UINT:
      answer = (this->scalar_.uint_ != 0);
      break;

    case LONG:
      answer = (this->scalar_.long_ != 0);
      break;

    case ULONG:
      answer = (this->scalar_.ulong_ != 0);
      break;

    case DOUBLE:
      answer = (std::fabs(this->scalar_.double_) >
                std::numeric_limits<double>::epsilon());
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

int TlSerializeData::getInt() const {
  int answer = 0;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? 1 : 0;
      break;

    case STRING:
      answer = std::atoi(this->str_.c_str());
      break;

    //     case CHAR:
    //         answer = this->scalar_.char_;
    //         break;

    //     case UCHAR:
    //         answer = this->scalar_.uchar_;
    //         break;

    case INT:
      answer = this->scalar_.int_;
      break;

    case UINT:
      answer = (int)this->scalar_.uint_;
      break;

    case LONG:
      answer = (int)this->scalar_.long_;
      break;

    case ULONG:
      answer = (int)this->scalar_.ulong_;
      break;

    case DOUBLE:
      answer = (int)this->scalar_.double_;
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

unsigned int TlSerializeData::getUInt() const {
  unsigned int answer = 0;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? 1 : 0;
      break;

    case STRING:
      answer = std::atoi(this->str_.c_str());
      break;

    //     case CHAR:
    //         answer = (unsigned int)this->scalar_.char_;
    //         break;

    //     case UCHAR:
    //         answer = this->scalar_.uchar_;
    //         break;

    case INT:
      answer = (unsigned int)this->scalar_.int_;
      break;

    case UINT:
      answer = this->scalar_.uint_;
      break;

    case LONG:
      answer = (unsigned int)this->scalar_.long_;
      break;

    case ULONG:
      answer = (unsigned int)this->scalar_.ulong_;
      break;

    case DOUBLE:
      answer = (unsigned int)this->scalar_.double_;
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

long TlSerializeData::getLong() const {
  long answer = 0;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? 1 : 0;
      break;

    case STRING:
      answer = std::atol(this->str_.c_str());
      break;

    //     case CHAR:
    //         answer = this->scalar_.char_;
    //         break;

    //     case UCHAR:
    //         answer = this->scalar_.uchar_;
    //         break;

    case INT:
      answer = (long)this->scalar_.int_;
      break;

    case UINT:
      answer = (long)this->scalar_.uint_;
      break;

    case LONG:
      answer = this->scalar_.long_;
      break;

    case ULONG:
      answer = (long)this->scalar_.ulong_;
      break;

    case DOUBLE:
      answer = (long)this->scalar_.double_;
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

unsigned long TlSerializeData::getULong() const {
  unsigned long answer = 0;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? 1 : 0;
      break;

    case STRING:
      answer = std::atol(this->str_.c_str());
      break;

    //     case CHAR:
    //         answer = (unsigned long)this->scalar_.char_;
    //         break;

    //     case UCHAR:
    //         answer = this->scalar_.uchar_;
    //         break;

    case INT:
      answer = (unsigned long)this->scalar_.int_;
      break;

    case UINT:
      answer = (unsigned long)this->scalar_.uint_;
      break;

    case LONG:
      answer = (unsigned long)this->scalar_.long_;
      break;

    case ULONG:
      answer = this->scalar_.ulong_;
      break;

    case DOUBLE:
      answer = (unsigned long)this->scalar_.double_;
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

double TlSerializeData::getDouble() const {
  double answer = 0.0;
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? 1.0 : 0.0;
      break;

    case STRING:
      answer = std::atof(this->str_.c_str());
      break;

    //     case CHAR:
    //         answer = (double)this->scalar_.char_;
    //         break;

    //     case UCHAR:
    //         answer = (double)this->scalar_.uchar_;
    //         break;

    case INT:
      answer = (double)this->scalar_.int_;
      break;

    case UINT:
      answer = (double)this->scalar_.uint_;
      break;

    case LONG:
      answer = (double)this->scalar_.long_;
      break;

    case ULONG:
      answer = (double)this->scalar_.ulong_;
      break;

    case DOUBLE:
      answer = this->scalar_.double_;
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

std::string TlSerializeData::getStr() const {
  std::string answer = "";
  switch (this->getType()) {
    case BOOLEAN:
      answer = (this->scalar_.int_ != 0) ? "true" : "false";
      break;

    case STRING:
      answer = this->str_;
      break;

    //     case CHAR:
    //         answer = TlUtils::xtos(this->scalar_.char_);
    //         break;

    //     case UCHAR:
    //         answer = TlUtils::xtos(this->scalar_.uchar_);
    //         break;

    case INT:
      answer = TlUtils::xtos(this->scalar_.int_);
      break;

    case UINT:
      answer = TlUtils::xtos(this->scalar_.uint_);
      break;

    case LONG:
      answer = TlUtils::xtos(this->scalar_.long_);
      break;

    case ULONG:
      answer = TlUtils::xtos(this->scalar_.ulong_);
      break;

    case DOUBLE:
      answer = TlUtils::xtos(this->scalar_.double_);
      break;

    default:
      // do nothing
      break;
  }

  return answer;
}

void TlSerializeData::merge(const TlSerializeData& rhs) {
  const DataType dataType = rhs.getType();
  //     if (dataType == ARRAY) {
  //         ArrayConstIterator pEnd = rhs.endArray();
  //         for (ArrayConstIterator p = rhs.beginArray(); p != pEnd; ++p) {
  //             this->pushBack(*p);
  //         }
  //     } else
  if (dataType == MAP) {
    MapConstIterator pEnd = rhs.endMap();
    for (MapConstIterator p = rhs.beginMap(); p != pEnd; ++p) {
      TlSerializeData key = p->first;
      TlSerializeData value = p->second;
      (*this)[key].merge(value);
    }
  } else {
    this->operator=(rhs);
  }
}

std::string TlSerializeData::str() const {
  std::string ans = "";

  switch (this->getType()) {
    case ARRAY: {
      ans += "[";
      ArrayConstIterator pEnd = this->endArray();
      for (ArrayConstIterator p = this->beginArray(); p != pEnd; ++p) {
        const std::string item = p->str();
        ans += TlUtils::format("\"%s\", ", item.c_str());
      }
      ans = ans.substr(0, ans.size() - 2);  // remove ", ".
      ans += "]";
    } break;

    case MAP: {
      ans += "{";
      MapConstIterator pEnd = this->endMap();
      for (MapConstIterator p = this->beginMap(); p != pEnd; ++p) {
        const std::string key = p->first.str();
        const std::string value = p->second.str();
        ans += TlUtils::format("\"%s\": \"%s\", ", key.c_str(), value.c_str());
      }
      ans = ans.substr(0, ans.size() - 2);  // remove ", ".
      ans += "}";
    } break;

    default: { ans = this->getStr(); } break;
  }

  return ans;
}
