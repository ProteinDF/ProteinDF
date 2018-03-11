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

#ifndef TLSTLUTILS_H
#define TLSTLUTILS_H

#include <list>
#include <map>
#include "TlSharedPointer.h"

class TlStlUtils {
 public:
  template <typename MapType, typename KeyArgType, typename ValueArgType>
  static typename MapType::iterator efficientAddOrUpdate(MapType& m,
                                                         const KeyArgType& k,
                                                         const ValueArgType& v);
};

// see Effective STL 3-24
template <typename MapType, typename KeyArgType, typename ValueArgType>
typename MapType::iterator TlStlUtils::efficientAddOrUpdate(
    MapType& m, const KeyArgType& k, const ValueArgType& v) {
  typename MapType::iterator lb = m.lower_bound(k);
  if (lb != m.end() && !(m.key_comp()(k, lb->first))) {
    lb->second = v;
    return lb;
  } else {
    typedef typename MapType::value_type MVT;
    return m.insert(lb, MVT(k, v));
  }
}

template <typename KeyType, typename DataType,
          typename KeyCompare = std::less<KeyType> >
class TlCache {
 private:
  // データ本体はCacheTableに保存される。
  // MRUリストはKeyListにキーのリストとして格納される。
  typedef std::list<KeyType> KeyList;
  typedef typename KeyList::iterator KeyListItr;
  typedef typename std::pair<TlSharedPointer<DataType>, KeyListItr> CacheType;
  typedef std::map<KeyType, CacheType, KeyCompare> CacheTable;

 public:
  // デフォルトのキャッシュ個数は10,000個
  TlCache(const std::size_t maxItems = 10000)
      : maxItems_(maxItems), numOfItems_(0) {}
  ~TlCache() {}

 public:
  void set(const KeyType& key, const DataType& data) {
    // 古いデータを消す
    while (this->getNumOfItems() >= this->getMaxItems()) {
      KeyListItr reqHistoryTail = this->reqHistory_.end();
      --reqHistoryTail;
      this->cacheTbl_.erase(*reqHistoryTail);
      this->reqHistory_.erase(reqHistoryTail);
      --(this->numOfItems_);
    }

    typename CacheTable::iterator pCache = this->cacheTbl_.lower_bound(key);
    if ((pCache != this->cacheTbl_.end()) &&
        (this->cacheTbl_.key_comp()(key, pCache->first) != true)) {
      pCache->second.first = TlSharedPointer<DataType>(new DataType(data));
    } else {
      // create new data
      CacheType cacheData = std::make_pair(
          TlSharedPointer<DataType>(new DataType(data)),
          this->reqHistory_.insert(this->reqHistory_.begin(), key));
      this->cacheTbl_.insert(pCache, std::make_pair(key, cacheData));
      ++(this->numOfItems_);
    }
  }

  TlSharedPointer<DataType> get(const KeyType& key) {
    KeyListItr reqHistoryHead = this->reqHistory_.begin();
    typename CacheTable::iterator pCache = this->cacheTbl_.lower_bound(key);
    if ((pCache != this->cacheTbl_.end()) &&
        (this->cacheTbl_.key_comp()(key, pCache->first) != true)) {
      // cache hit!
      this->reqHistory_.splice(reqHistoryHead, this->reqHistory_,
                               pCache->second.second);
      return pCache->second.first;
    } else {
      // 古いデータを消す
      while (this->getNumOfItems() >= this->getMaxItems()) {
        KeyListItr reqHistoryTail = this->reqHistory_.end();
        --reqHistoryTail;
        this->cacheTbl_.erase(*reqHistoryTail);
        this->reqHistory_.erase(reqHistoryTail);
        --(this->numOfItems_);

        // update iterator
        pCache = this->cacheTbl_.lower_bound(key);
      }

      // create new data
      CacheType cacheData =
          std::make_pair(TlSharedPointer<DataType>(new DataType()),
                         this->reqHistory_.insert(reqHistoryHead, key));
      this->cacheTbl_.insert(pCache, std::make_pair(key, cacheData));
      ++(this->numOfItems_);

      return cacheData.first;
    }
  }

  bool find(const KeyType& key) {
    bool answer = false;
    typename CacheTable::iterator pCache = this->cacheTbl_.find(key);
    if (pCache != this->cacheTbl_.end()) {
      KeyListItr reqHistoryHead = this->reqHistory_.begin();
      this->reqHistory_.splice(reqHistoryHead, this->reqHistory_,
                               pCache->second.second);
      answer = true;
    }
    return answer;
  }

  std::size_t getMaxItems() const { return this->maxItems_; }

  void setMaxItems(const std::size_t maxItems) { this->maxItems_ = maxItems; }

  std::size_t getNumOfItems() const { return this->numOfItems_; }

 private:
  CacheTable cacheTbl_;
  KeyList reqHistory_;

  std::size_t numOfItems_;
  std::size_t maxItems_;
};

#endif  // TLSTLUTILS_H
