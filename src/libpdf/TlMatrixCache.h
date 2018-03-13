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

#ifndef TLMATRIXCACHE_H
#define TLMATRIXCACHE_H

#include <iostream>
#include <list>
#include <map>
#include <string>
#include "TlFile.h"
#include "TlLogging.h"
#include "TlUtils.h"
#include "tl_matrix_object.h"

class TlMatrixCache {
 public:
  struct Stats {
   public:
    Stats()
        : setTotal(0),
          getTotal(0),
          setToCache(0),
          getFromCache(0),
          setToDisk(0),
          getFromDisk(0) {}

   public:
    int setTotal;
    int getTotal;
    int setToCache;
    int getFromCache;
    int setToDisk;
    int getFromDisk;
  };

 public:
  TlMatrixCache(std::size_t maxMemSize = 0,
                bool isForceLoadingMatrixFromDisk = false);
  ~TlMatrixCache();

 public:
  void setMaxMemSize(std::size_t memSize);
  std::size_t getFreeMemSize() const;

  /// 強制的に行列をディスクから読み取る場合はtrueに設定する
  void forceLoadingFromDisk(bool flag);

  template <typename MatrixType>
  bool set(const std::string& path, const MatrixType& matrix,
           bool needSave = false);

  template <typename MatrixType>
  MatrixType get(const std::string& path);

  /// 指定の行列オブジェクトをメモリから削除する
  ///
  /// 書き出しフラグが付いたものはファイルに書き出してからメモリから削除する
  void erase(const std::string& path);

  /// 書き出しフラグが付いた行列オブジェクトを強制的にファイルに書き出す
  void flush();

  static std::string reportStats();

 private:
  struct MatrixInfo {
   public:
    template <typename MatrixType>
    MatrixInfo(const std::string& path, const MatrixType& mat,
               bool needSave = false)
        : path_(path), pMatrix_(NULL), needSave_(needSave) {
      this->pMatrix_ = new MatrixType(mat);
    }

    ~MatrixInfo() {
      delete this->pMatrix_;
      this->pMatrix_ = NULL;
    }

   public:
    std::string getPath() const { return this->path_; }

    std::size_t getSize() const {
      return this->pMatrix_->getNumOfRows() * this->pMatrix_->getNumOfCols() *
             sizeof(double);
    }

    TlMatrixObject const* getMatPtr() const { return this->pMatrix_; }

    bool needSave() const { return this->needSave_; }

    bool save() {
      this->needSave_ = false;
      return this->pMatrix_->save(this->path_);
    }

    bool operator==(const MatrixInfo& rhs) const {
      return (this->path_ == rhs.path_);
    }

    bool operator!=(const MatrixInfo& rhs) const {
      return !(this->operator==(rhs));
    }

   private:
    // prohibit copy constructer
    MatrixInfo(const MatrixInfo& rhs) {}

    // prohibit operator=
    void operator=(const MatrixInfo& rhs) {}

   private:
    std::string path_;
    TlMatrixObject* pMatrix_;
    bool needSave_;
  };

  typedef std::list<MatrixInfo*> LRU_List;

 private:
  void cleanOut(const std::size_t needMemSize);

  LRU_List::iterator erase(const LRU_List::iterator& it);

  void flush(const std::string& path);
  void flush(LRU_List::iterator it);

 private:
  bool isDebugOut_;
  std::size_t maxMemSize_;
  std::size_t usedMemSize_;

  /// 強制的にディスクから読み取る場合はtrue
  bool isForceLoadingFromDisk_;

  // typedef std::map<std::string, MatrixInfo*> CACHE_DATA_TYPE;
  // CACHE_DATA_TYPE cache_;

  // LRU方式のリスト
  LRU_List lruList_;

  static std::map<std::string, Stats> stats_;
};

template <typename MatrixType>
bool TlMatrixCache::set(const std::string& path, const MatrixType& matrix,
                        bool needSave) {
  bool answer = false;
  const std::size_t bufSize = matrix.getMemSize();
  this->cleanOut(bufSize);

  ++(this->stats_[path].setTotal);

  if (bufSize < this->getFreeMemSize()) {
    if (this->isDebugOut_ == true) {
      std::cerr << TlUtils::format(
                       "TlMatrixCache::set() to memory: %s, size=%ld",
                       path.c_str(), bufSize)
                << std::endl;
    }

    this->erase(path);

    MatrixInfo* pMatInfo = new MatrixInfo(path, matrix, needSave);
    this->lruList_.push_front(pMatInfo);
    this->usedMemSize_ += pMatInfo->getSize();

    ++(this->stats_[path].setToCache);
    answer = true;

    if (needSave == true) {
      this->flush(path);
    }
  } else if (needSave == true) {
    if (this->isDebugOut_ == true) {
      std::cerr << TlUtils::format("TlMatrixCache::set() save: %s",
                                   path.c_str())
                << std::endl;
    }

    matrix.save(path);
    needSave = false;
    ++(this->stats_[path].setToDisk);
  }

  return answer;
}

template <typename MatrixType>
MatrixType TlMatrixCache::get(const std::string& path) {
  ++(this->stats_[path].getTotal);

  MatrixType matrix;

  bool isPrepareMatrix = false;
  if (this->isForceLoadingFromDisk_ == false) {
    LRU_List::iterator itEnd = this->lruList_.end();
    for (LRU_List::iterator it = this->lruList_.begin(); it != itEnd; ++it) {
      if (path == (*it)->getPath()) {
        matrix = *(dynamic_cast<MatrixType const*>((*it)->getMatPtr()));
        isPrepareMatrix = true;
        ++(this->stats_[path].getFromCache);

        // listの先頭に移動
        this->lruList_.insert(this->lruList_.begin(), *it);
        this->lruList_.erase(it);

        break;
      }
    }
  }

  if (isPrepareMatrix == false) {
    // cacheに無いのでloadしてくる
    if (TlFile::isExistFile(path) == true) {
      if (matrix.load(path) == true) {
        if (this->isDebugOut_ == true) {
          std::cerr << TlUtils::format("TlMatrixCache::get() from disk: %s",
                                       path.c_str())
                    << std::endl;
        }
        this->set(path, matrix, false);
        ++(this->stats_[path].getFromDisk);
      }
    } else {
      TlLogging& log = TlLogging::getInstance();
      log.critical(TlUtils::format("TlMatrixCache::get(): file not found: %s",
                                   path.c_str()));
#ifndef NDEBUG
      abort();
#endif  // NDEBUG
    }
  }

  return matrix;
}

#endif  // TLMATRIXCACHE_H
