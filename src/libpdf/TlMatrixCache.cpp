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

#include "TlMatrixCache.h"
#include "TlUtils.h"

// static member variables
std::map<std::string, TlMatrixCache::Stats> TlMatrixCache::stats_;

TlMatrixCache::TlMatrixCache(std::size_t maxMemSize,
                             bool isForceLoadingFromDisk)
    : maxMemSize_(maxMemSize),
      usedMemSize_(0),
      isForceLoadingFromDisk_(isForceLoadingFromDisk) {
  // this->isDebugOut_ = true;
}

TlMatrixCache::~TlMatrixCache() {
  LRU_List::iterator itEnd = this->lruList_.end();
  for (LRU_List::iterator it = this->lruList_.begin(); it != itEnd;) {
    it = this->erase(it);
  }
}

void TlMatrixCache::setMaxMemSize(std::size_t memSize) {
  if (memSize < this->maxMemSize_) {
    // shrink memory
    this->cleanOut(memSize - this->maxMemSize_);
  }

  this->maxMemSize_ = memSize;
}

std::size_t TlMatrixCache::getFreeMemSize() const {
  return this->maxMemSize_ - this->usedMemSize_;
}

void TlMatrixCache::forceLoadingFromDisk(bool flag) {
  this->isForceLoadingFromDisk_ = flag;
}

void TlMatrixCache::erase(const std::string& path) {
  LRU_List::iterator itEnd = this->lruList_.end();
  for (LRU_List::iterator it = this->lruList_.begin(); it != itEnd; ++it) {
    const std::string it_path = (*it)->getPath();
    if (it_path == path) {
      it = this->erase(it);
      break;
    }
  }
}

TlMatrixCache::LRU_List::iterator TlMatrixCache::erase(
    const LRU_List::iterator& it) {
  this->flush(it);

  this->usedMemSize_ -= (*it)->getSize();

  delete *it;
  *it = NULL;

  return this->lruList_.erase(it);
}

void TlMatrixCache::cleanOut(const std::size_t needMemSize) {
  if (needMemSize < this->maxMemSize_) {
    // 必要量が最大容量を下回っているときにしか、
    // メモリを空ける意味がない。
    while ((this->getFreeMemSize() < needMemSize) ||
           (this->lruList_.size() > 1)) {
      LRU_List::iterator it = --(this->lruList_.end());
      this->erase(it);
    }
  }
}

void TlMatrixCache::flush() {
  LRU_List::iterator itEnd = this->lruList_.end();
  for (LRU_List::iterator it = this->lruList_.begin(); it != itEnd; ++it) {
    this->flush(it);
  }
}

void TlMatrixCache::flush(const std::string& path) {
  LRU_List::iterator itEnd = this->lruList_.end();
  for (LRU_List::iterator it = this->lruList_.begin(); it != itEnd; ++it) {
    if ((*it)->getPath() == path) {
      this->flush(it);
    }
  }
}

void TlMatrixCache::flush(LRU_List::iterator it) {
  if ((*it)->needSave()) {
    const std::string path = (*it)->getPath();
    (*it)->save();
    ++(this->stats_[path].setToDisk);
  }
}

std::string TlMatrixCache::reportStats() {
  std::string answer = "";
  std::map<std::string, Stats>::const_iterator pEnd =
      TlMatrixCache::stats_.end();
  for (std::map<std::string, Stats>::const_iterator p =
           TlMatrixCache::stats_.begin();
       p != pEnd; ++p) {
    const std::string path = p->first;
    const Stats stats = p->second;
    std::string headStr = TlUtils::format("  %s: ", path.c_str());
    std::string dataStr = TlUtils::format(
        "set(%2d/%2d/%2d) get(%2d/%2d/%2d)", stats.setTotal, stats.setToCache,
        stats.setToDisk, stats.getTotal, stats.getFromCache, stats.getFromDisk);
    std::string ws = "";
    TlUtils::pad(ws, 80 - (headStr.size() + dataStr.size()), ' ');
    answer += headStr + ws + dataStr + "\n";
  }

  return answer;
}
