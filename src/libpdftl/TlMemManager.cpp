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

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <algorithm>
#include <cassert>
#include <iostream>

#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "TlFile.h"
#include "TlMemManager.h"
#include "TlSystem.h"
#include "TlUtils.h"

bool TlMemManager::isSetParam_ = false;
TlMemManager* TlMemManager::instance_ = NULL;
std::string TlMemManager::filePath_ = "memfile";
std::size_t TlMemManager::needMemSize_ =
    std::size_t(2 * 1024UL * 1024UL * 1024UL);  // 2 GB

void TlMemManager::setParam(const std::size_t needMemSize,
                            const std::string& baseFileName) {
  if (TlMemManager::isSetParam_ == true) {
    TlLogging& log = TlLogging::getInstance();
    log.warn("memory manager has already been setup.");
  }

  TlMemManager::needMemSize_ = needMemSize;
  TlMemManager::filePath_ =
      TlUtils::format("%s.%s%d.map", baseFileName.c_str(),
                      TlSystem::getHostName().c_str(), TlSystem::getPID());
  TlMemManager::isSetParam_ = true;
}

TlMemManager& TlMemManager::getInstance() {
  if (TlMemManager::instance_ == NULL) {
    if (TlMemManager::isSetParam_ == false) {
      TlLogging& log = TlLogging::getInstance();
      log.warn("memory manager is initialized using default parameter.");
      log.warn(
          "please call TlMemManager::setParam() function for performance.");
    }

    TlMemManager::instance_ = new TlMemManager();
  }

  return *(TlMemManager::instance_);
}

TlMemManager::TlMemManager() : log_(TlLogging::getInstance()) {
  this->mmapLength_ = 0;
  this->beginMmap_ = NULL;

  this->memMap();

  this->freeMemList_.clear();
  this->freeMemList_.push_back(MemItem(0, this->mmapLength_));
}

TlMemManager::~TlMemManager() {
  this->memUnmap();
  TlFile::remove(TlMemManager::filePath_);
}

char* TlMemManager::allocate(const std::size_t size) {
  char* pAnswer = NULL;

  // freeMemList_を空きメモリ量の少ない順にソート
  std::sort(this->freeMemList_.begin(), this->freeMemList_.end(),
            MemItemSortFunctor_size());

  // 空きメモリリストから必要量のメモリを確保する
  std::vector<MemItem>::iterator p =
      std::lower_bound(this->freeMemList_.begin(), this->freeMemList_.end(),
                       MemItem(0, size), MemItemSortFunctor_size());
  if (p != this->freeMemList_.end()) {
    // メモリ確保
    assert(size <= (p->end - p->begin));
    this->useList_[p->begin] = size;

    pAnswer = this->beginMmap_ + p->begin;
    p->begin += size;

    if (p->begin == p->end) {
      this->freeMemList_.erase(p);
    }
  } else {
    this->log_.critical(TlUtils::format(
        "TlMemManager::allocate(): cannot allocate mem. size=%lu", size));
    throw std::bad_alloc();
    std::abort();
  }

  this->checkFreeList();

  // this->log_.debug(TlUtils::format("TlMemManager::allocated: size=%16lu
  // [addr=%16lu, index=%16lu]",
  //                                  size,
  //                                  reinterpret_cast<uintptr_t>(pAnswer),
  //                                  (pAnswer - this->beginMmap_)));
  // this->debugOutUseMemList();

  return pAnswer;
}

void TlMemManager::deallocate(char* p) {
  // this->log_.debug(TlUtils::format("TlMemManager::deallocate() p=%16ld,
  // size=%16ld",
  //                                  reinterpret_cast<uintptr_t>(p), size));

  const std::size_t begin = p - this->beginMmap_;
  const std::size_t size = this->useList_[begin];
  const std::size_t end = begin + size;
  const MemItem mi(begin, end);

  this->useList_.erase(begin);

  // 解放領域を追加
  this->freeMemList_.push_back(mi);

  // freeMemList_をindex順にソート
  std::sort(this->freeMemList_.begin(), this->freeMemList_.end(),
            MemItemSortFunctor_index());

  std::vector<MemItem>::iterator pBegin = this->freeMemList_.begin();
  std::vector<MemItem>::iterator pEnd = this->freeMemList_.end();
  assert(pBegin->begin <= pBegin->end);
  std::size_t prevEnd = pBegin->end;
  std::vector<MemItem>::iterator curr = pBegin + 1;
  while (curr != pEnd) {
    assert(curr->begin < curr->end);

    if (prevEnd == curr->begin) {
      // 前のブロックに連結する
      std::size_t newEnd = curr->end;
      std::vector<MemItem>::iterator prev = curr - 1;
      prev->end = newEnd;
      curr = this->freeMemList_.erase(curr);
      prevEnd = newEnd;
      continue;
    }

    prevEnd = curr->end;
    ++curr;
  }

  this->checkFreeList();
  p = NULL;

  // this->log_.debug("TlMemManager::deallocated:");
  // this->debugOutUseMemList();
}

void TlMemManager::memMap() {
  // file descriptor
  const int fd = open(this->filePath_.c_str(), O_RDWR | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR);
  if (fd == -1) {
    int errorNo = errno;
    std::string errStr(strerror(errorNo));
    this->log_.critical(TlUtils::format(
        "open error: %s (%s)", this->filePath_.c_str(), errStr.c_str()));
    abort();
  }

  // page aligned.
  const std::size_t pageSize = sysconf(_SC_PAGE_SIZE);
  const size_t bufSize = TlMemManager::needMemSize_;
  TlMemManager::mmapLength_ = ((bufSize / pageSize) + 1) * pageSize;

  // fill up to file size.
  {
    lseek(fd, this->mmapLength_, SEEK_SET);
    const char c = '\0';
    write(fd, &c, sizeof(char));
  }

  // message
  this->log_.info(
      TlUtils::format("TlMemManager::memMap() file create: path=%s, size=%lu",
                      this->filePath_.c_str(), TlMemManager::mmapLength_));

  TlMemManager::beginMmap_ = (char*)mmap(
      0, this->mmapLength_, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if ((TlMemManager::beginMmap_ != NULL) && (this->beginMmap_ == MAP_FAILED)) {
    // error
    perror("mmap");
  }

  // we can close the file after mmap() was called.
  close(fd);
}

void TlMemManager::sync() {
  msync(this->beginMmap_, this->mmapLength_, MS_SYNC);
}

void TlMemManager::memUnmap() {
  this->sync();
  munmap(this->beginMmap_, this->mmapLength_);
  this->beginMmap_ = NULL;
}

void TlMemManager::checkFreeList() {
// debug用チェックルーチン
#ifndef NDEBUG
  std::sort(this->freeMemList_.begin(), this->freeMemList_.end(),
            MemItemSortFunctor_index());

  bool isWrong = false;
  std::size_t currentPos = 0;
  std::vector<MemItem>::const_iterator pEnd = this->freeMemList_.end();
  for (std::vector<MemItem>::const_iterator p = this->freeMemList_.begin();
       p != pEnd; ++p) {
    if ((p->begin < currentPos) || (p->end <= p->begin)) {
      isWrong = true;
      break;
    }

    currentPos = p->end;
  }

  if (isWrong == true) {
    this->log_.critical(
        "TlMemManager::checkFreeList(): check failed. someting wrong.");
    this->debugOutFreeMemList();
    abort();
  }
#endif  // NDEBUG
}

void TlMemManager::debugOutFreeMemList() const {
  std::sort(this->freeMemList_.begin(), this->freeMemList_.end(),
            MemItemSortFunctor_index());

  std::string str = "";
  std::vector<MemItem>::const_iterator pEnd = this->freeMemList_.end();
  for (std::vector<MemItem>::const_iterator p = this->freeMemList_.begin();
       p != pEnd; ++p) {
    str += TlUtils::format("s=%16lu, e=%16lu, f=%16lu\n", p->begin, p->end,
                           (p->end - p->begin));
  }
  str += "\n";

  this->log_.debug(str);
}

void TlMemManager::debugOutUseMemList() const {
  std::string str = "";
  for (std::map<std::size_t, std::size_t>::const_iterator it =
           this->useList_.begin();
       it != this->useList_.end(); ++it) {
    str += TlUtils::format("USE: %16lu -> %16lu\n", it->first, it->second);
  }

  this->log_.debug(str);
}
