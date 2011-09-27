#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <cassert>
#include <iostream>
#include <algorithm>

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

#include "TlMemManager.h"
#include "TlFile.h"
#include "TlUtils.h"

bool TlMemManager::isSetParam_ = false;
TlMemManager* TlMemManager::instance_ = NULL;
std::string TlMemManager::filePath_ = "memfile";
std::size_t TlMemManager::needMemSize_ = std::size_t(2 * 1024UL * 1024UL * 1024UL); // 2 GB

void TlMemManager::setParam(const std::size_t needMemSize,
                            const std::string& baseFileName)
{
    if (TlMemManager::isSetParam_ == true) {
        std::cerr << "memory manager has already been setup." << std::endl;
    }

    TlMemManager::needMemSize_ = needMemSize;
    const pid_t pid = getpid();
    TlMemManager::filePath_ = baseFileName + TlUtils::xtos(pid);
    TlMemManager::isSetParam_ = true;
}


TlMemManager& TlMemManager::getInstance()
{
    if (TlMemManager::instance_ == NULL) {
        if (TlMemManager::isSetParam_ == false) {
            std::cerr << "memory manager is initialized using default parameter." << std::endl;
            std::cerr << "please call TlMemManager::setParam() function for performance." << std::endl;
        }

        TlMemManager::instance_ = new TlMemManager();
    }

    return *(TlMemManager::instance_);
}


TlMemManager::TlMemManager()
{
    this->mmapLength_ = 0;
    this->beginMmap_ = NULL;

    this->memMap();

    this->freeMemList_.clear();
    this->freeMemList_.push_back(MemItem(0, this->mmapLength_));
}


TlMemManager::~TlMemManager()
{
    this->memUnmap();
    TlFile::remove(TlMemManager::filePath_);
}


char* TlMemManager::allocate(const std::size_t size)
{
    char* pAnswer = NULL;

    // freeMemList_を空きメモリ量の少ない順にソート
    std::sort(this->freeMemList_.begin(),
              this->freeMemList_.end(),
              MemItemSortFunctor_size());

    // 空きメモリリストから必要量のメモリを確保する
    std::vector<MemItem>::iterator p = std::lower_bound(this->freeMemList_.begin(),
                                                        this->freeMemList_.end(),
                                                        MemItem(0, size),
                                                        MemItemSortFunctor_size());
    if (p != this->freeMemList_.end()) {
        // メモリ確保
        assert(size <= (p->end - p->begin));
        this->useList_[p->begin] = size;

        pAnswer = this->beginMmap_ + p->begin;
        p->begin += size;

        if (p->begin == p->end) {
            this->freeMemList_.erase(p);
        }
    }

    this->checkFreeList();
    return pAnswer;
}


void TlMemManager::deallocate(char* p, const std::size_t size)
{
    const std::size_t begin = p - this->beginMmap_;
    const std::size_t end = begin + size;
    const MemItem mi(begin, end);

    assert(this->useList_[begin] == size);
    this->useList_.erase(begin);

    // 解放領域を追加
    this->freeMemList_.push_back(mi);

    // freeMemList_をindex順にソート
    std::sort(this->freeMemList_.begin(),
              this->freeMemList_.end(),
              MemItemSortFunctor_index());

    std::vector<MemItem>::iterator pBegin = this->freeMemList_.begin();
    std::vector<MemItem>::iterator pEnd = this->freeMemList_.end();
    assert(pBegin->begin <= pBegin->end);
    std::size_t prevEnd = pBegin->end;
    std::vector<MemItem>::iterator curr = pBegin +1;
    while (curr != pEnd) {
        assert(curr->begin < curr->end);

        if (prevEnd == curr->begin) {
            // 前のブロックに連結する
            std::size_t newEnd = curr->end;
            std::vector<MemItem>::iterator prev = curr -1;
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
}


void TlMemManager::memMap()
{
    // file descriptor
    const int fd = open(this->filePath_.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        int errorNo = errno;
        std::string errStr(strerror(errorNo));
        std::cerr << " open error: " << this->filePath_ << ". " << errStr << std::endl;
        abort();
    }

    // page aligned.
    const std::size_t pageSize = sysconf(_SC_PAGE_SIZE);
    const size_t bufSize = TlMemManager::needMemSize_;
    TlMemManager::mmapLength_ = ((bufSize / pageSize) +1) * pageSize;

    // fill up to file size.
    {
        lseek(fd, this->mmapLength_, SEEK_SET);
        const char c = '\0';
        write(fd, &c, sizeof(char));
    }

    // message
//   std::cerr << TlUtils::format("TlMemManager::memMap() file create.: path=%s, size=%ld",
//                 filePath.c_str(),
//                 TlMemManager::mmapLength_)
//      << std::endl;


    TlMemManager::beginMmap_ = (char*)mmap(0, this->mmapLength_,
                                           PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if ((TlMemManager::beginMmap_ != NULL) &&
            (this->beginMmap_ == MAP_FAILED)) {
        // error
        perror("mmap");
    }

    // we can close the file after mmap() was called.
    close(fd);
}


void TlMemManager::sync()
{
    msync(this->beginMmap_, this->mmapLength_, MS_SYNC);
}


void TlMemManager::memUnmap()
{
    this->sync();
    munmap(this->beginMmap_, this->mmapLength_);
    this->beginMmap_ = NULL;
}


void TlMemManager::debugOutFreeMemList() const
{
    std::sort(this->freeMemList_.begin(),
              this->freeMemList_.end(),
              MemItemSortFunctor_index());

    std::string str = "";
    std::vector<MemItem>::const_iterator pEnd = this->freeMemList_.end();
    for (std::vector<MemItem>::const_iterator p = this->freeMemList_.begin();
            p != pEnd; ++p) {

        str += TlUtils::format("s=%16ld, e=%16ld, f=%16ld\n",
                               p->begin, p->end, (p->end - p->begin));
    }
    str += "\n";
    std::cerr << str << std::endl;
}


void TlMemManager::checkFreeList()
{
    // debug用チェックルーチン
#ifndef NDEBUG
    std::sort(this->freeMemList_.begin(),
              this->freeMemList_.end(),
              MemItemSortFunctor_index());

    bool isWrong = false;
    std::size_t currentPos = 0;
    std::vector<MemItem>::const_iterator pEnd = this->freeMemList_.end();
    for (std::vector<MemItem>::const_iterator p = this->freeMemList_.begin(); p != pEnd; ++p) {
        if ((p->begin < currentPos) ||
                (p->end <= p->begin)) {
            isWrong = true;
            break;
        }

        currentPos = p->end;
    }

    if (isWrong == true) {
        std::cerr << "TlMemManager::checkFreeList(): check failed. someting wrong." << std::endl;
        //TlMemManager::debugOutFreeMemList();
        abort();
    }
#endif // NDEBUG
}

