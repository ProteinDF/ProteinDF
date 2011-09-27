#include <iostream>
#include "TlMatrixCache.h"
#include "TlUtils.h"

// static member variables
std::map<std::string, TlMatrixCache::Stats> TlMatrixCache::stats_;


TlMatrixCache::TlMatrixCache(std::size_t maxMemSize, bool isForceLoadingFromDisk)
    : maxMemSize_(maxMemSize), usedMemSize_(0),
      isForceLoadingFromDisk_(isForceLoadingFromDisk) {
    //this->isDebugOut_ = true;
}


TlMatrixCache::~TlMatrixCache()
{
    this->flush();
}


void TlMatrixCache::setMaxMemSize(std::size_t memSize)
{
    if (memSize < this->maxMemSize_) {
        // shrink memory
        this->cleanOut(memSize - this->maxMemSize_);
    }
    
    this->maxMemSize_ = memSize;
}


std::size_t TlMatrixCache::getFreeMemSize() const
{
    return this->maxMemSize_ - this->usedMemSize_;
}


void TlMatrixCache::forceLoadingFromDisk(bool flag)
{
    this->isForceLoadingFromDisk_ = flag;
}

void TlMatrixCache::erase(const std::string& path)
{
    assert(this->cache_.size() == this->lruList_.size());
    CACHE_DATA_TYPE::iterator p = this->cache_.find(path);
    if (p != this->cache_.end()) {
        const std::size_t bufSize = p->second.bufSize;
        this->flush(path);

        if (this->isDebugOut_ == true) {
            std::cerr << TlUtils::format("TlMatrixCache::erase() from mem: %s",
                                         path.c_str())
                      << std::endl;
        }
        delete p->second.pBuf;
        p->second.pBuf = NULL;
        this->cache_.erase(p);
        this->usedMemSize_ -= bufSize;

        std::list<std::string>::iterator lruListIt = std::find(this->lruList_.begin(),
                                                               this->lruList_.end(),
                                                               path);
        if (lruListIt != this->lruList_.end()) {
            this->lruList_.erase(lruListIt);
        }
    }
    assert(this->cache_.size() == this->lruList_.size());
}


void TlMatrixCache::cleanOut(const std::size_t needMemSize)
{
    assert(this->cache_.size() == this->lruList_.size());
    if (needMemSize < this->maxMemSize_) {
        // 必要量が最大容量を下回っているときにしか、
        // メモリを空ける意味がない。
        while ((this->getFreeMemSize() < needMemSize) ||
               (this->lruList_.size() > 1)) {
            const std::string path = this->lruList_.front();
            this->erase(path);
        }
    }
    assert(this->cache_.size() == this->lruList_.size());
}


void TlMatrixCache::flush()
{
    CACHE_DATA_TYPE::iterator itEnd = this->cache_.end();
    for (CACHE_DATA_TYPE::iterator it = this->cache_.begin(); it != itEnd; ++it) {
        this->flush(it->first);
    }
}


void TlMatrixCache::flush(const std::string& path)
{
    CACHE_DATA_TYPE::iterator it = this->cache_.find(path);
    if (it != this->cache_.end()) {
        if (it->second.needSave == true) {
            it->second.pBuf->save(path);
            it->second.needSave = false;

            // stats
            ++(this->stats_[it->first].setToDisk);
        }
    }
}


std::string TlMatrixCache::reportStats()
{
    std::string answer = "";
    std::map<std::string, Stats>::const_iterator pEnd = TlMatrixCache::stats_.end();
    for (std::map<std::string, Stats>::const_iterator p = TlMatrixCache::stats_.begin(); p != pEnd; ++p) {
        const std::string path = p->first;
        const Stats stats = p->second;
        std::string headStr = TlUtils::format("  %s: ",
                                              path.c_str());
        std::string dataStr = TlUtils::format("set(%2d/%2d/%2d) get(%2d/%2d/%2d)",
                                              stats.setTotal, stats.setToCache, stats.setToDisk,
                                              stats.getTotal, stats.getFromCache, stats.getFromDisk);
        std::string ws = "";
        TlUtils::pad(ws, 80 - (headStr.size() + dataStr.size()), ' ');
        answer += headStr + ws + dataStr + "\n";
    }

    return answer;
}

