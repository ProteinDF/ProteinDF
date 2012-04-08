#ifndef TLMATRIXCACHE_H
#define TLMATRIXCACHE_H

#include <map>
#include <list>
#include <string>
#include "TlMatrixObject.h"
#include "TlUtils.h"
#include "TlFile.h"
#include "TlLogging.h"

class TlMatrixCache {
public:
    struct Stats {
    public:
        Stats() : setTotal(0), getTotal(0),
                  setToCache(0), getFromCache(0), setToDisk(0), getFromDisk(0) {
        }
    public:
        int setTotal;
        int getTotal;
        int setToCache;
        int getFromCache;
        int setToDisk;
        int getFromDisk;
    };
    
public:
    TlMatrixCache(std::size_t maxMemSize = 0, bool isForceLoadingMatrixFromDisk = false);
    ~TlMatrixCache();

public:
    void setMaxMemSize(std::size_t memSize);
    std::size_t getFreeMemSize() const;

    /// 強制的に行列をディスクから読み取る場合はtrueに設定する
    void forceLoadingFromDisk(bool flag);

    template<typename MatrixType>
    bool set(const std::string& path,
             const MatrixType& matrix,
             bool needSave = false);

    template<typename MatrixType>
    MatrixType get(const std::string& path);

    /// 指定の行列オブジェクトをメモリから削除する
    ///
    /// 書き出しフラグが付いたものはファイルに書き出してからメモリから削除する
    void erase(const std::string& path);

    /// 書き出しフラグが付いた行列オブジェクトを強制的にファイルに書き出す
    void flush();

    static std::string reportStats();

private:
    void cleanOut(const std::size_t needMemSize);
    void flush(const std::string& path);
       
private:
    struct MatrixInfo {
    public:
        explicit MatrixInfo(TlMatrixObject* buf = NULL, std::size_t size = 0, bool ns = false)
            : pBuf(buf), bufSize(size), needSave(ns) {
        }
        MatrixInfo(const MatrixInfo& rhs)
            : pBuf(rhs.pBuf), bufSize(rhs.bufSize), needSave(rhs.needSave) {
        }
        
    public:
        TlMatrixObject* pBuf;
        std::size_t bufSize;
        bool needSave;
    };
    
private:
    bool isDebugOut_;
    std::size_t maxMemSize_;
    std::size_t usedMemSize_;

    /// 強制的にディスクから読み取る場合はtrue
    bool isForceLoadingFromDisk_;
    
    typedef std::map<std::string, MatrixInfo> CACHE_DATA_TYPE;
    CACHE_DATA_TYPE cache_;
    
    // LRU方式のリスト
    std::list<std::string> lruList_;

    static std::map<std::string, Stats> stats_;
};


template<typename MatrixType>
bool TlMatrixCache::set(const std::string& path,
                        const MatrixType& matrix,
                        bool needSave)
{
    bool answer = false;
    const std::size_t bufSize = matrix.getMemSize();
    this->cleanOut(bufSize);

    Stats stats = this->stats_[path];
    ++(stats.setTotal);
    
    if (bufSize < this->getFreeMemSize()) {
        if (this->isDebugOut_ == true) {
            std::cerr << TlUtils::format("TlMatrixCache::set() to memory: %s, size=%ld",
                                         path.c_str(), bufSize)
                      << std::endl;
        }
        TlMatrixObject* pNew = new MatrixType(matrix);
        MatrixInfo info(pNew, bufSize, needSave);

        // 既にキーが存在している場合はメモリ確保容量を返却する
        if (this->cache_.find(path) != this->cache_.end()) {
            this->usedMemSize_ -= this->cache_[path].bufSize;
        }
        this->cache_[path] = info;
        this->usedMemSize_ += bufSize;

        // LRUリストの最後尾に登録する
        std::list<std::string>::iterator lruListIt = std::find(this->lruList_.begin(),
                                                               this->lruList_.end(),
                                                               path);
        if (lruListIt != this->lruList_.end()) {
            this->lruList_.erase(lruListIt);
        }
        this->lruList_.push_back(path);
        assert(this->cache_.size() == this->lruList_.size());

        ++(stats.setToCache);
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
        ++(stats.setToDisk);
        matrix.save(path);
        needSave = false;
    }

    this->stats_[path] = stats;
    
    return answer;
}


template<typename MatrixType>
MatrixType TlMatrixCache::get(const std::string& path)
{
    Stats& stats = this->stats_[path];
    ++(stats.getTotal);
    
    MatrixType matrix;

    bool isPrepareMatrix = false;
    if (this->isForceLoadingFromDisk_ == false) {
        CACHE_DATA_TYPE::iterator p = this->cache_.find(path);
        if (p != this->cache_.end()) {
            matrix = *(dynamic_cast<MatrixType*>(p->second.pBuf));
            ++(stats.getFromCache);
            
            // LRUリストの最後尾に登録する
            std::list<std::string>::iterator lruListIt = std::find(this->lruList_.begin(),
                                                                   this->lruList_.end(),
                                                                   path);
            if (lruListIt != this->lruList_.end()) {
                this->lruList_.erase(lruListIt);
            }
            this->lruList_.push_back(path);
            isPrepareMatrix = true;
        }
    }

    if (isPrepareMatrix == false) {
        // cacheに無いのでloadしてくる
        if (TlFile::isExist(path) == true) {
            if (matrix.load(path) == true) {
                if (this->isDebugOut_ == true) {
                    std::cerr << TlUtils::format("TlMatrixCache::get() from disk: %s",
                                                 path.c_str())
                              << std::endl;
                }
                this->set(path, matrix, false);
                ++(stats.getFromDisk);
            }
        } else {
            TlLogging& log = TlLogging::getInstance();
            log.critical(TlUtils::format("TlMatrixCache::get(): file not found: %s",
                                         path.c_str()));
#ifndef NDEBUG
            abort();
#endif // NDEBUG
        }
    }

    //this->stats_[path] = stats;
    return matrix;
}


#endif // TLMATRIXCACHE_H
