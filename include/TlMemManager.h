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

#ifndef TLMEMMANAGER_H
#define TLMEMMANAGER_H

#include <vector>
#include <map>
#include "TlLogging.h"

class TlMemManager {
public:
    /// メモリブロック毎の情報を保持する構造体
    struct MemItem {
public:
        MemItem(std::size_t b =0, std::size_t e =0) : begin(b), end(e) {
        }

public:
        /// 開始位置
        std::size_t begin;

        /// 終了位置(このアドレスは含まない)
        std::size_t end;
    };

    struct MemItemSortFunctor_index {
        bool operator()(const MemItem& a, const MemItem& b) {
            return (a.begin < b.begin);
        }
    };

    struct MemItemSortFunctor_size {
        bool operator()(const MemItem& a, const MemItem& b) {
            return ((a.end - a.begin) < (b.end - b.begin));
        }
    };

public:
    static void setParam(std::size_t needMemSize, const std::string& baseFileName);
    static TlMemManager& getInstance();

    char* allocate(std::size_t size);
    void deallocate(char* p);

    void debugOutFreeMemList() const;
    void debugOutUseMemList() const;

private:
    void memMap();
    void sync();
    void memUnmap();

    void checkFreeList();

private:
    // for singleton
    TlMemManager();
    TlMemManager(const TlMemManager& rhs);
    ~TlMemManager();
    const TlMemManager& operator=(const TlMemManager& rhs);

private:
    TlLogging& log_;

    /// パラメータの設定が行われたかどうかを保持するフラグ
    /// setParam() が呼び出されるとtrue
    static bool isSetParam_;

    /// 唯一のオブジェクトを保持するポインタ
    static TlMemManager* instance_;

    /// 仮想メモリとして利用するファイル
    static std::string filePath_;

    /// 必要なメモリサイズ(byte単位)
    static std::size_t needMemSize_;

    /// 実際に確保したmmap長
    std::size_t mmapLength_;

    /// memmap先頭ポインタ
    char* beginMmap_;

    mutable std::vector<MemItem> freeMemList_;
    mutable std::map<std::size_t, std::size_t> useList_; // deallocate時に返すメモリ量が必要なため
};

#endif // TLMEMMANAGER_H

