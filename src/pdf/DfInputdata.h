#ifndef DFINPUTDATA_H
#define DFINPUTDATA_H

#include <string>
#include "TlSerializeData.h"

/// 入力データの解析を行う
class DfInputdata {
public:
    DfInputdata();
    virtual ~DfInputdata();
    
    // メインルーチン;
    virtual TlSerializeData main();

protected:
    void show(const TlSerializeData& data) const;

protected:
    mutable TlSerializeData data_;
};

#endif // DFINPUTDATA_H
