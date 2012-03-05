#ifndef TLVECTOROBJECT_H
#define TLVECTOROBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <fstream>

/** 行列クラスのインターフェースを規定する
 */
class TlVectorObject {
public:
    typedef signed int index_type; /// 行数、列数を表す整数型(範囲外を示す値として-1を取る場合がある)
    typedef long size_type; /// 数値配列の総数を表す整数型
    
public:
    /** 仮想デストラクタ
     */
    virtual ~TlVectorObject() {
    }

public:
    virtual void resize(size_type nSize) =0;

public:
    virtual void add(const size_type index, const double value) =0;

    /// 要素を返す(読み取り専用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    virtual double operator[](size_type index) const =0;

    /// 要素を返す(代入用)
    ///
    /// @param[in] index 要素番号
    /// @return 要素
    virtual double& operator[](size_type index) =0;

public:
    virtual bool load(std::ifstream& ifs) =0;
    virtual bool save(std::ofstream& ofs) const =0;

    virtual bool save(const std::string& str) const =0;
    virtual bool load(const std::string& str) =0;
};

#endif // TLVECTOROBJECT_H

