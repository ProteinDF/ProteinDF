#ifndef TLMSGPACK_H
#define TLMSGPACK_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <string>
#include <sstream>
#include "pdflib-int.h"
#include "TlSerializeData.h"

class TlMsgPack {
public:
    typedef int8_t INT8;
    typedef int16_t INT16;
    typedef int32_t INT32;
    typedef int64_t INT64;
    typedef uint8_t UINT8;
    typedef uint16_t UINT16;
    typedef uint32_t UINT32;
    typedef uint64_t UINT64;

public:
    explicit TlMsgPack(const TlSerializeData& data = TlSerializeData());
    TlMsgPack(const TlMsgPack& rhs);
    ~TlMsgPack();

    TlMsgPack& operator=(const TlMsgPack& rhs);

public:
    TlSerializeData getSerializeData() const;

    /// MsgPack形式のファイルを読み込む
    ///
    /// @retval true  ファイルの読み込みに成功した
    /// @retval false ファイルの読み込みに失敗した
    bool load(const std::string& path);

    void save(const std::string& path) const;

    void pack(const std::string& str);
    std::string dump() const;

protected:
    TlSerializeData loadBinary(std::istream& ifs);
    int unpack_positiveFixNum(unsigned char c);
    int unpack_negativeFixNum(unsigned char c);
    UINT8 unpack_uint8(std::istream& ifs);
    UINT16 unpack_uint16(std::istream& ifs);
    UINT32 unpack_uint32(std::istream& ifs);
    UINT64 unpack_uint64(std::istream& ifs);
    
    INT8 unpack_int8(std::istream& ifs);
    INT16 unpack_int16(std::istream& ifs);
    INT32 unpack_int32(std::istream& ifs);
    INT64 unpack_int64(std::istream& ifs);

    float unpack_float(std::istream& ifs);
    double unpack_double(std::istream& ifs);
    
    TlSerializeData unpack_fixraw(const char in, std::istream& ifs);
    TlSerializeData unpack_raw16(std::istream& ifs);
    TlSerializeData unpack_raw32(std::istream& ifs);

    TlSerializeData unpack_fixarray(const char in, std::istream& ifs);
    TlSerializeData unpack_array16(std::istream& ifs);
    TlSerializeData unpack_array32(std::istream& ifs);
    TlSerializeData unpack_fixmap(const char in, std::istream& ifs);
    TlSerializeData unpack_map16(std::istream& ifs);
    TlSerializeData unpack_map32(std::istream& ifs);

    std::string dump(const TlSerializeData& data) const;
    std::string dump_scalar(const TlSerializeData& data) const;
    std::string dump_array(const TlSerializeData& data) const;
    std::string dump_map(const TlSerializeData& data) const;

    std::string pack(bool value) const;
    std::string pack(UINT8 value) const;
    std::string pack(UINT16 value) const;
    std::string pack_uint32(UINT32 value) const;
    std::string pack_uint64(UINT64 value) const;
    std::string pack(INT8 value) const;
    std::string pack(INT16 value) const;
    std::string pack_int32(INT32 value) const;
    std::string pack_int64(INT64 value) const;
    std::string pack(double value) const;
    std::string pack(const char* pBuf, const int size) const;

    template<typename T>
    void write(std::ostringstream& os, T value) const {
        os.write((char*)&value, sizeof(T));
    }
    
protected:
    TlSerializeData data_;

    /// デバッグ用変数
    /// 現在の読み込み位置(byte)を記憶する
    std::size_t debugCurrentPos_;
};

#endif // TLMSGPACK_H
