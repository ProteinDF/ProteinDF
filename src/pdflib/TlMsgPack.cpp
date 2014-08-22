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
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <fstream>

#include "TlMsgPack.h"
#include "TlUtils.h"

TlMsgPack::TlMsgPack(const TlSerializeData& data)
        : data_(data)
{
}

TlMsgPack::TlMsgPack(const TlMsgPack& rhs)
        : data_(rhs.data_)
{
}

TlMsgPack::~TlMsgPack()
{
}

TlMsgPack& TlMsgPack::operator=(const TlMsgPack& rhs)
{
    if (this != &rhs) {
        this->data_ = rhs.data_;
    }

    return *this;
}

TlSerializeData TlMsgPack::getSerializeData() const
{
    return this->data_;
}

bool TlMsgPack::load(const std::string& path)
{
    std::ifstream ifs;
    ifs.open(path.c_str(), std::ios::in | std::ios::binary);
    if (!ifs) {
        return false;
    }

    this->debugCurrentPos_ = 0; // initialize
    this->data_ = this->loadBinary(ifs);
    return true;
}

void TlMsgPack::pack(const std::string& str)
{
    std::istringstream iss(str);
    this->debugCurrentPos_ = 0; // initialize
    this->data_ = this->loadBinary(iss);
}

TlSerializeData TlMsgPack::loadBinary(std::istream& ifs)
{
    TlSerializeData ans;

    if (ifs.eof() != true) {
        unsigned char c;
        ifs.read((char*)&c, 1);
        ++(this->debugCurrentPos_);
        
        switch (c) {
        case (unsigned char)(0xc0):
            // do nothing
            break;
        
        case (unsigned char)(0xc2):
            ans.set(false);
            break;
        
        case (unsigned char)(0xc3):
            ans.set(true);
            break;

        case (unsigned char)(0xca):
            ans = this->unpack_float(ifs);
            break;

        case (unsigned char)(0xcb):
            ans = this->unpack_double(ifs);
            break;

        case (unsigned char)(0xcc):
            ans = this->unpack_uint8(ifs);
            break;

        case (unsigned char)(0xcd):
            ans = this->unpack_uint16(ifs);
            break;

        case (unsigned char)(0xce):
            ans = this->unpack_uint32(ifs);
            break;

        case (unsigned char)(0xcf):
            ans = (unsigned long)this->unpack_uint64(ifs);
            break;

        case (unsigned char)(0xd0):
            ans = this->unpack_int8(ifs);
            break;

        case (unsigned char)(0xd1):
            ans = this->unpack_int16(ifs);
            break;

        case (unsigned char)(0xd2):
            ans = this->unpack_int32(ifs);
            break;

        case (unsigned char)(0xd3):
            ans = (long)this->unpack_int64(ifs);
            break;
        
        case (unsigned char)(0xda):
            ans = this->unpack_raw16(ifs);
            break;

        case (unsigned char)(0xdb):
            ans = this->unpack_raw32(ifs);
            break;

        case (unsigned char)(0xdc):
            ans = this->unpack_array16(ifs);
            break;

        case (unsigned char)(0xdd):
            ans = this->unpack_array32(ifs);
            break;

        case (unsigned char)(0xde):
            ans = this->unpack_map16(ifs);
            break;

        case (unsigned char)(0xdf):
            ans = this->unpack_map32(ifs);
            break;

        default:
            if (c <= (unsigned char)(0x7f)) {
                ans = this->unpack_positiveFixNum(c);
            } else if (((unsigned char)(0xe0) <= c) && (c <= (unsigned char)(0xff))) {
                ans = this->unpack_negativeFixNum(c);
            } else if (((unsigned char)(0xa0) <= c) && (c <= (unsigned char)(0xbf))) {
                ans = this->unpack_fixraw(c, ifs);
            } else if (((unsigned char)(0x90) <= c) && (c <= (unsigned char)(0x9f))) {
                ans = this->unpack_fixarray(c, ifs);
            } else if (((unsigned char)(0x80) <= c) && (c <= (unsigned char)(0x8f))) {
                ans = this->unpack_fixmap(c, ifs);
            } else {
                std::cerr << TlUtils::format("unknown id=%2x @ %ld",
                                             c,
                                             this->debugCurrentPos_)
                          << std::endl;
                abort();
            }
            break;
        }
    }

    return ans;
}


int TlMsgPack::unpack_positiveFixNum(unsigned char c)
{
    return (c & 127);
}


int TlMsgPack::unpack_negativeFixNum(unsigned char c)
{
    return -(c & 31);
}


TlMsgPack::UINT8 TlMsgPack::unpack_uint8(std::istream& ifs)
{
    unsigned char value;
    ifs.read((char*)&value, 1);
    ++(this->debugCurrentPos_);
    
    return value;
}

TlMsgPack::UINT16 TlMsgPack::unpack_uint16(std::istream& ifs)
{
    union {
        char buf[2];
        UINT16 value;
    } u;
    ifs.read(u.buf, 2);
    this->debugCurrentPos_ += 2;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 2);
    }

    return u.value;
}

TlMsgPack::UINT32 TlMsgPack::unpack_uint32(std::istream& ifs)
{
    union {
        char buf[4];
        UINT32 value;
    } u;
    ifs.read(u.buf, 4);
    this->debugCurrentPos_ += 4;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 4);
    }
    
    return u.value;
}

TlMsgPack::UINT64 TlMsgPack::unpack_uint64(std::istream& ifs)
{
    union {
        char buf[8];
        UINT64 value;
    } u;
    ifs.read(u.buf, 8);
    this->debugCurrentPos_ += 8;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 8);
    }
    
    return u.value;
}

TlMsgPack::INT8 TlMsgPack::unpack_int8(std::istream& ifs)
{
    char value;
    ifs.read((char*)&value, 1);
    ++(this->debugCurrentPos_);

    return value;
}

TlMsgPack::INT16 TlMsgPack::unpack_int16(std::istream& ifs)
{
    union {
        char buf[2];
        INT16 value;
    } u;
    ifs.read(u.buf, 2);
    this->debugCurrentPos_ += 2;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 2);
    }
    
    return u.value;
}

TlMsgPack::INT32 TlMsgPack::unpack_int32(std::istream& ifs)
{
    union {
        char buf[4];
        INT32 value;
    } u;
    ifs.read(u.buf, 4);
    this->debugCurrentPos_ += 4;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 4);
    }
    
    return u.value;
}

TlMsgPack::INT64 TlMsgPack::unpack_int64(std::istream& ifs)
{
    union {
        char buf[8];
        INT64 value;
    } u;
    ifs.read(u.buf, 8);
    this->debugCurrentPos_ += 8;

    if (TlUtils::isLittleEndian() == true) {
        TlUtils::changeEndian(&(u.buf[0]), 8);
    }
    
    return u.value;
}

float TlMsgPack::unpack_float(std::istream& ifs)
{
    assert(sizeof(float) == 4);
    float value;
    ifs.read((char*)&value, 4);
    this->debugCurrentPos_ += 4;

    if (TlUtils::isLittleEndian() == true) {
        value = TlUtils::changeEndian<float>(value);
    }

    return value;
}

double TlMsgPack::unpack_double(std::istream& ifs)
{
    assert(sizeof(double) == 8);
    double value;
    ifs.read((char*)&value, 8);
    this->debugCurrentPos_ += 8;

    if (TlUtils::isLittleEndian() == true) {
        value = TlUtils::changeEndian<double>(value);
    }

    return value;
}

TlSerializeData TlMsgPack::unpack_fixraw(const char in, std::istream& ifs)
{
    const std::size_t size = (in & 31);

    TlSerializeData ans;
    if (size != 0) {
        char* pBuf = new char[size];
        ifs.read(pBuf, size);
        this->debugCurrentPos_ += size;

        const std::string str(pBuf, size);
        ans.set(str);
        
        delete[] pBuf;
        pBuf = NULL;
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_raw16(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint16(ifs);

    TlSerializeData ans;
    if (size != 0) {
        char* pBuf = new char[size];
        ifs.read(pBuf, size);
        this->debugCurrentPos_ += size;
        
        const std::string str(pBuf, size);
        ans.set(str);
        
        delete[] pBuf;
        pBuf = NULL;
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_raw32(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint32(ifs);

    TlSerializeData ans;
    if (size != 0) {
        char* pBuf = new char[size];
        ifs.read(pBuf, size);
        this->debugCurrentPos_ += size;
        
        const std::string str(pBuf, size);

        ans.set(str);
        
        delete[] pBuf;
        pBuf = NULL;
    }
    
    return ans;
}

TlSerializeData TlMsgPack::unpack_fixarray(const char in, std::istream& ifs)
{
    const std::size_t size = (in & 15);

    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData tmp = this->loadBinary(ifs);
        ans.pushBack(tmp);
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_array16(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint16(ifs);

    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData tmp = this->loadBinary(ifs);
        ans.pushBack(tmp);
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_array32(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint32(ifs);

    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData tmp = this->loadBinary(ifs);
        ans.pushBack(tmp);
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_fixmap(const char in, std::istream& ifs)
{
    const std::size_t size = (in & 15);

    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData key = this->loadBinary(ifs);
        const TlSerializeData value = this->loadBinary(ifs);
        ans.add(key, value);
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_map16(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint16(ifs);

    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData key = this->loadBinary(ifs);
        const TlSerializeData value = this->loadBinary(ifs);
        ans.add(key, value);
    }

    return ans;
}

TlSerializeData TlMsgPack::unpack_map32(std::istream& ifs)
{
    const std::size_t size = this->unpack_uint32(ifs);
    
    TlSerializeData ans;
    for (std::size_t i = 0; i < size; ++i) {
        const TlSerializeData key = this->loadBinary(ifs);
        const TlSerializeData value = this->loadBinary(ifs);
        ans.add(key, value);
    }

    return ans;
}

void TlMsgPack::save(const std::string& path) const
{
    std::ofstream ofs;
    ofs.open(path.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    ofs << this->dump();
    ofs.close();
}

std::string TlMsgPack::dump() const
{
    return this->dump(this->data_);
}

std::string TlMsgPack::dump(const TlSerializeData& data) const
{
    std::ostringstream os;

    switch (data.getType()) {
    case TlSerializeData::ARRAY:
        os << this->dump_array(data);
        break;

    case TlSerializeData::MAP:
        os << this->dump_map(data);
        break;

    default:
        os << this->dump_scalar(data);
        break;
    }

    return os.str();
}

std::string TlMsgPack::dump_scalar(const TlSerializeData& data) const
{
    std::ostringstream os;

    switch (data.getType()) {
    case TlSerializeData::BOOLEAN:
        os << this->pack(data.getBoolean());
        break;
        
    case TlSerializeData::STRING:
        {
            const std::string str = data.getStr();
            os << this->pack(str.c_str(), str.size());
        }
        break;

    case TlSerializeData::INT:
        {
            int value = data.getInt();
#if COMPILE_VALUE_SIZEOF_INT == 4
	      os << this->pack_int32(value);
#else
	      os << this->pack_int64(value);
#endif
        }
        break;

    case TlSerializeData::LONG:
        {
	    long value = data.getLong();
#if COMPILE_VALUE_SIZEOF_LONG == 4
	      os << this->pack_int32(value);
#else
	      os << this->pack_int64(value);
#endif
        }
        break;

    case TlSerializeData::UINT:
        {
            unsigned int value = data.getUInt();
#if COMPILE_VALUE_SIZEOF_INT == 4
	      os << this->pack_uint32(value);
#else
	      os << this->pack_uint64(value);
#endif
        }
        break;

    case TlSerializeData::ULONG:
        {
            unsigned long value = data.getULong();
#if COMPILE_VALUE_SIZEOF_LONG == 4
	      os << this->pack_uint32(value);
#else
	      os << this->pack_uint64(value);
#endif
	}
	break;

    case TlSerializeData::DOUBLE:
        {
            double value = data.getDouble();
            os << this->pack(value);
        }
        break;

    case TlSerializeData::NONE:
        {
            std::ostringstream os_tmp;
            this->write(os_tmp, char(0xc0));
            os << os_tmp.str();
        }
        break;

    default:
        // something wrong.
        abort();
        break;
    }

    return os.str();
}

std::string TlMsgPack::dump_array(const TlSerializeData& data) const
{
    assert(data.getType() == TlSerializeData::ARRAY);
    std::ostringstream os;

    assert(sizeof(int) == 4);
    const int size = data.getSize();
    this->write(os, char(0xdd));
    this->write(os, TlUtils::toBigEndian(size));

    for (TlSerializeData::ArrayConstIterator p = data.beginArray(); p != data.endArray(); ++p) {
        os << this->dump(*p);
    }

    return os.str();
}

std::string TlMsgPack::dump_map(const TlSerializeData& data) const
{
    assert(data.getType() == TlSerializeData::MAP);
    std::ostringstream os;

    assert(sizeof(int) == 4);
    const int size = data.getSize();
    this->write(os, char(0xdf));
    this->write(os, TlUtils::toBigEndian(size));

    for (TlSerializeData::MapConstIterator p = data.beginMap(); p != data.endMap(); ++p) {
        os << this->dump(p->first);
        os << this->dump(p->second);
    }

    return os.str();
}

std::string TlMsgPack::pack(const bool value) const
{
    std::ostringstream os;
    if (value == true) {
        this->write(os, char(0xc3));
    } else {
        this->write(os, char(0xc2));
    }
    
    return os.str();
}

std::string TlMsgPack::pack(const UINT8 value) const
{
    std::ostringstream os;

    this->write(os, char(0xcc));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack(const UINT16 value) const
{
    std::ostringstream os;

    this->write(os, char(0xcd));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack_uint32(const UINT32 value) const
{
    std::ostringstream os;

    this->write(os, char(0xce));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack_uint64(const UINT64 value) const
{
    std::ostringstream os;

    this->write(os, char(0xcf));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack(const INT8 value) const
{
    std::ostringstream os;

    this->write(os, char(0xd0));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack(const INT16 value) const
{
    std::ostringstream os;

    this->write(os, char(0xd1));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack_int32(const INT32 value) const
{
    std::ostringstream os;

    this->write(os, char(0xd2));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack_int64(const INT64 value) const
{
    std::ostringstream os;

    this->write(os, char(0xd3));
    this->write(os, TlUtils::toBigEndian(value));

    return os.str();
}

std::string TlMsgPack::pack(const double value) const
{
    assert(sizeof(double) == 8);
    std::ostringstream os;

    this->write(os, char(0xcb));
    this->write(os, TlUtils::toBigEndian(value));
    
    return os.str();
}

std::string TlMsgPack::pack(const char* pBuf, const int size) const
{
    assert(pBuf != NULL);
    std::ostringstream os;

    assert(sizeof(int) == 4);
    this->write(os, char(0xdb)); // raw32
    this->write(os, TlUtils::toBigEndian(size));
    os.write(pBuf, size);

    return os.str();
}
