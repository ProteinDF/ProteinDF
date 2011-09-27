#include <cassert>
#include <iostream>
#include <fstream>
#include <stack>

#include "TlPseudoYaml.h"
#include "TlUtils.h"

TlPseudoYaml::TlPseudoYaml(const TlSerializeData& data)
        : data_(data)
{
}

TlPseudoYaml::TlPseudoYaml(const TlPseudoYaml& rhs)
        : data_(rhs.data_)
{
}

TlPseudoYaml::~TlPseudoYaml()
{
}

TlPseudoYaml& TlPseudoYaml::operator=(const TlPseudoYaml& rhs)
{
    if (this != &rhs) {
        this->data_ = rhs.data_;
    }

    return *this;
}

std::string TlPseudoYaml::str() const
{
    return this->str(this->data_, 0);
}


std::string TlPseudoYaml::str(const TlSerializeData& data, const int level) const
{
    std::string ans = "";

    switch (data.getType()) {
    case TlSerializeData::ARRAY:
        ans += this->str_array(data, level);
        break;

    case TlSerializeData::MAP:
        ans += this->str_map(data, level);
        break;

    default: {
        ans += data.getStr() + "\n";
    }
        break;
    }
    
    return ans;
}

std::string TlPseudoYaml::str_array(const TlSerializeData& data, const int level) const
{
    const std::string indent = TlUtils::repeat("  ", level);
    std::string ans = "";

    for (TlSerializeData::ArrayConstIterator p = data.beginArray(); p != data.endArray(); ++p) {
        const TlSerializeData::DataType childType = p->getType();
        ans += indent + "- ";

        if ((childType == TlSerializeData::ARRAY) || (childType == TlSerializeData::MAP)) {
            ans += "\n";
        }
        const std::string contents = this->str(*p, level +1);
        ans += contents;
    }

    return ans;
}


std::string TlPseudoYaml::str_map(const TlSerializeData& data, const int level) const
{
    const std::string indent = TlUtils::repeat("  ", level);
    std::string ans = "";

    for (TlSerializeData::MapConstIterator p = data.beginMap(); p != data.endMap(); ++p) {
        const TlSerializeData::DataType childType = p->second.getType();
        std::string key = p->first.getStr();
        //p->first.get(&key);
        ans += indent + key + ": ";

        if ((childType == TlSerializeData::ARRAY) || (childType == TlSerializeData::MAP)) {
            ans += "\n";
        }
        const std::string contents = this->str(p->second, level +1);
        ans += contents;
    }

    return ans;
}


TlSerializeData TlPseudoYaml::getSerializeObject() const
{
    return this->data_;
}

void TlPseudoYaml::load(const std::string& path)
{
    TlSerializeData rootObj;

    std::ifstream ifs;
    ifs.open(path.c_str(), std::ios::in);
    //std::cerr << "loading...: " << path << std::endl;

    std::stack<parentInfo> parentStack;
    const parentInfo root(0, &rootObj);
    parentStack.push(root);

    //int indentLevel = 0;
    TlSerializeData* pCurrent = &rootObj;
    while (ifs.eof() == false) {
        std::string line = "";
        std::getline(ifs, line);

        const int numOfIndentSpace = this->getIndentSpace(line);
        const std::string content = line.substr(numOfIndentSpace);

        while (parentStack.size() > 0) {
            if (numOfIndentSpace >= parentStack.top().numOfIndentSpace) {
                // 必要なインデントと同じもしくは大きいので、
                // 親の変更は無し。
                break;
            } else {
                // インデントが小さくなったので、
                // 現在の親(pCurrent)をスタック上の親に登録する。
                assert(numOfIndentSpace < parentStack.top().numOfIndentSpace);
                const std::string key = parentStack.top().key;

                TlSerializeData* parent = parentStack.top().pObj;
                if (key.empty() == true) {
                    parent->pushBack(*pCurrent);
                } else {
                    parent->add(TlSerializeData(key), *pCurrent);
                }

                // 現在の親を、その更に親に変更する。
                delete pCurrent;
                pCurrent = parent;

                parentStack.pop();
            }
        }
        assert(pCurrent != NULL);

        std::string key = "";
        std::string value = "";
        ContentType type = this->judgeContentType(content, &key, &value);
        switch (type) {
        case ARRAY_STR:
            pCurrent->pushBack(TlSerializeData(value));
            break;

        case ARRAY_CONTAINER: {
            // 現在の親をスタックに乗せる。
            parentInfo tmp(numOfIndentSpace +1, pCurrent, "");
            parentStack.push(tmp);
            // 新たな親をcurrentにする。
            TlSerializeData* pObj = new TlSerializeData;
            pCurrent = pObj;
        }
            break;

        case HASH_STR:
            pCurrent->add(TlSerializeData(key), TlSerializeData(value));
            break;

        case HASH_CONTAINER: {
            // 現在の親をスタックに乗せる。
            parentInfo tmp(numOfIndentSpace +1, pCurrent, key);
            parentStack.push(tmp);
            // 新たな親をcurrentにする。
            TlSerializeData* pObj = new TlSerializeData;
            pCurrent = pObj;
        }
            break;
            
        default:
            break;
        }
    }

    ifs.close();

    // finalize
    while (parentStack.size() > 1) {
        const std::string key = parentStack.top().key;
        TlSerializeData* parent = parentStack.top().pObj;
        if (key.empty() == true) {
            parent->pushBack(*pCurrent);
        } else {
            parent->add(TlSerializeData(key), *pCurrent);
        }
        pCurrent = parent;
        parentStack.pop();
    }
    this->data_ = rootObj;
}


std::size_t TlPseudoYaml::getIndentSpace(const std::string& str)
{
    std::size_t ans = 0;
    const std::size_t length = str.size();
    while ((ans < length) && (str[ans] == ' ')) {
        ++ans;
    }

    return ans;
}


TlPseudoYaml::ContentType TlPseudoYaml::judgeContentType(const std::string& str,
                                                         std::string* pKey, std::string* pValue)
{
    assert(pKey != NULL);
    assert(pValue != NULL);

    ContentType type = UNKNOWN;

    if (str[0] == '-') {
        type = ARRAY_CONTAINER;

        if (str.size() > 1) {
            std::string value = str.substr(1);
            TlUtils::rtrim_ws(value);

            if (value.size() > 0) {
                type = ARRAY_STR;
                *pValue = value;
            }
        }
    } else {
        std::string::size_type pos = str.find(':');
        if (pos != std::string::npos) {
            std::string key = str.substr(0, pos);
            TlUtils::rtrim_ws(key);
            std::string value = str.substr(pos +1);
            TlUtils::trim_ws(value);
            TlUtils::rtrim_ws(value);

            if (key.empty() == false) {
                if (value.empty() == true) {
                    type = HASH_CONTAINER;
                    *pKey = key;
                } else {
                    type = HASH_STR;
                    *pKey = key;
                    *pValue = value;
                }
            }
        }
    }

    return type;
}





