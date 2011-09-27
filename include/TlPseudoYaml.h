#ifndef TLPSEUDOYAML_H
#define TLPSEUDOYAML_H

#include <string>
#include <vector>
#include <map>

#include "TlSerializeData.h"

class TlPseudoYaml {
public:
    enum ContentType {
        UNKNOWN,
        ARRAY_STR,
        ARRAY_CONTAINER,
        HASH_STR,
        HASH_CONTAINER
    };

    struct parentInfo {
    public:
        explicit parentInfo(int n = 0, TlSerializeData* p = NULL, const std::string& k = "") 
            : numOfIndentSpace(n), pObj(p), key(k) {
    };
        
    public:
        int numOfIndentSpace;
        TlSerializeData* pObj;
        std::string key;
    };


public:
    explicit TlPseudoYaml(const TlSerializeData& data = TlSerializeData());
    TlPseudoYaml(const TlPseudoYaml& rhs);
    ~TlPseudoYaml();
    
public:
    TlPseudoYaml& operator=(const TlPseudoYaml& rhs);

    TlSerializeData getSerializeObject() const;
    void load(const std::string& path);
    std::string str() const;
    
protected:
    std::string str(const TlSerializeData& data, const int level) const;
    std::string str_array(const TlSerializeData& data, const int level) const;
    std::string str_map(const TlSerializeData& data, const int level) const;
    
    std::size_t getIndentSpace(const std::string& str);
    ContentType judgeContentType(const std::string& str,
                                 std::string* pKey, std::string* pValue);

protected:
    TlSerializeData data_;
};

#endif // TLPSEUDOYAML_H
