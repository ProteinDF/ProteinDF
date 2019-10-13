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

#ifndef TLPSEUDOYAML_H
#define TLPSEUDOYAML_H

#include <map>
#include <string>
#include <vector>

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
        explicit parentInfo(int n = 0, TlSerializeData* p = NULL,
                            const std::string& k = "")
            : numOfIndentSpace(n), pObj(p), key(k){};

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
    ContentType judgeContentType(const std::string& str, std::string* pKey,
                                 std::string* pValue);

   protected:
    TlSerializeData data_;
};

#endif  // TLPSEUDOYAML_H
