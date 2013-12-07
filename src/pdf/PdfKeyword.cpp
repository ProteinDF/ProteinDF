#include <iostream>
#include <cassert>
#include "PdfKeyword.h"
#include "TlUtils.h"
#include "TlMsgPack.h"
#include "TlLogging.h"

PdfKeyword::PdfKeyword()
 : log_(TlLogging::getInstance()) {
    this->initialize();
}


PdfKeyword::~PdfKeyword()
{
}


TlSerializeData PdfKeyword::getDefault() const
{
    TlSerializeData data;

    KeywordListType::const_iterator pEnd = this->kwdList_.end();
    for (KeywordListType::const_iterator p = this->kwdList_.begin(); p != pEnd; ++p) {
        const std::string keyword = p->keyword;
        const std::string value = p->defaultValue;
        data[keyword] = value;
    }

    return data;
}


void PdfKeyword::initialize()
{
    this->kwdList_.clear();
    KeywordInfo item;

    std::string kwdPath = "pdfkwd.mpac";
    const char* PdfHome = std::getenv("PDF_HOME");
    if (PdfHome != NULL) {
        kwdPath = TlUtils::format("%s/data/%s", PdfHome, kwdPath.c_str());
    }

    TlMsgPack mpac;
    if (mpac.load(kwdPath) == true) {
        const TlSerializeData kwds = mpac.getSerializeData();
        TlSerializeData::ArrayConstIterator itEnd = kwds.endArray();
        for (TlSerializeData::ArrayConstIterator it = kwds.beginArray(); it != itEnd; ++it) {
            KeywordInfo item;
            item.keyword      = (*it)["keyword"].getStr();
            item.explanation  = (*it)["property"]["description"].getStr();
            item.defaultValue = (*it)["property"]["default"].getStr();
            item.syntax       = (*it)["property"]["syntax"].getStr();
            if ((*it)["property"]["hidden"].getBoolean() == true) {
                item.type = KWD_INTERNAL | KWD_DEBUG;
            } else {
                item.type = KWD_DEFAULT;
            }

            this->kwdList_.push_back(item);
        }
    }

    // alias
    this->setupAliasList();
}


void PdfKeyword::setupAliasList()
{
    std::string aliasKwdPath = "alias.mpac";
    const char* PdfHome = std::getenv("PDF_HOME");
    if (PdfHome != NULL) {
        aliasKwdPath = TlUtils::format("%s/data/%s", PdfHome, aliasKwdPath.c_str());
    }

    TlMsgPack mpac;
    if (mpac.load(aliasKwdPath) == true) {
        const TlSerializeData kwds = mpac.getSerializeData();
        TlSerializeData::ArrayConstIterator itEnd = kwds.endArray();
        for (TlSerializeData::ArrayConstIterator it = kwds.beginArray(); it != itEnd; ++it) {
            std::string oldKwd = (*it)["old_kwd"].getStr();
            std::string newKwd = (*it)["new_kwd"].getStr();

            this->kwdAlias_[oldKwd] = newKwd;
        }
    }
}


void PdfKeyword::convertAlias(TlSerializeData* pData)
{
    assert(pData != NULL);
    TlLogging& log = TlLogging::getInstance();
    
    TlSerializeData::MapIterator p = pData->beginMap();
    TlSerializeData::MapIterator pEnd = pData->endMap();
    while (p != pEnd) {
        const std::string keyword = p->first.getStr();
        const TlSerializeData value = p->second;

        AliasContainerType::const_iterator it = this->kwdAlias_.find(keyword);
        if (it != this->kwdAlias_.end()) {
            const std::string newKeyword = it->second;

            if (pData->hasKey(newKeyword) == true) {
                log.warn(TlUtils::format(" kwd[%s] overwrite by alias[%s]:",
                                         newKeyword.c_str(),
                                         keyword.c_str()));
                log.warn(TlUtils::format("  from [%s] to [%s]",
                                         (*pData)[newKeyword].getStr().c_str(),
                                         value.getStr().c_str()));
            }

            (*pData)[newKeyword] = value;
            pData->erase(keyword);
            p = pData->beginMap();
        } else {
            ++p;
        }
    }
}


std::string PdfKeyword::getCSV(bool showHiddenItem) const
{
    std::string output = "keyword, explanation, default, syntax\n";

    const int numOfItems = this->kwdList_.size();
    for (int i = 0; i < numOfItems; ++i) {
        const KeywordInfo& item = this->kwdList_[i];
        if (((item.type & KWD_HIDDEN) != 0) &&
            (showHiddenItem == false)) {
            continue;
        }
        output += TlUtils::format("\"%s\", \"%s\", \"%s\", \"%s\"\n",
                                  item.keyword.c_str(),
                                  item.explanation.c_str(),
                                  item.defaultValue.c_str(),
                                  item.syntax.c_str());
    }

    return output;
}

std::string PdfKeyword::getCSV_jp(bool showHiddenItem) const
{
    std::string output = "keyword, explanation, default, syntax\n";
    
    const int numOfItems = this->kwdList_.size();
    for (int i = 0; i < numOfItems; ++i) {
        const KeywordInfo& item = this->kwdList_[i];
        if (((item.type & KWD_HIDDEN) != 0) &&
            (showHiddenItem == false)) {
            continue;
        }
        output += TlUtils::format("\"%s\", \"%s\", \"%s\", \"%s\"\n",
                                  item.keyword.c_str(),
                                  item.desc_jp.c_str(),
                                  item.defaultValue.c_str(),
                                  item.syntax.c_str());
    }
    
    return output;
}

std::string PdfKeyword::get_reST_jp(bool showHiddenItem) const
{
    std::string output = "";
    
    const int numOfItems = this->kwdList_.size();
    for (int i = 0; i < numOfItems; ++i) {
        const KeywordInfo& item = this->kwdList_[i];
        if (((item.type & KWD_HIDDEN) != 0) &&
            (showHiddenItem == false)) {
            continue;
        }
        
        output += TlUtils::format("* %s\n\n", item.keyword.c_str());
        output += ".. csv-table::\n";
        output += "  :widths: 20,80\n";
        output += "  :stub-columns: 1\n";
        output += "\n";
        output += TlUtils::format("  \"パラメータ\", \"%s\"\n", item.keyword.c_str());
        output += TlUtils::format("  \"説明\", \"%s\"\n", item.desc_jp.c_str());
        output += TlUtils::format("  \"規定値\", \"%s\"\n", item.defaultValue.c_str());
        output += TlUtils::format("  \"取りうる値\", \"%s\"\n", item.syntax.c_str());
        output += "\n\n";
    }
    
    return output;
}

TlSerializeData PdfKeyword::getSerializeData() const
{
    TlSerializeData data;

    KeywordListType::const_iterator itEnd = this->kwdList_.end();
    for (KeywordListType::const_iterator it = this->kwdList_.begin(); it != itEnd; ++it) {
        TlSerializeData property;
        property["description"] = it->explanation;
        //property["description_jp"] = it->explanationJ;
        property["default"] = it->defaultValue;
        property["syntax"] = it->syntax;
        property["hidden"] = (it->type != KWD_DEFAULT);

        TlSerializeData item;
        item["keyword"] = it->keyword;
        item["property"] = property;
        
        data.pushBack(item);
    }
    
    return data;
}


void PdfKeyword::checkInputParam(const TlSerializeData& param) const
{
    if (param.getType() != TlSerializeData::MAP) {
        std::cerr << TlUtils::format("[WARN] input param type mismatch. file: %s, line: %d",
                                     __FILE__, __LINE__)
                  << std::endl;
    }

    TlSerializeData::MapConstIterator mapItEnd = param.endMap();
    for (TlSerializeData::MapConstIterator mapIt = param.beginMap(); mapIt != mapItEnd; ++mapIt) {
        const std::string keyword = mapIt->first.getStr();
        const bool hasKeyword = this->hasKeyword(keyword);
        if (hasKeyword == false) {
            std::cerr << TlUtils::format("[WARN] unknown keyword: %s", keyword.c_str())
                      << std::endl;
        }
    }
}


bool PdfKeyword::hasKeyword(const std::string& keyword) const
{
    if (this->kwdDb_.empty() == true) {
        this->makeDB();
    }

    KeywordDbType::const_iterator it = this->kwdDb_.find(keyword);
    return (it != this->kwdDb_.end());
}


void PdfKeyword::makeDB() const
{
    this->kwdDb_.clear();

    const std::size_t numOfKwds = this->kwdList_.size();
    for (std::size_t i = 0; i < numOfKwds; ++i) {
        const KeywordInfo& info = this->kwdList_[i];
        const std::string& keyword = info.keyword;

        this->kwdDb_[keyword] = info;
    }

    if (this->kwdDb_.size() != numOfKwds) {
        this->log_.debug(TlUtils::format("[WARN] size mismatch: file: %s, line: %d.", __FILE__, (int)__LINE__));
    }
}

