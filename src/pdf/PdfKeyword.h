#ifndef PDFKEYWORD_H
#define PDFKEYWORD_H

#include <string>
#include <map>
#include "TlUtils.h"
#include "TlParameter.h"
#include "TlSerializeData.h"

/// ProteinDFのキーワード群を保持するクラス
class PdfKeyword {
public:
    /// キーワードタイプに使われる値
    ///
    /// maskに利用されるため、定数は2の累乗であること
    enum KeywordType {
        KWD_DEFAULT  = 0,
        KWD_DEBUG    = 1,
        KWD_HIDDEN   = 2
    };
    
    /// キーワード情報を保持する構造体
    struct KeywordInfo {
    public:
        std::string keyword;       /// キーワード
        KeywordType type;         /// キーワードタイプ
        std::string explanation;   /// 説明(英語)
        std::string explanationJ;  /// 説明(日本語)
        std::string defaultValue;  /// デフォルト値
        std::string syntax;        /// 文法
    };

public:
    typedef std::vector<KeywordInfo> KeywordListType;
    typedef std::map<std::string, std::string> AliasContainerType;
    typedef std::map<std::string, KeywordInfo> KeywordDbType;
    
public:
    /// デフォルトコンストラクタ
    PdfKeyword();

    /// デストラクタ
    ~PdfKeyword();

public:
    /// デフォルト値を設定する
    TlSerializeData getDefault() const;

    /// 入力パラメータのキーワードが含まれているかをチェックする
    void checkInputParam(const TlSerializeData& param) const;
    
    /// データの内容をテキスト出力する
    ///
    /// outオブジェクトは <<演算子を定義してある必要がある.
    template<typename T>
    void printDefault(T& out) const;

    //
    void setupAliasList();
    void convertAlias(TlSerializeData* pData);
    
    template<typename T>
    void print(T& out, const TlSerializeData& data) const;

    std::string getCSV(bool showHiddenItem = false) const;
    std::string getCSV_J(bool showHiddenItem = false) const;

    TlSerializeData getSerializeData() const;
    
protected:
    /// 初期化
    void initialize();

    bool hasKeyword(const std::string& keyword) const;
    void makeDB() const;
    
private:
    /// prohibition of constructer
    PdfKeyword(const PdfKeyword&);
    PdfKeyword& operator=(const PdfKeyword&);

protected:
    KeywordListType kwdList_;
    AliasContainerType kwdAlias_;
    mutable KeywordDbType kwdDb_;
};


template<typename T>
void PdfKeyword::printDefault(T& out) const
{
    const TlSerializeData data = this->getDefault();
    const int numOfKeywords = this->kwdList_.size();
    int index = 0;
    for (int i = 0; i < numOfKeywords; ++i) {
        if ((this->kwdList_[i].type & KWD_HIDDEN) == 0) {
            const std::string keyword = this->kwdList_[i].keyword;
                
            // 出力
            const std::string index_str = TlUtils::format("% 3d: ", index +1);
            std::string explanation = "      " + TlUtils::textWrap(this->kwdList_[i].explanation, 74);
            TlUtils::replace(explanation, "\n", "\n      "); // 行の先頭に空白を入れる
            const std::string value = "      default: " + data[keyword].getStr();
            const std::string syntax = "      syntax: " + this->kwdList_[i].syntax;
            
            out << index << " "
                << keyword
                << "\n";
            out << explanation << "\n";
            out << value << "\n";
            out << syntax << "\n";

            ++index;
        }
    }
    out.flush();
}


template<typename T>
void PdfKeyword::print(T& out, const TlSerializeData& data) const
{
    const int numOfKeywords = this->kwdList_.size();
    int index = 0;
    for (int i = 0; i < numOfKeywords; ++i) {
        if ((this->kwdList_[i].type & KWD_HIDDEN) == 0) {
            const std::string keyword = this->kwdList_[i].keyword;

            // 出力
            const std::string index_str = TlUtils::format("% 3d: ", index +1);
            const std::string value = "[" + data[keyword].getStr() + "]";
            
            std::string line = index_str + " " + keyword + " ";
            if ((line.length() + value.length()) < 80) {
                TlUtils::pad(line, (80 - value.length()), ' ');
            }
            line += value;
            
            out << line << "\n";
            ++index;
        }
    }
    out.flush();
}



// template<typename T>
// void PdfKeyword::printGroup(const std::string& sGroupName, T& out) const
// {
//     PdfKeyword::Groups::const_iterator pGroup = this->m_data.find(sGroupName);
//     if (pGroup == this->m_data.end()) {
//         return;
//     }
//     const PdfKeyword::Keywords& keywords = pGroup->second;

//     int index = 1;
//     for (PdfKeyword::Keywords::const_iterator p = keywords.begin();  p != keywords.end(); ++p) {
//         const KeywordInfo info = p->second;

//         // No. (5)
//         const std::string sIndex = TlUtils::format("% 2d: ", index);

//         // Keyword (--)
//         std::string sKeyword = p->first;
//         TlUtils::pad(sKeyword, 40, ' ');

//         // Rating (8)
//         std::string sRating;
//         switch (info.nRating) {
//         case KWD_SHORT:
//             sRating = "SHORT ";
//             break;
//         case KWD_MEDIUM:
//             sRating = "MEDIUM";
//             break;
//         case KWD_LONG:
//             sRating = "LONG  ";
//             break;
//         default:
//             sRating = "UNDEF ";
//             break;
//         }
//         sRating = "[" + sRating + "]";

//         // Type (10)
//         std::string sType;
//         switch (info.nType) {
//         case KWD_FORMATED:
//             sType = "FORMATED";
//             break;
//         case KWD_CONTROL:
//             sType = "CONTROL ";
//             break;
//         case KWD_JUNCTION:
//             sType = "JUNCTION";
//             break;
//         case KWD_COMPDEF:
//             sType = "COMPDEF ";
//         default:
//             sType = "DEFAULT ";
//         }
//         sType = "[" + sType + "]";

//         // explanation
//         std::string sExplanation = std::string("  ") + TlUtils::textWrap(info.sExplanation, 108);
//         TlUtils::replace(sExplanation, "\n", "\n  "); // 行の先頭に空白を入れる

//         // default
//         std::string sDefaultValue = info.sDefaultValue;

//         // syntax
//         std::string sSyntax = info.sSyntax;

//         // reference
//         std::string sReference = TlUtils::textWrap(info.sParallelKeywords, 106);
//         TlUtils::replace(sReference, "\n", "\n              "); // 行の先頭に空白を入れる

//         // 出力
//         out << sIndex   << " "
//         << sKeyword << " "
//         << sRating  << " "
//         << sType    << "\n";

//         out << " default    : " << sDefaultValue << "\n";
//         out << " syntax     : " << sSyntax       << "\n";
//         if (!sReference.empty()) {
//             out << " reference  : " << sReference    << "\n";
//         }
//         out << sExplanation  << "\n";
//         out << "\n";

//         //
//         ++index;
//     }
// }


// template<typename T>
// void PdfKeyword::print(T& out) const
// {
//     for (Groups::const_iterator p = this->m_data.begin(); p != this->m_data.end(); p++) {
//         const std::string sGroup = p->first;
//         out << "// ---------------------------------------------------------------------\n";
//         out << "// >>>> " << sGroup << "\n";
//         this->printGroup<T>(sGroup, out);
//         out << "\n";
//     }
//     out.flush();
// }

// template<typename T>
// void PdfKeyword::print(const std::string& sGroup, T& out) const
// {
//     Groups::const_iterator p = this->m_data.find(sGroup);
//     if (p != this->m_data.end()) {
//         out << "// ---------------------------------------------------------------------\n";
//         out << "// >>>> " << sGroup << "\n";
//         this->printGroup<T>(sGroup, out);
//         out << "\n";
//     }
//     out.flush();
// }


#endif // PDFKEYWORD_H
