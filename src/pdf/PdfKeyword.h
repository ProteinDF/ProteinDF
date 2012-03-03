#ifndef PDFKEYWORD_H
#define PDFKEYWORD_H

#include <string>
#include <map>
#include "TlUtils.h"
#include "TlParameter.h"
#include "TlSerializeData.h"

#define PDFKWD_PRINT_WIDTH (72)

/// ProteinDFのキーワード群を保持するクラス
class PdfKeyword {
public:
    /// キーワードタイプに使われる値
    ///
    /// maskに利用されるため、定数は2の累乗であること
    enum KeywordType {
        KWD_DEFAULT  = 0, // デフォルト値が設定されるキーワード
        KWD_HIDDEN   = 1,
        KWD_INTERNAL = 2, // 内部で利用されるキーワード
        KWD_DEBUG    = 4
    };
    
    /// キーワード情報を保持する構造体
    struct KeywordInfo {
    public:
        std::string keyword;       /// キーワード
        unsigned int type;         /// キーワードタイプ
        std::string explanation;   /// 説明(英語)
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
            std::string explanation = "      " + TlUtils::textWrap(this->kwdList_[i].explanation, 72);
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
            // if ((line.length() + value.length()) < PDFKWD_PRINT_WIDTH) {
            //     TlUtils::pad(line, (PDFKWD_PRINT_WIDTH - value.length()), ' ');
            // }
            line += value;
            
            out << line << "\n";
            ++index;
        }
    }
    out.flush();
}


#endif // PDFKEYWORD_H
