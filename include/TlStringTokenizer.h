#ifndef TLSTRINGTOKENIZER_H
#define TLSTRINGTOKENIZER_H

#include <string>
#include <iostream>

/// std::string オブジェクトをトークンに分割するクラス
class TlStringTokenizer {
public:
    /// コンストラクタ
    ///
    /// @param[in] str 分割される文字列
    /// @param[in] delimiter デリミタとなる文字列. 空文字列の場合はホワイトスペースが使われる.
    explicit TlStringTokenizer(const std::string& str, const std::string& delimiter = "");

    /// デストラクタ
    ~TlStringTokenizer();

public:
    /// トークンの数を返す
    ///
    /// @return トークンの数
    std::size_t countTokens();

    /// 次のトークンの有無を返す
    ///
    /// @retval true 次のトークンあり
    /// @retval false 次のトークンなし
    bool hasMoreTokens();

    /// 次のトークンを返す
    ///
    /// @return 次のトークン. 無い場合は空文字列を返す.
    std::string nextToken();

private:
    std::string m_delimiter; /// デリミタを保持する
    std::string m_str;       /// 分割される文字列
    int m_count;             /// トークンの数
    std::string::size_type m_begin; /// トークン処理を行う始めの文字カウント
    std::string::size_type m_end;   /// トークン処理を行う最後の文字カウント
};

#endif // TLSTRINGTOKENIZER_H
