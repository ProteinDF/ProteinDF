#include "TlStringTokenizer.h"

TlStringTokenizer::TlStringTokenizer(const std::string& str, const std::string& delimiter)
        : m_str(str), m_count(-1), m_begin(0), m_end(0)
{

    if (delimiter == "") {
        this->m_delimiter = " \f\n\r\t\v"; // デフォルトはホワイトスペース
    } else {
        this->m_delimiter = delimiter;
    }

    // 1つめのトークンをポイントする
    this->m_begin = str.find_first_not_of(this->m_delimiter);
    this->m_end   = str.find_first_of(this->m_delimiter, this->m_begin);
}

TlStringTokenizer::~TlStringTokenizer()
{
}

std::size_t TlStringTokenizer::countTokens()
{
    if (this->m_count >= 0) {
        // 既に計算されているので制御を戻す
        return m_count;
    }

    std::string::size_type n = 0;
    std::string::size_type i = 0;

    for (;;) {
        // 1つめのトークンに進む
        if ((i = this->m_str.find_first_not_of(this->m_delimiter, i)) == std::string::npos) {
            break;
        }

        // 次の区切り文字に進む
        i = this->m_str.find_first_of(this->m_delimiter, i+1);
        n++;

        if (i == std::string::npos) {
            break;
        }
    }

    return (this->m_count = n);
}

bool TlStringTokenizer::hasMoreTokens()
{
    return (this->m_begin != this->m_end);
}

std::string TlStringTokenizer::nextToken()
{
    std::string s = "";

    if (this->m_begin != std::string::npos) {
        if (this->m_end != std::string::npos) {
            s = this->m_str.substr(this->m_begin, (this->m_end - this->m_begin));
            this->m_begin = this->m_str.find_first_not_of(this->m_delimiter, this->m_end);
            this->m_end   = this->m_str.find_first_of(this->m_delimiter, this->m_begin);
        } else {
            s = this->m_str.substr(this->m_begin, (this->m_str.length() - this->m_begin));
            this->m_begin = this->m_str.find_first_not_of(this->m_delimiter, this->m_end);
        }
    }

    return s;
}

