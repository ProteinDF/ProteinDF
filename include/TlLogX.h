#ifndef TLLOGX_H
#define TLLOGX_H

#include <string>
#include <fstream>

class CnTimeX;

/** ログをファイルに出力するクラス
 */
class TlLogFileX {
public:
    /** コンストラクタ
     *
     *  @param[in] sFileName 出力するファイル名
     */
    TlLogFileX(const std::string& sFileName);

    /** デストラクタ
     */
    ~TlLogFileX();

public:
    /** ファイルが開かれているかを返す
     *
     *  @retval true ファイルは開かれている
     *  @retval false ファイルは開かれていない
     */
    bool isOpen();

    /** ファイルを開く
     */
    void open();

    ///
    void setFilePath(const std::string& path);
   
    /** フラッシュする
     */
    void flush();

    /** iostream のような使い方ができる出力用 <<演算子
     */
    template<typename T>
    friend TlLogFileX& operator<<(TlLogFileX& lhs, const T& t) {
        //assert(lhs.isOpen());
        if (!lhs.isOpen()) {
            lhs.open();
        }

        lhs.m_FileStream << t;
        return lhs;
    }

    /** iostream のような使い方ができる出力用 <<演算子
     *
     *  std::endlが使えるように特殊化
     */
    void operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&)) {
        if (!this->isOpen()) {
            this->open();
        }

        m_FileStream << (*pf);
    }

private:
    /// ファイル名を保持する
    std::string m_sFileName;

    /// 出力用ファイルストリーム
    std::fstream m_FileStream;
};


/**
 *  ログクラス
 *  singleton pattern
 */
class TlLogX {
public:
    /** ログレベル
     */
    enum Level {
        DEBUG = 0,
        INFO,
        WARN,
        ERROR,
        FATAL
    };

public:
    static TlLogX& getInstance();
    void setLevel(Level nLevel);

    void setFilePath(const std::string filePath) {
        this->m_stdout.setFilePath(filePath);
    }
    
public:
    void debug(const std::string& str);
    void info(const std::string& str);
    void warn(const std::string& str);
    void error(const std::string& str);
    void fatal(const std::string& str);

    void flush();

public:
    template <typename T>
    TlLogX& operator<<(const T& t);

    // std::endlが使えるようにするために必要なメンバ関数
    void operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&));

public:
    /** セクションタイトル表示
     */
//   void outputStartTitle(const std::string& stepName, const char lineChar = '-');
//   void outputEndTitle(const std::string& stepName ="", const char lineChar = '-');
//   void outputEndTitle(const std::string& stepName, const CnTimeX& rTime, const char lineChar = '-');

private:
    TlLogX();
    ~TlLogX();

    TlLogX(const TlLogX&);
    void operator=(const TlLogX&);

protected:
    Level m_nLevel;

    TlLogFileX m_stdout;
    TlLogFileX m_debug;

private:
    static TlLogX* m_pInstance;
};

////////////////////////////////////////////////////////////////////////
// template
//
template <typename T>
TlLogX& TlLogX::operator<<(const T& t)
{
    if (this->m_nLevel <= TlLogX::INFO) {
        this->m_stdout << t;
        this->m_stdout.flush();
    }

    return (*this);
}

#endif // TLLOGX_H
