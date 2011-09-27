#ifndef FILEX_H
#define FILEX_H

#include <ios>
#include <iostream>
#include <fstream>
#include <string>

class FileX {
public:
    explicit FileX(const std::string& rFileName = "",
                   std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out);
    virtual ~FileX();

public:
    // カレントディレクトリ以下にファイルsFileNameが存在していればtrue、なければfalseを返す。
    static bool isExist(const std::string& sFileName);

public:
    virtual void open();
    virtual void close();
    bool isOpen();

public:
    // バッファに残っているデータを強制的に書き込みます.
    void flush();

    int width(int w);
    int width() const;

    void unsetf(std::ios_base::fmtflags f);
    std::ios_base::fmtflags setf(std::ios_base::fmtflags f, std::ios_base::fmtflags m);
    std::ios_base::fmtflags setf(std::ios_base::fmtflags f);
    std::streamsize precision(std::streamsize n);
    std::streamsize precision() const;
    std::ios_base::fmtflags flags(std::ios_base::fmtflags f);
    std::ios_base::fmtflags flags();
    char fill(char c);
    char fill() const;

public:
    // operator<<() for output
    template<typename T>
    friend FileX& operator<<(FileX& lhs, const T& t) {
        lhs.m_FileStream << t;
        return lhs;
    }

    // std::endlが使えるようにするために必要なメンバ関数
    void operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&)) {
        m_FileStream << (*pf);
    }

    // operator>>() for input
    template<typename T>
    friend FileX& operator>>(FileX& lhs, T& t) {
        lhs.m_FileStream >> t;
        return lhs;
    }

    // member variables
protected:
    //std::string m_sCurrentAbsoluteDirectoryPath;
    std::string m_sFileName;
    std::ios_base::openmode m_mode;

    std::fstream m_FileStream;
};

#endif // FILEX_H
