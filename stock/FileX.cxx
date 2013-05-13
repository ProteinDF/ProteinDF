#include <cassert>
#include  "FileX.h"


FileX::FileX(const std::string& rFileName, std::ios_base::openmode mode)
        : m_sFileName(rFileName), m_mode(mode)
{
}

FileX::~FileX()
{
}

/**
 * カレントディレクトリにsFileNameが存在していればtrue、なければfalseを返す。
 */
bool FileX::isExist(const std::string& sFileName)
{
    bool bAnswer = true;
    std::ifstream fs;
    fs.open(sFileName.c_str());
    if (fs.is_open()) {
        bAnswer = true;
    } else {
        bAnswer = false;
    }
    fs.close();

    return bAnswer;
}

void FileX::open()
{
    this->m_FileStream.open(this->m_sFileName.c_str(), this->m_mode);
}

void FileX::close()
{
    if (this->isOpen()) {
        this->m_FileStream.close();
    }
}

bool FileX::isOpen()
{
    return this->m_FileStream.is_open();
}

// バッファに残っているデータを強制的に書き込みます.
void FileX::flush()
{
    this->m_FileStream.flush();
};

int FileX::width(int w)
{
    return this->m_FileStream.width(w);
}

int FileX::width() const
{
    return this->m_FileStream.width();
}

void FileX::unsetf(std::ios_base::fmtflags f)
{
    return this->m_FileStream.unsetf(f);
}

std::ios_base::fmtflags FileX::setf(std::ios_base::fmtflags f,
                                    std::ios_base::fmtflags m)
{
    return this->m_FileStream.setf(f, m);
}

std::ios_base::fmtflags FileX::setf(std::ios_base::fmtflags f)
{
    return this->m_FileStream.setf(f);
}

std::streamsize FileX::precision(std::streamsize n)
{
    return this->m_FileStream.precision(n);
}

std::streamsize FileX::precision() const
{
    return this->m_FileStream.precision();
}

std::ios_base::fmtflags FileX::flags(std::ios_base::fmtflags f)
{
    return this->m_FileStream.flags(f);
}

std::ios_base::fmtflags FileX::flags()
{
    return this->m_FileStream.flags();
}

char FileX::fill(char c)
{
    return this->m_FileStream.fill(c);
}

char FileX::fill() const
{
    return this->m_FileStream.fill();
}
