#include <iostream>
#include "TlLogX.h"
#include "TlUtils.h"

TlLogX* TlLogX::m_pInstance = NULL;

TlLogX::TlLogX()
    : m_nLevel(TlLogX::INFO), m_stdout("fl_Out_Std"), m_debug("pdf_debug.txt")
{
}

TlLogX::~TlLogX()
{
    this->flush();
}

TlLogX& TlLogX::getInstance()
{
    if (TlLogX::m_pInstance == NULL) {
        TlLogX::m_pInstance = new TlLogX();
    }

    return *(TlLogX::m_pInstance);
}

void TlLogX::setLevel(TlLogX::Level nLevel)
{
    this->m_nLevel = nLevel;
}

void TlLogX::debug(const std::string& str)
{
    if (this->m_nLevel == TlLogX::DEBUG) {
        this->m_debug << str;
        this->flush();
    }
}

void TlLogX::info(const std::string& str)
{
    if (this->m_nLevel <= TlLogX::INFO) {
        this->m_stdout << str;
        this->flush();
    }
}

void TlLogX::warn(const std::string& str)
{
    if (this->m_nLevel <= TlLogX::WARN) {
        this->m_stdout << str;
        this->flush();
    }
}

void TlLogX::error(const std::string& str)
{
    if (this->m_nLevel <= TlLogX::ERROR) {
        this->m_stdout << str;
        std::cerr << str;
        this->flush();
    }
}

void TlLogX::fatal(const std::string& str)
{
    this->m_stdout << str;
    std::cerr << str;

    this->flush();
}

void TlLogX::flush()
{
    this->m_debug.flush();
    this->m_stdout.flush();
    std::cerr.flush();
    std::cout.flush();
}

void TlLogX::operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&))
{
    if (this->m_nLevel <= TlLogX::INFO) {
        this->m_stdout << *pf;
    }
}

// void TlLogX::outputStartTitle(const std::string& stepName, const char lineChar){
//   std::string line = "";
//   TlUtils::pad(line, 72, lineChar);

//   const std::string timeString = TlUtils::format("[%s %s]", CnTimeX::getNowDate().c_str(), CnTimeX::getNowTime().c_str());

//   std::string title = ">>>> " + stepName + " ";
//   TlUtils::pad(title, (72 - timeString.length()), ' ');
//   title += timeString;

//   // 出力
//   *this << line << "\n";
//   *this << title;
//   *this << "\n";
//   //*this << line << "\n\n";
// }

// void TlLogX::outputEndTitle(const std::string& stepName, const char lineChar){
//   std::string line = "";
//   TlUtils::pad(line, 72, lineChar);

//   const std::string timeString = TlUtils::format("[%s %s]", CnTimeX::getNowDate().c_str(), CnTimeX::getNowTime().c_str());

//   std::string title = "<<<< " + stepName + " ";
//   TlUtils::pad(title, (72 - timeString.length()), ' ');
//   title += timeString;

//   // 出力
//   *this << "\n";
//   *this << title << "\n";
//   *this << line << "\n\n";
// }

// void TlLogX::outputEndTitle(const std::string& stepName, const CnTimeX& rTime, const char lineChar){
//   std::string line = "";
//   TlUtils::pad(line, 72, lineChar);

//   const std::string timeString = TlUtils::format("[%s %s]", CnTimeX::getNowDate().c_str(), CnTimeX::getNowTime().c_str());

//   std::string title = "<<<< " + stepName + " ";
//   TlUtils::pad(title, (72 - timeString.length()), ' ');
//   title += timeString;

//   // 出力
//   *this << "\n";
//   *this << title << "\n";
//   *this << "CPU TIME  : " << rTime.getCpuTime() << "\n";
//   *this << "ELAPS TIME: " << rTime.getElapsTime() << "\n";
//   *this << line << "\n\n";
// }


// =====================================================================
// TlLogFileX
TlLogFileX::TlLogFileX(const std::string& sFileName)
    : m_sFileName(sFileName)
{
}

TlLogFileX::~TlLogFileX()
{
    if (this->isOpen()) {
        this->m_FileStream.flush();
        this->m_FileStream.close();
    }
}

bool TlLogFileX::isOpen()
{
    return this->m_FileStream.is_open();
}

void TlLogFileX::open()
{
    this->m_FileStream.open(this->m_sFileName.c_str(), std::ios_base::out | std::ios_base::app);
}

void TlLogFileX::flush()
{
    this->m_FileStream.flush();
}


void TlLogFileX::setFilePath(const std::string& path)
{
    if (this->isOpen() == true) {
        this->m_FileStream.flush();
        this->m_FileStream.close();
    }

    this->m_sFileName = path;
}

