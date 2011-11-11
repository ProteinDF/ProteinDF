#include <iostream>
#include "TlLogging.h"
#include "TlTime.h"
#include "TlUtils.h"

TlLogging* TlLogging::pInstance_ = NULL;
const char* TlLogging::pLevelStr_[] = {
    "DEBUG",
    "INFO",
    "WARN",
    "ERROR",
    "CRITIACAL"
};

TlLogging::TlLogging()
{
    this->filePath_ = "output.log";
    this->procID_ = 0;
    this->setLevel(INFO, WARN);
}

TlLogging::~TlLogging()
{
    this->output_.close();
}

TlLogging& TlLogging::getInstance()
{
    if (TlLogging::pInstance_ == NULL) {
        TlLogging::pInstance_ = new TlLogging();
    }

    return *(TlLogging::pInstance_);
}

void TlLogging::setFilePath(const std::string& filePath)
{
    if (filePath != this->filePath_) {
        this->close();
    }

    this->filePath_ = filePath;
}

void TlLogging::setProcID(const int procID)
{
    this->procID_ = procID;
    this->setLevel();
}

void TlLogging::setLevel(const TlLogging::Level masterLevel,
                         const TlLogging::Level workerLevel)
{
    this->masterLevel_ = masterLevel;
    this->workerLevel_ = workerLevel;
    this->setLevel();
}

void TlLogging::setLevel()
{
    this->level_ = (this->procID_ == 0) ? this->masterLevel_ : this->workerLevel_;
}

void TlLogging::debug(const std::string& msg)
{
    this->output(DEBUG, msg);
}

void TlLogging::info(const std::string& msg)
{
    this->output(INFO, msg);
}

void TlLogging::warn(const std::string& msg)
{
    this->output(WARN, msg);
}

void TlLogging::error(const std::string& msg)
{
    this->output(ERROR, msg);
}

void TlLogging::critical(const std::string& msg)
{
    this->output(CRITICAL, msg);
}

void TlLogging::output(const int level,
                       const std::string& msg)
{
    if (level >= this->level_) {
        const std::string prettyMsg = this->format(level, msg);

        this->open();
        this->output_ << prettyMsg << std::flush;
        
        if (level >= CRITICAL) {
            std::cerr << prettyMsg << std::flush;
        }
    }
}

void TlLogging::open()
{
    if (this->output_.is_open() != true) {
        this->output_.open(this->filePath_.c_str(),
                           std::ios::out | std::ios::app);
    }
}

void TlLogging::close()
{
    if (this->output_.is_open() == true) {
        this->output_.close();
    }    
}


std::string TlLogging::format(const int level, const std::string& input)
{
    // prepare header
    const std::string header = TlUtils::format("[%d:%s:%s] ",
                                               this->procID_,
                                               TlTime::getNow().c_str(),
                                               TlLogging::pLevelStr_[level]);
    const int headerSize = header.size();
    
    // text wrap
    std::string body = TlUtils::textWrap(input, 136 - headerSize);
    TlUtils::rtrim_ws(body);
    
    std::string answer = TlUtils::replace(body, "\n", "\n" + header);

    answer = header + answer + "\n";
    return answer;
}

