// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include "TlLogging.h"

#include <iostream>

#include "TlTime.h"
#include "TlUtils.h"

TlLogging* TlLogging::pInstance_ = NULL;
const char* TlLogging::pLevelStr_[] = {"DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"};

TlLogging::TlLogging() {
    this->filePath_ = "output.log";
    this->procID_ = 0;
    this->setLevel(TL_INFO, TL_WARN);
}

TlLogging::~TlLogging() {
    this->output_.close();
}

TlLogging& TlLogging::getInstance() {
    if (TlLogging::pInstance_ == NULL) {
        TlLogging::pInstance_ = new TlLogging();
    }

    return *(TlLogging::pInstance_);
}

void TlLogging::setFilePath(const std::string& filePath) {
    if (filePath != this->filePath_) {
        this->close();
    }

    this->filePath_ = filePath;
}

void TlLogging::setProcID(const int procID) {
    this->procID_ = procID;
    this->setLevel();
}

void TlLogging::setLevel(const TlLogging::Level masterLevel, const TlLogging::Level workerLevel) {
    this->masterLevel_ = masterLevel;
    this->workerLevel_ = workerLevel;
    this->setLevel();
}

void TlLogging::setLevel() {
    this->level_ = (this->procID_ == 0) ? this->masterLevel_ : this->workerLevel_;
}

void TlLogging::debug(const std::string& msg) {
    this->output(TL_DEBUG, msg);
}

void TlLogging::info(const std::string& msg) {
    this->output(TL_INFO, msg);
}

void TlLogging::warn(const std::string& msg) {
    this->output(TL_WARN, msg);
}

void TlLogging::error(const std::string& msg) {
    this->output(TL_ERROR, msg);
}

void TlLogging::critical(const std::string& msg) {
    this->output(TL_CRITICAL, msg);
}

void TlLogging::output(const int level, const std::string& msg) {
    if (level >= this->level_) {
        const std::string prettyMsg = this->format(level, msg);

        this->open();
        this->output_ << prettyMsg << std::flush;

        if (level >= TL_CRITICAL) {
            std::cerr << prettyMsg << std::flush;
        }
    }
}

void TlLogging::open() {
    if (this->output_.is_open() != true) {
        this->output_.open(this->filePath_.c_str(), std::ios::out | std::ios::app);
    }
}

void TlLogging::close() {
    if (this->output_.is_open() == true) {
        this->output_.close();
    }
}

std::string TlLogging::format(const int level, const std::string& input) {
    // prepare header
    const std::string header =
        TlUtils::format("[%d:%s:%s] ", this->procID_, TlTime::getNow().c_str(), TlLogging::pLevelStr_[level]);
    const int headerSize = header.size();

    // text wrap
    std::string body = TlUtils::textWrap(input, 136 - headerSize);
    TlUtils::rtrim_ws(body);

    std::string answer = TlUtils::replace(body, "\n", "\n" + header);

    answer = header + answer + "\n";
    return answer;
}
