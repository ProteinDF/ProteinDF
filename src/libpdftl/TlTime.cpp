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

#include "TlTime.h"
#include "TlUtils.h"

const TlTime g_GlobalTime(true);
const double TlTime::BILLION = 1.0E9;

TlTime::TlTime(bool isAutoStart) {
    this->reset();
    if (isAutoStart) {
        this->start();
    }
}

TlTime::~TlTime() {}

std::string TlTime::createDateString(const std::time_t& rTime) {
    char dayBuff[12];
    std::strftime(dayBuff, 12, "%Y/%m/%d", std::localtime(&rTime));

    return std::string(dayBuff);
}

std::string TlTime::createTimeString(const std::time_t& rTime) {
    char timeBuff[9];
    std::strftime(timeBuff, 9, "%H:%M:%S", std::localtime(&rTime));

    return std::string(timeBuff);
}

std::string TlTime::getNow() {
    return (TlTime::getNowDate() + " " + TlTime::getNowTime());
}

// 現在の日付を返す
std::string TlTime::getNowDate() {
    std::time_t now = std::time(NULL);

    return TlTime::createDateString(now);
}

// 現在の時間を返す
std::string TlTime::getNowTime() {
    std::time_t now = std::time(NULL);

    return TlTime::createTimeString(now);
}

// 基準となる時刻からのCPU時間を返す
double TlTime::getCpuTime() const {
    double answer = 0.0;
    double thisCpuTime = 0.0;

    if (this->isRunning()) {
#if defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
        {
            struct rusage ru;
            (void)getrusage(RUSAGE_SELF, &ru);
            thisCpuTime = this->timeVal2double(ru.ru_utime) +
                          this->timeVal2double(ru.ru_stime);
        }
#else
        { thisCpuTime = std::difftime(std::clock(), 0); }
#endif  // defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
        answer = thisCpuTime - this->startCpuTime_;
    } else {
        answer = this->accumCpuTime_;
    }

    return answer;
}

// 基準となる時刻からの経過時間を返す
double TlTime::getElapseTime() const {
    double answer = 0.0;
    double thisElapseTime = 0.0;

    if (this->isRunning()) {
#if defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
        {
            struct timeval tv;
            (void)gettimeofday(&tv, NULL);
            thisElapseTime = TlTime::timeVal2double(tv);
        }
#else
        { thisElapseTime = std::difftime(std::time(NULL), 0); }
#endif  // defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
        answer = thisElapseTime - this->startElapseTime_;
    } else {
        answer = this->accumElapseTime_;
    }

    return answer;
}

void TlTime::sleep(const unsigned long x) {
#ifdef HAVE_USLEEP
    usleep(x * 1000);
#else
    std::clock_t s = std::clock();
    std::clock_t c;
    do {
        c = std::clock();
        if (c == (std::clock_t)-1) {
            return;  // error
        }
    } while ((1000UL * (c - s) / CLOCKS_PER_SEC) <= x);
#endif  // HAVE_USLEEP
}

// std::string TlTime::getReferenceDate() const
// {
//     return this->createDateString(this->startTime_);
// }

// std::string TlTime::getReferenceTime() const
// {
//     return this->createTimeString(this->startClock_);
// }
