#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif // HAVE_UNISTD_H

#include <ctime>

#include "TlTime.h"
#include "TlUtils.h"

const TlTime g_GlobalTime(true);
const double TlTime::BILLION = 1.0E9;

TlTime::TlTime(bool isAutoStart)
{
    this->reset();
    if (isAutoStart) {
        this->start();
    }
}


TlTime::~TlTime()
{
}


std::string TlTime::createDateString(const std::time_t& rTime)
{
    char dayBuff[12];
    std::strftime(dayBuff, 12, "%Y/%m/%d", std::localtime(&rTime));

    return std::string(dayBuff);
}


std::string TlTime::createTimeString(const std::time_t& rTime)
{
    char timeBuff[9];
    std::strftime(timeBuff, 9, "%H:%M:%S", std::localtime(&rTime));

    return std::string(timeBuff);
}


std::string TlTime::getNow()
{
    return (TlTime::getNowDate() + " " + TlTime::getNowTime());
}


// 現在の日付を返す
std::string TlTime::getNowDate()
{
    std::time_t now = std::time(NULL);

    return TlTime::createDateString(now);
}


// 現在の時間を返す
std::string TlTime::getNowTime()
{
    std::time_t now = std::time(NULL);

    return TlTime::createTimeString(now);
}


// 基準となる時刻からのCPU時間を返す
double TlTime::getCpuTime() const
{
    double answer = 0.0;

#ifdef HAVE_TIME_H
    if (this->isRunning() == true) {
        struct timespec stop;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
        answer = (stop.tv_sec - this->startCpuTime_.tv_sec) 
            + (double)(stop.tv_nsec - this->startCpuTime_.tv_nsec) / TlTime::BILLION;
    } else {
        answer = this->accumCpuTime_.tv_sec 
            + (double)(this->accumCpuTime_.tv_nsec) / TlTime::BILLION;
    }
#else
    std::clock_t clocks;
    if (this->isRunning() == true) {
        std::clock_t now = std::clock();
        clocks = now - this->startClock_;
    } else {
        clocks = this->cumulativeClock_;
    }
    answer = static_cast<double>(clocks / CLOCKS_PER_SEC);
#endif // HAVE_TIME_H

    return answer;
}


// 基準となる時刻からの経過時間を返す
double TlTime::getElapseTime() const
{
    double answer = 0.0;

#ifdef HAVE_TIME_H
    if (this->isRunning() == true) {
        struct timespec stop;
        clock_gettime(CLOCK_MONOTONIC, &stop);
        answer = (stop.tv_sec - this->startElapseTime_.tv_sec) 
            + (double)(stop.tv_nsec - this->startElapseTime_.tv_nsec) / TlTime::BILLION;
    } else {
        answer = this->accumElapseTime_.tv_sec 
            + (double)(this->accumElapseTime_.tv_nsec) / TlTime::BILLION;
    }
#else
    if (this->isRunning() == true) {
        answer = std::difftime(std::time(NULL), this->startTime_);
    } else {
        answer = std::difftime(this->cumulativeTime_, 0);
    }
#endif // HAVE_TIME_H

    return answer;
}


void TlTime::sleep(const unsigned long x)
{
#ifdef HAVE_USLEEP
    usleep(x * 1000);
#else
    std::clock_t s = std::clock();
    std::clock_t c;
    do {
        c = std::clock();
        if (c == (std::clock_t)-1) {
            return; // error
        }
    } while ((1000UL * (c - s) / CLOCKS_PER_SEC) <= x);
#endif // HAVE_USLEEP
}


// std::string TlTime::getReferenceDate() const
// {
//     return this->createDateString(this->startTime_);
// }


// std::string TlTime::getReferenceTime() const
// {
//     return this->createTimeString(this->startClock_);
// }



