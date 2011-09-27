#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif // HAVE_UNISTD_H

#include <ctime>

#include "TlTime.h"
#include "TlUtils.h"

const TlTime g_GlobalTime;

TlTime::TlTime()
    : cumulativeTime_(0), cumulativeClock_(0) {
    this->start();
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


std::string TlTime::start()
{
    this->isRunning_  = true;
    this->startTime_  = std::time(NULL);
    this->startClock_ = std::clock();

    return (TlUtils::format("%s %s",
                            this->createDateString(this->startTime_).c_str(),
                            this->createTimeString(this->startTime_).c_str()));
}


std::string TlTime::stop()
{
    std::string answer = "";

    if (this->isRunning() == true) {
        std::time_t endTime = std::time(NULL);
        std::clock_t endClock = std::clock();
        this->cumulativeTime_ += endTime - this->startTime_;
        this->cumulativeClock_ += endClock - this->startClock_;
        this->isRunning_ = false;
        
        answer = TlUtils::format("%s %s",
                                 this->createDateString(this->startTime_).c_str(),
                                 this->createTimeString(this->startTime_).c_str());
    }

    return answer;
}

void TlTime::reset()
{
    this->startTime_ = 0;
    this->startClock_ = 0;
}


// 基準となる時刻からのCPU時間を返す
double TlTime::getCpuTime() const
{
    std::clock_t clocks;
    if (this->isRunning() == true) {
        std::clock_t now = std::clock();
        clocks = now - this->startClock_;
    } else {
        clocks = this->cumulativeClock_;
    }

    return static_cast<double>(clocks / CLOCKS_PER_SEC);
}


// 基準となる時刻からの経過時間を返す
double TlTime::getElapseTime() const
{
    double answer;
    if (this->isRunning() == true) {
        answer = std::difftime(std::time(NULL), this->startTime_);
    } else {
        answer = std::difftime(this->cumulativeTime_, 0);
    }

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


std::string TlTime::getReferenceDate() const
{
    return this->createDateString(this->startTime_);
}


std::string TlTime::getReferenceTime() const
{
    return this->createTimeString(this->startClock_);
}



