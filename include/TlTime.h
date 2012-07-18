#ifndef TLTIME_H
#define TLTIME_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <ctime>
#include <string>

/**
 *  時間計測クラス
 *
 */
class TlTime {
public:
    explicit TlTime(bool isAutoStart = false);
    ~TlTime();

public:
    /**
     * 現在の日時と時刻を返す。
     *
     * @return 現在の日時と時刻
     */
    static std::string getNow();

    /**
     * 現在の日時を返す。
     *
     * @return 現在の日時
     */
    static std::string getNowDate();

    /**
     * 現在の時刻を返す。
     *
     * @return 現在の時刻
     */
    static std::string getNowTime();

    /// 基準時刻を現在時刻(呼び出し時刻)に設定する。
    void start();

    /// 現在時刻(呼び出し時刻)でタイマーを止める。
    void stop();

    void reset();
    
    bool isRunning() const {
        return this->isRunning_;
    }   
    
    /**
     *  基準時刻からのCPU時間を返す。
     */
    double getCpuTime() const;

    /**
     *  基準時刻からの経過時間を返す。
     */
    double getElapseTime() const;

public:
    /// xミリ秒待つ
    static void sleep(unsigned long x);

private:
    /**
     *  指定日付の文字列を返す。
     *
     *  @param [in] rTime 指定する時刻
     *  @return 日付の文字列
     */
    static std::string createDateString(const std::time_t& rTime);

    /**
     *  指定時刻の文字列を返す。
     *
     *  @param [in] rTime 指定する時刻
     *  @return 時刻の文字列
     */
    static std::string createTimeString(const std::time_t& rTime);

    // std::string getReferenceDate() const;
    // std::string getReferenceTime() const;

private:
    static const double BILLION;
    bool isRunning_;

#ifdef HAVE_TIME_H
    struct timespec accumElapseTime_;
    struct timespec accumCpuTime_;
    struct timespec startElapseTime_;
    struct timespec startCpuTime_;
#else
    std::time_t cumulativeTime_;
    std::clock_t cumulativeClock_;
    std::time_t startTime_;
    std::clock_t startClock_;
#endif // HAVE_TIME_H
};

extern const TlTime g_GlobalTime;

inline void TlTime::start()
{
#pragma omp critical(TlTime)
    {
        this->isRunning_  = true;
#ifdef HAVE_TIME_H
        clock_gettime(CLOCK_MONOTONIC, &(this->startElapseTime_));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(this->startCpuTime_));
#else
        this->startTime_  = std::time(NULL);
        this->startClock_ = std::clock();
#endif // HAVE_TIME_H
    }
}


inline void TlTime::stop()
{
#pragma omp critical(TlTime)
    {
        if (this->isRunning() == true) {
            this->isRunning_ = false;
#ifdef HAVE_TIME_H
            struct timespec stopElapseTime;
            struct timespec stopCpuTime;
            clock_gettime(CLOCK_MONOTONIC, &stopElapseTime);
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stopCpuTime);
            this->accumElapseTime_.tv_sec += stopElapseTime.tv_sec - this->startElapseTime_.tv_sec;
            this->accumElapseTime_.tv_nsec += stopElapseTime.tv_nsec - this->startElapseTime_.tv_nsec;
            this->accumCpuTime_.tv_sec += stopCpuTime.tv_sec - this->startCpuTime_.tv_sec;
            this->accumCpuTime_.tv_nsec += stopCpuTime.tv_nsec - this->startCpuTime_.tv_nsec;
#else
            std::time_t endTime = std::time(NULL);
            std::clock_t endClock = std::clock();
            this->cumulativeTime_ += endTime - this->startTime_;
            this->cumulativeClock_ += endClock - this->startClock_;
#endif // HAVE_TIME_H
        }
    }
}


inline void TlTime::reset()
{
#pragma omp critical(TlTime)
    {
#ifdef HAVE_TIME_H
        clock_gettime(CLOCK_MONOTONIC, &(this->startElapseTime_));
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(this->startCpuTime_));
        this->accumElapseTime_.tv_sec = 0;
        this->accumElapseTime_.tv_nsec = 0;
        this->accumCpuTime_.tv_sec = 0;
        this->accumCpuTime_.tv_nsec = 0;
#else
        this->cumulativeTime_ = 0;
        this->cumulativeClock_ = 0;
#endif // HAVE_TIME_H
    }
}



#endif // TLTIME_H

