#ifndef TLTIME_H
#define TLTIME_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <ctime>
#include <string>

#if defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
#include <sys/resource.h>
#include <sys/time.h>
#endif // defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif // HAVE_UNISTD_H

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

#if defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
    static double timeVal2double(struct timeval& tv) {
        return (double)tv.tv_sec + tv.tv_usec * 0.000001;
    }
#endif // defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)

private:
    static const double BILLION;
    bool isRunning_;

    double accumElapseTime_;
    double accumCpuTime_;
    double startElapseTime_;
    double startCpuTime_;
};


extern const TlTime g_GlobalTime;


inline void TlTime::start()
{
    this->isRunning_  = true;

#if defined(HAVE_SYS_RESOURCE_H) && defined(HAVE_SYS_TIME_H)
    {
        struct timeval tv;
        (void)gettimeofday(&tv, NULL);
        this->startElapseTime_ = this->timeVal2double(tv);
        
        struct rusage ru;
        (void)getrusage(RUSAGE_SELF, &ru);
        this->startCpuTime_ = this->timeVal2double(ru.ru_utime) + this->timeVal2double(ru.ru_stime);
    }
#else 
    this->startTime_  = std::difftime(std::time(NULL), 0);
    this->startClock_ = std::difftime(std::clock(), 0);
#endif // HAVE_TIME_H
}


inline void TlTime::stop()
{
    if (this->isRunning() == true) {
        this->accumElapseTime_ += this->getCpuTime();
        this->accumCpuTime_ += this->getElapseTime();
        this->isRunning_ = false;
    }
}


inline void TlTime::reset()
{
    this->isRunning_ = false;
    this->accumElapseTime_ = 0.0;
    this->accumCpuTime_ = 0.0;
    this->startElapseTime_ = 0.0;
    this->startCpuTime_ = 0.0;
}



#endif // TLTIME_H

