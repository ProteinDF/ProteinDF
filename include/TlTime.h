#ifndef TLTIME_H
#define TLTIME_H

#include <ctime>
#include <string>

/**
 *  時間計測クラス
 *
 */
class TlTime {
public:
    TlTime();
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

    /**
     *  基準時刻を現在時刻(呼び出し時刻)に設定する。
     *  @return 基準時刻のタイムスタンプを返す。
     */
    std::string start();

    /**
     *  現在時刻(呼び出し時刻)でタイマーを止める。
     *  @return 基準時刻のタイムスタンプを返す。
     */
    std::string stop();

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

    std::string getReferenceDate() const;
    std::string getReferenceTime() const;

private:
    bool isRunning_;
    std::time_t cumulativeTime_;
    std::clock_t cumulativeClock_;
    
    std::time_t startTime_;
    std::clock_t startClock_;
};

extern const TlTime g_GlobalTime;

#endif // TLTIME_H

