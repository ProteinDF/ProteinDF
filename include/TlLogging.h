#ifndef TLLOGGING_H
#define TLLOGGING_H

#include <string>
#include <fstream>

/// ログを扱うクラス
/// 
/// MPI並列プログラム用に出力レベルを2つ(master用, worker用)用意しています。
/// 出力レベル以上の出力が無ければ無視されます。
class TlLogging {
public:
    enum Level {
        DEBUG,
        INFO,
        WARN,
        ERROR,
        CRITICAL
    };

public:
    static TlLogging& getInstance();

public:
    void setFilePath(const std::string& filePath);
    void setProcID(const int procID);

    void setLevel(const TlLogging::Level masterLevel,
                  const TlLogging::Level workerLevel = WARN);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg);
    void error(const std::string& msg);
    void critical(const std::string& msg);
    
private:
    TlLogging();
    ~TlLogging();

    void setLevel();
    
    void open();
    void close();

    void output(const int level,
                const std::string& msg);
    std::string format(const int level,
                       const std::string& input);
    
private:
    static TlLogging* pInstance_;
    static const char* pLevelStr_[];
    
    std::string filePath_;
    
    /// 出力用ファイルストリーム
    std::fstream output_;

    /// 出力レベル
    Level masterLevel_;
    Level workerLevel_;
    Level level_;
    
    /// プロセスIDを保持
    int procID_;
};

#endif // TLLOGGING_H
