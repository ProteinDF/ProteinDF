#ifndef TLSYSTEM_H
#define TLSYSTEM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <string>

class TlSystem {
public:
    static int getPID();
    static int getPPID();
    static std::string getEnv(const std::string& key);
    
    /// 使用された常駐セットサイズの最大値(MB)を返す
    static double getMaxRSS();

private:
    static int pid_;
    static int ppid_;
};

#endif // TLSYSTEM_H
