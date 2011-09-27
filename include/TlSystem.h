#ifndef TLSYSTEM_H
#define TLSYSTEM_H

#include <string>

class TlSystem {
public:
    static int getPID();
    static int getPPID();
    static std::string getEnv(const std::string& key);
    
private:
    static int pid_;
    static int ppid_;
};

#endif // TLSYSTEM_H
