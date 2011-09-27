#include "TlSystem.h"
#include <sys/types.h>
#include <unistd.h>
#include <cstdlib>

int TlSystem::pid_ = -1;
int TlSystem::ppid_ = -1;

int TlSystem::getPID()
{
    if (TlSystem::pid_ == -1) {
        TlSystem::pid_ = getpid();
    }
    return TlSystem::pid_;
}


int TlSystem::getPPID()
{
    if (TlSystem::ppid_ == -1) {
        TlSystem::ppid_ = getppid();
    }
    return TlSystem::ppid_;
}


std::string TlSystem::getEnv(const std::string& key)
{
    std::string ans(std::getenv(key.c_str()));

    return ans;
}
