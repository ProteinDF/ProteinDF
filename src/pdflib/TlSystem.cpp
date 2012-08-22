#include "TlSystem.h"
#include <sys/types.h>
#include <cstdlib>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif // HAVE_UNISTD_H

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif // HAVE_SYS_RESOURCE_H

int TlSystem::pid_ = -1;
int TlSystem::ppid_ = -1;

int TlSystem::getPID()
{
#ifdef HAVE_UNISTD_H
    if (TlSystem::pid_ == -1) {
        TlSystem::pid_ = getpid();
    }
#endif // HAVE_UNISTD_H

    return TlSystem::pid_;
}


int TlSystem::getPPID()
{
#ifdef HAVE_UNISTD_H
    if (TlSystem::ppid_ == -1) {
        TlSystem::ppid_ = getppid();
    }
#endif // HAVE_UNISTD_H
    return TlSystem::ppid_;
}


std::string TlSystem::getEnv(const std::string& key)
{
    std::string ans(std::getenv(key.c_str()));

    return ans;
}


double TlSystem::getMaxRSS()
{
    double answer = 0.0;

#ifdef HAVE_SYS_RESOURCE_H
    struct rusage resourceUsage;
    getrusage(RUSAGE_SELF, &resourceUsage);
    answer = resourceUsage.ru_maxrss / 1024.0;
#endif // HAVE_SYS_RESOURCE_H

    return answer;
}


