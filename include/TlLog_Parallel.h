#ifndef TLLOG_PARALLEL_H
#define TLLOG_PARALLEL_H

#include "TlLog.h"
#include "TlCommunicate.h"

class TlLog_Parallel : public TlLog {
public:
    TlLog_Parallel(TlLog::TlLogLevel level = TlLog::m_defaultLevel);
    virtual ~TlLog_Parallel();

    // operation
public:
    template <typename T>
    TlLog_Parallel& operator<<(const T& t);

    // std::endlが使えるようにするために必要なメンバ関数
    void operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&));
};

////////////////////////////////////////////////////////////////////////
// template
//
template <typename T>
TlLog_Parallel& TlLog_Parallel::operator<<(const T& t)
{
    TlCommunicate* pComm = TlCommunicate::getInstance();

    if (pComm->isMaster()) {
        TlLog::operator<<(t);
    }

    return *this;
}

#endif // TLLOG_PARALLEL_H
