#include "TlLog_Parallel.h"

TlLog_Parallel::TlLog_Parallel(TlLog::TlLogLevel level) : TlLog(level)
{
}

TlLog_Parallel::~TlLog_Parallel()
{
}

void TlLog_Parallel::operator<<(std::basic_ostream<char>& (*pf)(std::basic_ostream<char>&))
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster()) {
        TlLog::operator<<(*pf);
    }
}

