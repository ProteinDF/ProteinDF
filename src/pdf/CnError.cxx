#include <iostream>
#include <cstdlib>
#include "CnError.h"

#include "TlLogging.h"
#include "TlUtils.h"
#include "TlTime.h"

// global object of CnError Class
CnError CnErr;

CnError::CnError()
{
}

CnError::~CnError()
{
}

void CnError::abort(const std::string& msg)
{
    TlLogging& log = TlLogging::getInstance();

    log.critical("**** STOP PROGRAM ****");

    if (msg.empty() != true) {
        log.critical(msg);
    }

    std::abort();
}

void CnError::abort(const std::string& ClassName, const std::string& ObjName, const std::string& MemFunc, const std::string& Message)
{
    const std::string sMsg = TlUtils::format("%-15s::%-15s(%-35s)\n : %s\n",
                                             ClassName.c_str(), ObjName.c_str(), MemFunc.c_str(), Message.c_str());
    this->abort(sMsg);
}

