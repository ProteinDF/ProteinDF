#include <iostream>
#include <cstdlib>
#include "CnError.h"

#include "TlLogX.h"
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

void CnError::abort(const std::string& sMsg)
{
    TlLogX& Log = TlLogX::getInstance();

    Log.fatal("**** ERROR WAS DETECTED ****\n");
    Log.fatal("**** TERMINATE PROGRAM  ****\n");

    if (sMsg != "") {
        Log.fatal(sMsg + "\n");
    }

    Log.fatal(TlTime::getNow() + "\n");
    Log.fatal("**** DETAILS WERE DESCRIBED IN 'fl_Out_Std' ****\n");

    std::abort();
}

void CnError::abort(const std::string& ClassName, const std::string& ObjName, const std::string& MemFunc, const std::string& Message)
{
    const std::string sMsg = TlUtils::format("%-15s::%-15s(%-35s)\n : %s\n",
                                             ClassName.c_str(), ObjName.c_str(), MemFunc.c_str(), Message.c_str());
    this->abort(sMsg);
}

