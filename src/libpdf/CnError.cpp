// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include "CnError.h"
#include <cstdlib>
#include <iostream>

#include "TlLogging.h"
#include "TlTime.h"
#include "TlUtils.h"

// global object of CnError Class
CnError CnErr;

CnError::CnError() {}

CnError::~CnError() {}

void CnError::abort(const std::string& msg) {
    TlLogging& log = TlLogging::getInstance();

    log.critical("**** STOP PROGRAM ****");

    if (msg.empty() != true) {
        log.critical(msg);
    }

    std::abort();
}

void CnError::abort(const std::string& ClassName, const std::string& ObjName,
                    const std::string& MemFunc, const std::string& Message) {
    const std::string sMsg =
        TlUtils::format("%-15s::%-15s(%-35s)\n : %s\n", ClassName.c_str(),
                        ObjName.c_str(), MemFunc.c_str(), Message.c_str());
    this->abort(sMsg);
}
