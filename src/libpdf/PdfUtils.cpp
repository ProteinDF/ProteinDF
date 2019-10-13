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

#include "PdfUtils.h"
#include "TlUtils.h"

bool PdfUtils::isComment(const std::string& str) {
    std::string check = str;
    const int length = str.size();

    TlUtils::trim_ws(check);
    if ((length >= 1) && (check.substr(0, 1) == "#")) {
        return true;
    }
    if (length >= 2) {
        const std::string head = check.substr(0, 2);
        if ((head == "//") || (head == "--")) {
            return true;
        }
    }

    return false;
}
