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

#include <iostream>
#include <string>
#include "DfLocalize.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

int main(int argc, char *argv[])
{
    std::string pdfparam_path = "pdfparam.mpac"; 

    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfparam_path);
        param = mpac.getSerializeData();
    }
    
    DfLocalize lo(&param);

    lo.localize();

    {
        TlMsgPack mpac(param);
        mpac.save(pdfparam_path);
    }

    return 0;
}




