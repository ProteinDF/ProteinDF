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

#ifndef CLASS_TLRESIDUE
#define CLASS_TLRESIDUE

#include <string>

class TlResidue {
   public:
    TlResidue();
    ~TlResidue() {}

    // return number of amino acid date on information table
    int number();

    // return p th amino acid label on information table
    std::string label(int p);

    static bool isResidue(const std::string& str);

   private:
    static const char* residueTable[];
};

#endif
