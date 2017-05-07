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

#include <string>

#ifndef TLPERIODICTABLE_H
#define TLPERIODICTABLE_H

#include <string>

class TlPrdctbl {
public:
    TlPrdctbl();
    ~TlPrdctbl();

public:
    static int getAtomicNumber(const std::string& str);
    static std::string getSymbol(int n);
    static double getVdwRadii(const int atomicNumber);
    static double getBraggSlaterRadii(const int atomicNumber);
    static double getCovalentRadii(const int atomicNumber);
    
private:
    static const char* periodicTable[];
    static const double vdwRadii[];
    static const double BraggSlaterRadii[];
    static const double covalentRadii[];

public:
    static const double UNDEFINED;
};

#endif // TLPERIODICTABLE_H
