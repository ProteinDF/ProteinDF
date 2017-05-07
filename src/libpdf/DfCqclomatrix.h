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

#ifndef DFCQCLOMATRIX_H
#define DFCQCLOMATRIX_H

#include "DfObject.h"

//  function : enaluate X-matrix by diagonalizeing orbital based S-matrix
//
//  input : S (hermite,positive define matrix)
//
//  output : Cqclo (Unitary ,column/row oriented,unpacked)

class DfCqclomatrix : public DfObject {
public:
    // flGbi は書き換えられる
    DfCqclomatrix(TlSerializeData* pPdfParam);
    ~DfCqclomatrix();

    void main();

private:
    void main(std::string type);

private:
    int number_fragment;           // Number of fragments
};

#endif // DFCQCLOMATRIX_H

