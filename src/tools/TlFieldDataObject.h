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

#ifndef TLFILEDDATAOBJECT_H
#define TLFILEDDATAOBJECT_H

#include <string>
#include <vector>
#include "TlPosition.h"

class TlFieldDataObject {
   protected:
    enum FieldType { UNIFORM, RECTILINEAR, IRREGULAR };

   public:
    TlFieldDataObject();
    TlFieldDataObject(const TlFieldDataObject& rhs);
    virtual ~TlFieldDataObject();

   public:
    void setComment(const std::string& comment);
    std::string getComment() const;

    void setLabel(const std::string& label);
    std::string getLabel() const;

   protected:
    static const char* DataTypeStrings[];
    static const char* FieldTypeStrings[];

    /// コメント
    std::string comment_;

    /// ラベル
    std::string label_;
};

#endif  // TLFIELDDATAOBJECT_H
