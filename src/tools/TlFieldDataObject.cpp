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

#include <cassert>
#include <fstream>
#include <iostream>

#include "TlFieldDataObject.h"
#include "TlUtils.h"

const char* TlFieldDataObject::DataTypeStrings[] = {"byte", "short", "integer",
                                                    "float", "double"};

const char* TlFieldDataObject::FieldTypeStrings[] = {"uniform", "rectilinear",
                                                     "irregular"};

TlFieldDataObject::TlFieldDataObject() : comment_(""), label_("") {}

TlFieldDataObject::TlFieldDataObject(const TlFieldDataObject& rhs)
    : comment_(rhs.comment_), label_(rhs.label_) {}

TlFieldDataObject::~TlFieldDataObject() {}

void TlFieldDataObject::setComment(const std::string& comment) {
    this->comment_ = comment;
}

std::string TlFieldDataObject::getComment() const { return this->comment_; }

void TlFieldDataObject::setLabel(const std::string& label) {
    this->label_ = label;
}

std::string TlFieldDataObject::getLabel() const { return this->label_; }
