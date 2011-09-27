#include <cassert>
#include <iostream>
#include <fstream>

#include "TlFieldDataObject.h"
#include "TlUtils.h"


const char* TlFieldDataObject::DataTypeStrings[] = {
    "byte", "short", "integer", "float", "double"
};

const char* TlFieldDataObject::FieldTypeStrings[] = {
    "uniform", "rectilinear", "irregular"
};


TlFieldDataObject::TlFieldDataObject()
    : comment_(""), label_("") 
{
}


TlFieldDataObject::TlFieldDataObject(const TlFieldDataObject& rhs)
    : comment_(rhs.comment_), label_(rhs.label_)
{
}


TlFieldDataObject::~TlFieldDataObject()
{
}

    
void TlFieldDataObject::setComment(const std::string& comment)
{
    this->comment_ = comment;
}


std::string TlFieldDataObject::getComment() const
{
    return this->comment_;
}


void TlFieldDataObject::setLabel(const std::string& label)
{
    this->label_ = label;
}


std::string TlFieldDataObject::getLabel() const
{
    return this->label_;
}


