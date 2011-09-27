#ifndef TLFILEDDATAOBJECT_H
#define TLFILEDDATAOBJECT_H

#include <string>
#include <vector>
#include "TlPosition.h"


class TlFieldDataObject {
protected:
    enum FieldType {
        UNIFORM,
        RECTILINEAR,
        IRREGULAR
    };

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


#endif // TLFIELDDATAOBJECT_H
