#ifndef CNFILE_H
#define CNFILE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>
#include "DfObject.h"

class TlMatrixObject;

class CnFile {
   public:
    CnFile();
    virtual ~CnFile();

   public:
    // virtual void getPMatrix(const DfObject::RUN_TYPE, int itr,
    // TlMatrixObject* pMatrix);
    virtual void saveMatrix(const std::string& path,
                            const TlMatrixObject& matrix);
};

#endif  // CNFILE_H
