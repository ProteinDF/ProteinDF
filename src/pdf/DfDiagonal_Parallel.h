#ifndef DFDIAGONAL_PARALLEL_H
#define DFDIAGONAL_PARALLEL_H

#include "DfDiagonal.h"

class DfDiagonal_Parallel : public DfDiagonal {
public:
    DfDiagonal_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDiagonal_Parallel();

    virtual void DfDiagMain();
    virtual void DfDiagQclo(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut);

protected:
    void DfDiagMain_SCALAPACK();
    void DfDiagQclo_SCALAPACK(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut);
};

#endif // DFDIAGONAL_PARALLEL_H
