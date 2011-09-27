#ifndef DFTRANSATOB_PARALLEL_H
#define DFTRANSATOB_PARALLEL_H

#include "DfTransatob.h"

class DfTransatob_Parallel : public DfTransatob {
public:
    DfTransatob_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfTransatob_Parallel();

    virtual void DfTrsatobMain();
    virtual void DfTrsatobQclo(const std::string& fragname, int norbcut);

protected:
    void DfTrsatobMain_SCALAPACK();
    void DfTrsatobQclo_SCALAPACK(const std::string& fragname, int norbcut);

    virtual void logger(const std::string& str) const;
};

#endif // DFTRANSATOB_PARALLEL_H
