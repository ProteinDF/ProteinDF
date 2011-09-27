#ifndef DFTRANSFMATRIX_PARALLEL_H
#define DFTRANSFMATRIX_PARALLEL_H

#include "DfTransFmatrix.h"

class DfTransFmatrix_Parallel : public DfTransFmatrix {
public:
    DfTransFmatrix_Parallel(TlSerializeData* pPdfParam, bool bExecDiis);
    ~DfTransFmatrix_Parallel();

public:
    virtual void DfTrsFmatMain();
    virtual void DfTrsFmatQclo(const std::string& fragname, int norbcut);

protected:
    virtual void logger(const std::string& str) const;

    void DfTrsFmatMain_SCALAPACK();
    void DfTrsFmatQclo_SCALAPACK(const std::string& fragname, int norbcut);
};

#endif // DFTRANSFMATRIX_PARALLEL_H
