#ifndef DFXMATRIX_PARALLEL_H
#define DFXMATRIX_PARALLEL_H

#include "DfXMatrix.h"

class DfXMatrix_Parallel : public DfXMatrix {
public:
    DfXMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix_Parallel();

public:
    virtual void main();

protected:
    virtual void logger(const std::string& str) const;
    virtual void saveNumOfIndependentBasis();

    void exec_LAPACK();
    void exec_ScaLAPACK();

protected:
    bool fastTrancate_;
};

#endif // DFXMATRIX_PARALLEL_H
