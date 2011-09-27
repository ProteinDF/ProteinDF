#ifndef DFINVMATRIX_PARALLEL_H
#define DFINVMATRIX_PARALLEL_H

#include "DfInvMatrix.h"

class DfInvMatrix_Parallel : public DfInvMatrix {
public:
    DfInvMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInvMatrix_Parallel();

    virtual void DfInvMain();

protected:
    virtual void logger(const std::string& str) const;
//     virtual void inverseMatrix(const std::string& sInFile, const std::string& sOutFile);
//     void inverseMatrix_LAPACK(const std::string& sInFile, const std::string& sOutFile);
//     void inverseMatrix_ScaLAPACK(const std::string& sInFile, const std::string& sOutFile);
};

#endif // DFINVMATRIX_PARALLEL_H
