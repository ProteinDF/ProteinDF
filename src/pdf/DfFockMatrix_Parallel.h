#ifndef DFFOCKMATRIX_PARALLEL_H
#define DFFOCKMATRIX_PARALLEL_H

#include "DfFockMatrix.h"

class DfFockMatrix_Parallel : public DfFockMatrix {
public:
    DfFockMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfFockMatrix_Parallel();

protected:
    virtual void logger(const std::string& str) const;
    
protected:
    virtual void mainDIRECT_RKS();
    // void mainDIRECT_RKS_LAPACK();
    // void mainDIRECT_RKS_ScaLAPACK();

    virtual void mainDIRECT_UKS();
    virtual void mainDIRECT_ROKS();

    virtual void setXC_RI(RUN_TYPE nRunType, TlSymmetricMatrix& F);
    virtual void setXC_DIRECT(RUN_TYPE nRunType, TlSymmetricMatrix& F);

    virtual void setCoulomb(METHOD_TYPE nRunType, TlSymmetricMatrix& F);
    void setCoulomb(const METHOD_TYPE nMethodType, TlDistributeSymmetricMatrix& F);

    virtual TlSymmetricMatrix getFpqMatrix(RUN_TYPE nRunType, int nIteration);

    virtual TlVector getRho(RUN_TYPE nRunType, int nIteration);
    virtual TlVector getMyu(RUN_TYPE nRunType, int nIteration);
    //virtual void saveFpqMatrix(RUN_TYPE nRunType, const TlSymmetricMatrix& F);

};

#endif // DFFOCKMATRIX_PARALLEL_H

