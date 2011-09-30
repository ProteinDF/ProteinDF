#ifndef DFDENSITYFITTINGX_PARALLEL_H
#define DFDENSITYFITTINGX_PARALLEL_H

#include "DfDensityFittingObject.h"
//#include "DfDensityFitting.h"
#include "DfEriX_Parallel.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"

// LAPACK版
class DfDensityFittingX_Parallel : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX_Parallel> {
public:
    DfDensityFittingX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX_Parallel();

    void exec();
    
protected:
    virtual TlVector getNalpha();
    virtual TlSymmetricMatrix getSabinv();
    virtual TlVector calcTAlpha_DIRECT(const TlSymmetricMatrix& P);
    virtual TlVector getTalpha(RUN_TYPE runType, const int iteration);
    virtual void getTalpha_ROKS(TlVector* pT_alphaA, TlVector* pT_alphaB);
    virtual TlSymmetricMatrix getDiffDensityMatrix(RUN_TYPE runType);
    virtual TlSymmetricMatrix getP1pq(const int nIteration);
    virtual TlSymmetricMatrix getP2pq(const int nIteration);
    virtual double getLamda(const TlVector& SabinvN, const TlVector& t_alpha,
                            const TlVector& N_alpha, const double dNumOfElec);
    virtual void saveRho(const TlVector& rRho, RUN_TYPE runType);

    virtual void logger(const std::string& str) const;
};


#endif // DFDENSITYFITTINGX_PARALLEL_H
