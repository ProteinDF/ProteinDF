#ifndef DFDENSITYFITTING_PARALLEL_H
#define DFDENSITYFITTING_PARALLEL_H

#include "DfDensityFittingObject.h"
#include "DfDensityFitting.h"
#include "DfEri_Parallel.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"

// LAPACK版
class DfDensityFitting_Parallel : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri_Parallel> {
public:
    DfDensityFitting_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDensityFitting_Parallel();

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


// ScaLAPACK版
class DfDensityFitting_ScaLAPACK
    : public DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEri_Parallel> {

public:
    DfDensityFitting_ScaLAPACK(TlSerializeData* pPdfParam);
    virtual ~DfDensityFitting_ScaLAPACK();

    void exec();
    
protected:
    virtual void logger(const std::string& str) const;
};


#endif // DFDENSITYFITTING_PARALLEL_H
