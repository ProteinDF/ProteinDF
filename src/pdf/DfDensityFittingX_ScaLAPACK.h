#ifndef DFDENSITYFITTINGX_SCALAPACK_H
#define DFDENSITYFITTINGX_SCALAPACK_H

#include "DfDensityFittingObject.h"
#include "DfEriX_Parallel.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"


// ScaLAPACK版
class DfDensityFittingX_ScaLAPACK
    : public DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEriX_Parallel> {

public:
    DfDensityFittingX_ScaLAPACK(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX_ScaLAPACK();

    void exec();
    
protected:
    virtual void logger(const std::string& str) const;
};


#endif // DFDENSITYFITTINGX_SCALAPACK_H
