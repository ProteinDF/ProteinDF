#ifndef DFDENSITYFITTINGX_SCALAPACK_H
#define DFDENSITYFITTINGX_SCALAPACK_H

#include "DfDensityFittingObject.h"
#include "DfEriX_Parallel.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"

// 特殊化
template<>
TlDistributeVector DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEriX_Parallel>::calcTAlpha_DIRECT
(const TlDistributeSymmetricMatrix& P)
{
    TlDistributeVector t_alpha(this->m_nNumOfAux);

    DfEriX_Parallel dfEri(this->pPdfParam_);
    dfEri.getJ_D(P, &t_alpha);

    return t_alpha;
}


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
