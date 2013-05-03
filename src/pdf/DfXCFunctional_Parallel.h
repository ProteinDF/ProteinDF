#ifndef DFXCFUNCTIONAL_PARALLEL_H
#define DFXCFUNCTIONAL_PARALLEL_H

#include "DfXCFunctional.h"
#include "DfCalcGridX_Parallel.h"
#include "CnError.h"

// 交換・相関項の計算を行う
class DfXCFunctional_Parallel : public DfXCFunctional {
public:
    DfXCFunctional_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXCFunctional_Parallel();

public:
    virtual void buildXcMatrix();

    virtual double getGrimmeDispersionEnergy();
    
protected:
    virtual void logger(const std::string& str) const;

    void buildXC_LAPACK();
    void buildXC_ScaLAPACK();

protected:
    virtual DfEriX* getDfEriXObject();
};


#endif // DFXCFUNCTIONAL_PARALLEL_H
