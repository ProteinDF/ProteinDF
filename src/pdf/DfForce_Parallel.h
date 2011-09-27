#ifndef DFFORCE_PARALLEL_H
#define DFFORCE_PARALLEL_H

#include "DfForce.h"

class DfForce_Parallel : public DfForce {
public:
    DfForce_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfForce_Parallel();

protected:
    virtual void logger(const std::string& str) const;
    
protected:
    virtual void calcForceFromCoulomb_RIJ(const RUN_TYPE runType);
    void calcForceFromCoulomb_RIJ_DC(const RUN_TYPE runType);

    virtual void calcForceFromK(const RUN_TYPE runType);
    void calcForceFromK_DC(const RUN_TYPE runType);
};


#endif // DFFORCE_PARALLEL_H
