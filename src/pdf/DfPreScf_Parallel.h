#ifndef DFPRESCF_PARALLEL_H
#define DFPRESCF_PARALLEL_H

#include "DfPreScf.h"

/// DfScfクラスの前処理を行うクラス
class DfPreScf_Parallel : public DfPreScf {
public:
    DfPreScf_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfPreScf_Parallel();
    
protected:
    virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);
    virtual void createOccupation(const RUN_TYPE runType);

protected:
    void createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType);
    void createInitialGuessUsingLCAO_onDisk(const RUN_TYPE runType);

    TlDistributeMatrix getLCAO_onScaLAPACK(const RUN_TYPE runType);
    
    virtual void logger(const std::string& str) const;
};

#endif // DFPRESCF_PARALLEL_H
