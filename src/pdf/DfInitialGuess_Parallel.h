#ifndef DFINITIALGUESS_PARALLEL_H
#define DFINITIALGUESS_PARALLEL_H

#include "DfInitialGuess.h"
#include "TlDistributeMatrix.h"

class DfInitialGuess_Parallel : public DfInitialGuess {
public:
    DfInitialGuess_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess_Parallel();
    
protected:
    virtual void createRho();
    
    virtual TlVector createOccupation(RUN_TYPE runType);

    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();
    virtual void createInitialGuessUsingHarris();

    virtual void createInitialGuessUsingLCAO(RUN_TYPE runType);

    /// 占有軌道情報を取得する
    virtual TlVector getOccupation(const RUN_TYPE runType);

    /// 占有軌道情報を保存する
    virtual void saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation);

protected:
    void createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType);
    void createInitialGuessUsingLCAO_onLAPACK(const RUN_TYPE runType);

    TlDistributeMatrix getLCAO_onScaLAPACK(const RUN_TYPE runType);

    virtual DfDmatrix* getDfDmatrixObject(TlSerializeData* param);
};

#endif // DFINITIALGUESS_PARALLEL_H
