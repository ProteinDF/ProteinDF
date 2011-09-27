#ifndef DFCLEANUP_H
#define DFCLEANUP_H

#include "DfObject.h"

class DfCleanup : public DfObject {
public:
    DfCleanup(TlSerializeData* pPdfParam);
    virtual ~DfCleanup();

public:
    virtual void cleanup();

protected:
    void cleanup(RUN_TYPE runType, int iteration);
    void cleanupFxc(RUN_TYPE runType, int iteration);
    void cleanupHFx(RUN_TYPE runType, int iteration);
    void cleanupFpq(RUN_TYPE runType, int iteration);
    void cleanupFprime(RUN_TYPE runType, int iteration);
    void cleanupCprime(RUN_TYPE runType, int iteration);
    void cleanupC(RUN_TYPE runType, int iteration);
    void cleanupP(RUN_TYPE runType, int iteration);
};

#endif // DFCLEANUP_H
