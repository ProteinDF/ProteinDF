#ifndef DFDMATRIX_PARALLEL_H
#define DFDMATRIX_PARALLEL_H

#include "DfDmatrix.h"
#include "TlVector.h"

class DfDmatrix_Parallel : public DfDmatrix {
public:
    DfDmatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDmatrix_Parallel();

protected:
    virtual void main(DfObject::RUN_TYPE runType);
    void main_SCALAPACK(DfObject::RUN_TYPE runType);

    virtual void checkOccupation(const TlVector& prevOcc, const TlVector& currOcc);
    virtual void printOccupation(const TlVector& occ);

    virtual TlVector getOccupation(DfObject::RUN_TYPE runType);
};

#endif // DFDMATRIX_PARALLEL_H
