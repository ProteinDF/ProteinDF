#ifndef DFINPUTDATA_PARALLEL_H
#define DFINPUTDATA_PARALLEL_H

#include "DfInputdata.h"

class DfInputdata_Parallel : public DfInputdata {
public:
    DfInputdata_Parallel();
    virtual ~DfInputdata_Parallel();
    
public:
    virtual TlSerializeData main();
};

#endif // DFINPUTDATA_PARALLEL_H

