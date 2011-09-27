#ifndef DFRLMO_H
#define DFRLMO_H

#include <vector>
#include "DfObject.h"
#include "TlMatrix.h"

class DfRLMO : public DfObject {
public:
    DfRLMO(TlSerializeData* pPdfParam);
    virtual ~DfRLMO();
    
public:
    void exec(const std::vector<index_type>& startBlockAOs);

protected:
    TlSymmetricMatrix getX();

    TlMatrix getT(const TlSymmetricMatrix& D,
                  const std::vector<index_type>& startBlockAOs);
};

#endif // DFRLMO_H
