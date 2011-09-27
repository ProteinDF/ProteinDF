#ifndef TLESPFIELD_H
#define TLESPFIELD_H

#include <vector>
#include "TlPosition.h"
#include "TlSymmetricMatrix.h"
#include "TlSerializeData.h"

class TlFieldData_Uniform;

class TlEspField {
public:
    TlEspField(const TlSerializeData& param);
    ~TlEspField();

    std::vector<double> makeEspFld(const TlSymmetricMatrix& P,
                                   const std::vector<TlPosition>& grids);

protected:
    TlSerializeData param_;
};

#endif // TLESPFIELD_H




