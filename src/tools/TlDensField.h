#ifndef TLDENSFIELD_H
#define TLDENSFIELD_H

#include "TlPosition.h"
#include "TlSymmetricMatrix.h"
#include "TlSerializeData.h"

class TlDensField {
public:
    TlDensField(const TlSerializeData& param);
    ~TlDensField();

    std::vector<double> makeDensFld(const TlSymmetricMatrix& P,
                                    const std::vector<TlPosition>& grids);

protected:
    TlSerializeData param_;
};

#endif // TLDENSFIELD_H




