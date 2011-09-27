#ifndef TLROTATELCAO_H
#define TLROTATELCAO_H

#include "TlOrbitalInfo.h"
#include "TlMatrix.h"

class TlRotateLCAO
{
public:
    TlRotateLCAO(const TlOrbitalInfo& orbInfo);
    ~TlRotateLCAO();

public:
    TlMatrix exec(const TlMatrix& lcao, const TlMatrix& rot);

protected:
    void rotateLCAO_typeS(const TlMatrix& lcao, const TlMatrix& rot,
                          const int row, const int maxCol,
                          TlMatrix* ioMatrix);
    void rotateLCAO_typeP(const TlMatrix& lcao, const TlMatrix& rot,
                          const int row, const int maxCol,
                          TlMatrix* ioMatrix);
    void rotateLCAO_typeD(const TlMatrix& lcao, const TlMatrix& rot,
                          const int row, const int maxCol,
                          TlMatrix* ioMatrix);
    
protected:
    TlOrbitalInfo orbInfo_;
    
};

#endif // TLROTATELCAO_H
