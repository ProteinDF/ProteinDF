#ifndef TlMMAPMATRIX_CSFD_H
#define TlMMAPMATRIX_CSFD_H

#include "TlMmapMatrixObject.h"

class TlMmapMatrix_CSFD : public TlMmapMatrixObject
{
public:
    explicit TlMmapMatrix_CSFD(const std::string& filePath, const index_type row =1, const index_type col =1);
    virtual ~TlMmapMatrix_CSFD();

protected:
    virtual TlMmapMatrix_CSFD* copy(const std::string& path) const;
    virtual size_type getIndex(const index_type row, const index_type col) const;
};

#endif // TlMMAPMATRIX_CSFD_H
