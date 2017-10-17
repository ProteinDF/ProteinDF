#ifndef TlMMAPMATRIX_RSFD_H
#define TlMMAPMATRIX_RSFD_H

#include "TlMmapMatrixObject.h"

class TlMmapMatrix_RSFD : public TlMmapMatrixObject
{
public:
    explicit TlMmapMatrix_RSFD(const std::string& filePath, const index_type row =1, const index_type col =1);
    virtual ~TlMmapMatrix_RSFD();

protected:
    virtual TlMmapMatrix_RSFD* copy(const std::string& path) const;
    virtual size_type getIndex(const index_type row, const index_type col) const;
};

#endif // TlMMAPMATRIX_RSFD_H
