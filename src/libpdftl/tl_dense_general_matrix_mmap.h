#ifndef TL_DENSE_GENERAL_MATRIX_MMAP_H
#define TL_DENSE_GENERAL_MATRIX_MMAP_H

#include "tl_dense_matrix_mmap_object.h"

class TlDenseGeneralMatrix_mmap : public TlDenseMatrixMmapObject {
   public:
    TlDenseGeneralMatrix_mmap(const std::string& filePath, const index_type row,
                              const index_type col);
    TlDenseGeneralMatrix_mmap(const std::string& filePath);
    virtual ~TlDenseGeneralMatrix_mmap();

   public:
    void resize(const index_type newRow, const index_type newCol);

   protected:
    // virtual TlDenseGeneralMatrix_mmap* copy(const std::string& path) const;
    virtual size_type getIndex(const index_type row,
                               const index_type col) const;
    virtual TlMatrixObject::size_type getNumOfElements() const;
};

#endif  // TL_DENSE_GENERAL_MATRIX_MMAP_H
