#ifndef TL_MATRIX_BLAS_ROW_MAJOR_H
#define TL_MATRIX_BLAS_ROW_MAJOR_H

#include "tl_dense_general_matrix_lapack.h"

class TlMatrixBlasRowMajor : public TlDenseGeneralMatrix_Lapack {
   public:
    explicit TlMatrixBlasRowMajor(const TlMatrixObject::index_type row = 1,
                                  const TlMatrixObject::index_type col = 1);
    virtual ~TlMatrixBlasRowMajor();

   public:
    virtual bool save(const std::string& filePath) const;
};

#endif  // TL_MATRIX_BLAS_ROW_MAJOR_H
