#ifndef TL_DENSE_GENERAL_MATRIX_MMAP_H
#define TL_DENSE_GENERAL_MATRIX_MMAP_H

#include "tl_dense_general_matrix_object.h"
#include "tl_dense_matrix_mmap_object.h"

class TlDenseGeneralMatrix_mmap : public TlDenseMatrixMmapObject {
   public:
    TlDenseGeneralMatrix_mmap(const std::string& filePath, const index_type row,
                              const index_type col);
    TlDenseGeneralMatrix_mmap(const std::string& filePath);
    virtual ~TlDenseGeneralMatrix_mmap();

   public:
    void resize(const index_type newRow, const index_type newCol);

   public:
    /// ブロック行列を返す
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] rowDistance 取得する行数
    /// @param[in] colDistance 取得する列数
    /// @return row_distance × col_distance 次元のブロック行列
    void block(const TlMatrixObject::index_type row,
               const TlMatrixObject::index_type col,
               const TlMatrixObject::index_type rowDistance,
               const TlMatrixObject::index_type colDistance,
               TlDenseGeneralMatrixObject* pOut) const;

    /// 行列要素を指定された位置に上書きする
    ///
    /// @param[in] row 始点となる行
    /// @param[in] col 始点となる列
    /// @param[in] ref 参照行列
    void block(const TlMatrixObject::index_type row,
               const TlMatrixObject::index_type col,
               const TlDenseGeneralMatrixObject& ref);

   protected:
    // virtual TlDenseGeneralMatrix_mmap* copy(const std::string& path) const;
    virtual size_type getIndex(const index_type row,
                               const index_type col) const;
    virtual TlMatrixObject::size_type getNumOfElements() const;
};

#endif  // TL_DENSE_GENERAL_MATRIX_MMAP_H
