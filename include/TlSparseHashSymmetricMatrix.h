#ifndef TLSPARSEHASHSYMMETRICMATRIX_H
#define TLSPARSEHASHSYMMETRICMATRIX_H

#include "TlSparseHashMatrix.h"

class TlSparseHashSymmetricMatrix : public TlSparseHashMatrix { 
public:
  /// 行列オブジェクトを作成する
  ///
  /// @param[in] size 作成する行列の次元数
  explicit TlSparseHashSymmetricMatrix(int size =1);

  /// コピーコンストラクタ
  TlSparseHashSymmetricMatrix(const TlSparseHashSymmetricMatrix& rhs);

  /// デストラクタ
  virtual ~TlSparseHashSymmetricMatrix();

  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  ///
  /// @param[in] row 行数
  /// @param[in] col 列数
  virtual void resize(int row, int col){
    assert(row == col);
    this->resize(row);
  }

  /// 行列のサイズを変更する
  ///
  /// 行列を大きくする場合、追加される要素は0で初期化される。
  /// 行列を小さくする場合、切り詰められる要素は破棄される。
  /// @param[in] size 要素数
  virtual void resize(int size);

  /// 要素を返す(代入可能)
  ///
  /// 内部では、行列要素を(2次元配列ではなく)
  /// 1次元配列として保持しているので、
  /// 他のメンバ関数内でもこのメンバ関数が呼ばれる。
  ///
  ///  @param[in] row 行数
  /// @param[in] col 列数
  /// @return 要素
  virtual double& operator()(int row, int col);

};

inline double& TlSparseHashSymmetricMatrix::operator()(int row, int col) {
  if (row < col){
    std::swap(row, col);
  }

  return TlSparseHashMatrix::operator()(row, col);
}

#endif // TLSPARSEHASHSYMMETRICMATRIX_H
