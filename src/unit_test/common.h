#ifndef COMMON_H
#define COMMON_H

#include "TlMatrix.h"

// 以下の要素を設定した行列を返す
// [ 0  1  2 ]
// [ 3  4  5 ]
// [ 6  7  8 ]
TlMatrix getTlMatrixA();

/// ランダムな行列を返す
TlMatrix getTlMatrix(int row, int col);

#endif  // COMMON_H
