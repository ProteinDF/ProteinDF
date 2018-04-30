#include "common.h"

TlMatrix getTlMatrixA() {
  TlMatrix a(3, 3);
  a(0, 0) = 0.0;
  a(0, 1) = 1.0;
  a(0, 2) = 2.0;
  a(1, 0) = 3.0;
  a(1, 1) = 4.0;
  a(1, 2) = 5.0;
  a(2, 0) = 6.0;
  a(2, 1) = 7.0;
  a(2, 2) = 8.0;

  return a;
}

TlMatrix getTlMatrix(int row, int col) {
  TlMatrix m(row, col);

  std::srand((unsigned int)time(NULL));
  for (int r = 0; r < row; ++r) {
    for (int c = 0; c < col; ++c) {
      m.set(r, c, double(double(std::rand()) / double(RAND_MAX)));
    }
  }

  return m;
}
