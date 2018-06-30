#ifndef MATRIX_COMMON_H
#define MATRIX_COMMON_H

template <typename T>
T getTlMatrix(int row, int col) {
  T m(row, col);

  std::srand((unsigned int)time(NULL));
  for (int r = 0; r < row; ++r) {
    for (int c = 0; c < col; ++c) {
      m.set(r, c, double(double(std::rand()) / double(RAND_MAX)));
    }
  }

  return m;
}

//          1 th       2 th       3 th
// -----------------------------------------
// 1     0.000000   1.000000   2.000000
// 2     3.000000   4.000000   5.000000
// 3     6.000000   7.000000   8.000000
template <typename T>
T getMatrixA() {
  T a(3, 3);
  a.set(0, 0, 0.0);
  a.set(0, 1, 1.0);
  a.set(0, 2, 2.0);
  a.set(1, 0, 3.0);
  a.set(1, 1, 4.0);
  a.set(1, 2, 5.0);
  a.set(2, 0, 6.0);
  a.set(2, 1, 7.0);
  a.set(2, 2, 8.0);

  return a;
};

template <typename T>
T getMatrixB() {
  T b(3, 3);
  b.set(0, 0, 0.0);
  b.set(1, 0, 1.0);
  b.set(2, 0, 2.0);
  b.set(0, 1, 3.0);
  b.set(1, 1, 4.0);
  b.set(2, 1, 5.0);
  b.set(0, 2, 6.0);
  b.set(1, 2, 7.0);
  b.set(2, 2, 8.0);

  return b;
}

// [ 1   2  3 ]
// [ 2  -1  1 ]
// [ 4   3  2 ]
template <typename T>
T getMatrixC() {
  T b(3, 3);
  b.set(0, 0, 1.0);
  b.set(1, 0, 2.0);
  b.set(2, 0, 3.0);
  b.set(0, 1, 2.0);
  b.set(1, 1, -1.0);
  b.set(2, 1, 1.0);
  b.set(0, 2, 4.0);
  b.set(1, 2, 3.0);
  b.set(2, 2, 2.0);

  return b;
}

// [ 0  1  3 ]
// [ -  2  4 ]
// [ -  -  5 ]
template <typename T>
T getSymMatrixA() {
  T a(3);
  a.set(0, 0, 0.0);
  a.set(0, 1, 1.0);
  a.set(1, 1, 2.0);
  a.set(0, 2, 3.0);
  a.set(1, 2, 4.0);
  a.set(2, 2, 5.0);

  return a;
}

// [ 0  1  2 ]
// [ -  3  4 ]
// [ -  -  5 ]
template <typename T>
T getSymMatrixB() {
  T b(3);
  b.set(0, 0, 0.0);
  b.set(0, 1, 1.0);
  b.set(1, 1, 3.0);
  b.set(0, 2, 2.0);
  b.set(1, 2, 4.0);
  b.set(2, 2, 5.0);

  return b;
}

//
// [ 0.937162
// [ 0.064600 0.233206
// [ 0.880494 0.228902 1.820559
// [ 0.633540 0.053748 1.080290 0.731896
template <typename T>
T getSymMatrixC() {
  T c(4);
  c.set(0, 0, 0.937162);
  c.set(1, 0, 0.064600);
  c.set(1, 1, 0.233206);
  c.set(2, 0, 0.880494);
  c.set(2, 1, 0.228902);
  c.set(2, 2, 1.820559);
  c.set(3, 0, 0.633540);
  c.set(3, 1, 0.053748);
  c.set(3, 2, 1.080290);
  c.set(3, 3, 0.731896);

  return c;
}

#endif  // MATRIX_COMMON_H
