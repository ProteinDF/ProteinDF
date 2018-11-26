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

template <typename T>
T getGeneralMatrix(int row, int col) {
  T m(row, col);
  int count = 1;
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      m.set(i, j, double(count));
      ++count;
    }
  }
  return m;
}

template <typename T>
T getGeneralSparseMatrix(int row, int col, int base = 7) {
  T m(row, col);
  int count = 1;
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      if ((count % base) == 0) {
        m.set(i, j, double(count));
      }
      ++count;
    }
  }
  return m;
}

template <typename T>
T getSymmetricMatrix(int dim) {
  T m(dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      m.set(i, j, double(i+1));
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

// [ 1  2  3  4]
// [ 5  6  7  8]
// [ 9 10 11 12]
template <typename T>
T getMatrixD() {
  T d(3, 4);
  d.set(0, 0, 1.0);  
  d.set(0, 1, 2.0);  
  d.set(0, 2, 3.0);  
  d.set(0, 3, 4.0);  
  d.set(1, 0, 5.0);  
  d.set(1, 1, 6.0);  
  d.set(1, 2, 7.0);  
  d.set(1, 3, 8.0);  
  d.set(2, 0, 9.0);  
  d.set(2, 1, 10.0);  
  d.set(2, 2, 11.0);  
  d.set(2, 3, 12.0);  

  return d;
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

template <typename T>
T getSymMatrixD() {
  T d(4);
  d.set(0, 0, 1.0);
  d.set(1, 0, 2.0);
  d.set(1, 1, 2.0);
  d.set(2, 0, 3.0);
  d.set(2, 1, 3.0);
  d.set(2, 2, 3.0);
  d.set(3, 0, 4.0);
  d.set(3, 1, 4.0);
  d.set(3, 2, 4.0);
  d.set(3, 3, 4.0);

  return d;
}


#endif  // MATRIX_COMMON_H
