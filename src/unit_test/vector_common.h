#ifndef VECTOR_COMMON_H
#define VECTOR_COMMON_H

// {  0, 1, 2 } を返す
template <typename T>
T getVectorA() {
  T a(3);
  a.set(0, 0.0);
  a.set(1, 1.0);
  a.set(2, 2.0);

  return a;
}

// {  2, 4, 6 } を返す
template <typename T>
T getVectorB() {
  T b(3);
  b.set(0, 2.0);
  b.set(1, 4.0);
  b.set(2, 6.0);

  return b;
}

// {  1, 7, 3, 2, 0, 5 } を返す
template <typename T>
T getVectorC() {
  T c(6);
  c.set(0, 1.0);
  c.set(1, 7.0);
  c.set(2, 3.0);
  c.set(3, 2.0);
  c.set(4, 0.0);
  c.set(5, 5.0);

  return c;
}

template<typename T>
T getVector(int size) {
  T v(size);
  for (int i = 0; i < size; ++i) {
    v.set(i, double(i) * 0.01);
  }

  return v;
}

#endif  // VECTOR_COMMON_H
