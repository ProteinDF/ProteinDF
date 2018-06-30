#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "gtest/gtest.h"
#include "tl_dense_vector_blas.h"
#include "vector_common.h"

static const std::string vct_path = "temp.vct";
static const std::string h5_path = "temp.vct.h5";

TEST(TlVector_BLAS, constructer) {
  TlVector_BLAS a(3);

  EXPECT_EQ(TlVectorAbstract::size_type(3), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(0.0, a[1]);
  EXPECT_DOUBLE_EQ(0.0, a[2]);
}

TEST(TlVector_BLAS, copyConstructer) {
  // a = {  0, 1, 2 }
  TlVector_BLAS a(3);
  a[0] = 0.0;
  a[1] = 1.0;
  a[2] = 2.0;

  TlVector_BLAS b = a;
  EXPECT_EQ(TlVectorAbstract::size_type(3), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b[0]);
  EXPECT_DOUBLE_EQ(1.0, b[1]);
  EXPECT_DOUBLE_EQ(2.0, b[2]);
}

TEST(TlVector_BLAS, operatorCopy) {
  // a = { 0, 1, 2 }
  TlVector_BLAS a(3);
  for (int i = 0; i < 3; ++i) {
    a[i] = static_cast<double>(i);
  }

  TlVector_BLAS b;
  b = a;

  EXPECT_EQ(TlVectorAbstract::size_type(3), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b[0]);
  EXPECT_DOUBLE_EQ(1.0, b[1]);
  EXPECT_DOUBLE_EQ(2.0, b[2]);
}

TEST(TlVector_BLAS, resize) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();

  a.resize(5);
  EXPECT_EQ(TlVectorAbstract::size_type(5), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(1.0, a[1]);
  EXPECT_DOUBLE_EQ(2.0, a[2]);
  EXPECT_DOUBLE_EQ(0.0, a[3]);
  EXPECT_DOUBLE_EQ(0.0, a[4]);

  a.resize(2);
  EXPECT_EQ(TlVectorAbstract::size_type(2), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(1.0, a[1]);
}

TEST(TlVector_BLAS, operator_add) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  TlVector_BLAS b = getVectorB<TlVector_BLAS>();

  TlVector_BLAS c = a + b;

  EXPECT_EQ(TlVectorAbstract::size_type(3), c.getSize());
  EXPECT_DOUBLE_EQ(2.0, c[0]);
  EXPECT_DOUBLE_EQ(5.0, c[1]);
  EXPECT_DOUBLE_EQ(8.0, c[2]);
}

TEST(TlVector_BLAS, operator_iadd) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  TlVector_BLAS b = getVectorB<TlVector_BLAS>();

  b += a;

  EXPECT_DOUBLE_EQ(2.0, b[0]);
  EXPECT_DOUBLE_EQ(5.0, b[1]);
  EXPECT_DOUBLE_EQ(8.0, b[2]);
}

TEST(TlVector_BLAS, operator_sub) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  TlVector_BLAS b = getVectorB<TlVector_BLAS>();

  // {  0, 1, 2 }
  // {  2, 4, 6 }

  TlVector_BLAS c = a - b;

  EXPECT_EQ(TlVectorAbstract::size_type(3), c.getSize());
  EXPECT_DOUBLE_EQ(-2.0, c[0]);
  EXPECT_DOUBLE_EQ(-3.0, c[1]);
  EXPECT_DOUBLE_EQ(-4.0, c[2]);
}

TEST(TlVector_BLAS, operator_isub) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  TlVector_BLAS b = getVectorB<TlVector_BLAS>();

  b -= a;
  EXPECT_DOUBLE_EQ(2.0, b[0]);
  EXPECT_DOUBLE_EQ(3.0, b[1]);
  EXPECT_DOUBLE_EQ(4.0, b[2]);
}

TEST(TlVector_BLAS, operator_mul_double) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();

  TlVector_BLAS b = a * 2.0;
  EXPECT_EQ(TlVectorAbstract::size_type(3), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b[0]);
  EXPECT_DOUBLE_EQ(2.0, b[1]);
  EXPECT_DOUBLE_EQ(4.0, b[2]);

  TlVector_BLAS c = 3.0 * a;
  EXPECT_EQ(TlVectorAbstract::size_type(3), c.getSize());
  EXPECT_DOUBLE_EQ(0.0, c[0]);
  EXPECT_DOUBLE_EQ(3.0, c[1]);
  EXPECT_DOUBLE_EQ(6.0, c[2]);
}

TEST(TlVector_BLAS, operator_imul_double) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  a *= 2.0;

  EXPECT_EQ(TlVectorAbstract::size_type(3), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(2.0, a[1]);
  EXPECT_DOUBLE_EQ(4.0, a[2]);
}

TEST(TlVector_BLAS, operator_mul_TlVector_BLAS) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  TlVector_BLAS b = getVectorB<TlVector_BLAS>();

  // {  0, 1, 2 }
  // {  2, 4, 6 }

  double c = a * b;
  EXPECT_DOUBLE_EQ(16.0, c);
}

TEST(TlVector_BLAS, getMaxAbsoluteElement) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  const double c = a.getMaxAbsoluteElement();
  EXPECT_DOUBLE_EQ(2.0, c);
}

TEST(TlVector_BLAS, push_back) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();

  a.push_back(3.0);

  EXPECT_EQ(TlVectorAbstract::size_type(4), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(1.0, a[1]);
  EXPECT_DOUBLE_EQ(2.0, a[2]);
  EXPECT_DOUBLE_EQ(3.0, a[3]);
}

TEST(TlVector_BLAS, save) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  a.save(vct_path);
}

TEST(TlVector_BLAS, load) {
  TlVector_BLAS a;
  a.load(vct_path);

  EXPECT_EQ(TlVectorAbstract::size_type(3), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(1.0, a[1]);
  EXPECT_DOUBLE_EQ(2.0, a[2]);
}

#ifdef HAVE_HDF5
TEST(TlVector_BLAS, saveHdf5) {
  TlVector_BLAS a = getVectorA<TlVector_BLAS>();
  a.saveHdf5(h5_path, "vector");
}

TEST(TlVector_BLAS, loadHdf5) {
  TlVector_BLAS a;
  a.loadHdf5(h5_path, "vector");

  EXPECT_EQ(TlVectorAbstract::size_type(3), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a[0]);
  EXPECT_DOUBLE_EQ(1.0, a[1]);
  EXPECT_DOUBLE_EQ(2.0, a[2]);
}
#endif  // HAVE_HDF5
