#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "gtest/gtest.h"
#include "tl_dense_vector_eigen.h"
#include "tl_dense_vector_viennacl.h"
#include "vector_common.h"

static const std::string vct_path = "temp.vct";
static const std::string h5_path = "temp.vct.h5";

TEST(TlDenseVector_ViennaCL, constructer) {
  TlDenseVector_ViennaCL a(3);

  EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a.get(0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2));
}

TEST(TlDenseVector_ViennaCL, copyConstructer) {
  // a = {  0, 1, 2 }
  TlDenseVector_ViennaCL a(3);
  a.set(0, 0.0);
  a.set(1, 1.0);
  a.set(2, 2.0);

  TlDenseVector_ViennaCL b = a;
  EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b.get(0));
  EXPECT_DOUBLE_EQ(1.0, b.get(1));
  EXPECT_DOUBLE_EQ(2.0, b.get(2));
}

// TEST(TlDenseVector_ViennaCL, operatorCopy) {
//   // a = { 0, 1, 2 }
//   TlDenseVector_ViennaCL a(3);
//   for (int i = 0; i < 3; ++i) {
//     a[i] = static_cast<double>(i);
//   }
//
//   TlDenseVector_ViennaCL b;
//   b = a;
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(1.0, b[1]);
//   EXPECT_DOUBLE_EQ(2.0, b[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, resize) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//
//   a.resize(5);
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(5), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
//   EXPECT_DOUBLE_EQ(0.0, a[3]);
//   EXPECT_DOUBLE_EQ(0.0, a[4]);
//
//   a.resize(2);
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(2), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
// }

// TEST(TlDenseVector_ViennaCL, operator_add) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   TlDenseVector_ViennaCL b = getVectorB<TlDenseVector_ViennaCL>();
//
//   TlDenseVector_ViennaCL c = a + b;
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(2.0, c[0]);
//   EXPECT_DOUBLE_EQ(5.0, c[1]);
//   EXPECT_DOUBLE_EQ(8.0, c[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_iadd) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   TlDenseVector_ViennaCL b = getVectorB<TlDenseVector_ViennaCL>();
//
//   b += a;
//
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(5.0, b[1]);
//   EXPECT_DOUBLE_EQ(8.0, b[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_sub) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   TlDenseVector_ViennaCL b = getVectorB<TlDenseVector_ViennaCL>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   TlDenseVector_ViennaCL c = a - b;
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(-2.0, c[0]);
//   EXPECT_DOUBLE_EQ(-3.0, c[1]);
//   EXPECT_DOUBLE_EQ(-4.0, c[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_isub) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   TlDenseVector_ViennaCL b = getVectorB<TlDenseVector_ViennaCL>();
//
//   b -= a;
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(3.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_mul_double) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//
//   TlDenseVector_ViennaCL b = a * 2.0;
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(2.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
//
//   TlDenseVector_ViennaCL c = 3.0 * a;
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(0.0, c[0]);
//   EXPECT_DOUBLE_EQ(3.0, c[1]);
//   EXPECT_DOUBLE_EQ(6.0, c[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_imul_double) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   a *= 2.0;
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(2.0, a[1]);
//   EXPECT_DOUBLE_EQ(4.0, a[2]);
// }
//
// TEST(TlDenseVector_ViennaCL, operator_mul_TlDenseVector_ViennaCL) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   TlDenseVector_ViennaCL b = getVectorB<TlDenseVector_ViennaCL>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   double c = a * b;
//   EXPECT_DOUBLE_EQ(16.0, c);
// }
//
// TEST(TlDenseVector_ViennaCL, getMaxAbsoluteElement) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   const double c = a.getMaxAbsoluteElement();
//   EXPECT_DOUBLE_EQ(2.0, c);
// }
//
// TEST(TlDenseVector_ViennaCL, push_back) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//
//   a.push_back(3.0);
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(4), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
//   EXPECT_DOUBLE_EQ(3.0, a[3]);
// }
//
TEST(TlDenseVector_ViennaCL, reverse) {
  const int dim = 100;
  TlDenseVector_Eigen ref = getVector<TlDenseVector_Eigen>(dim);
  TlDenseVector_ViennaCL a = ref;
  a.reverse();

  EXPECT_EQ(dim, a.getSize());
  for (int i = 0; i < dim; ++i) {
    EXPECT_DOUBLE_EQ(ref.get(dim - i - 1), a.get(i));
  }
}

// TEST(TlDenseVector_ViennaCL, save) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   a.save(vct_path);
// }
//
// TEST(TlDenseVector_ViennaCL, load) {
//   TlDenseVector_ViennaCL a;
//   a.load(vct_path);
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
// }
//
// #ifdef HAVE_HDF5
// TEST(TlDenseVector_ViennaCL, saveHdf5) {
//   TlDenseVector_ViennaCL a = getVectorA<TlDenseVector_ViennaCL>();
//   a.saveHdf5(h5_path, "vector");
// }
//
// TEST(TlDenseVector_ViennaCL, loadHdf5) {
//   TlDenseVector_ViennaCL a;
//   a.loadHdf5(h5_path, "vector");
//
//   EXPECT_EQ(TlDenseVector_ViennaCL::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
// }
// #endif  // HAVE_HDF5
