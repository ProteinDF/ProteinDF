#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "gtest/gtest.h"
#include "tl_dense_vector_eigen.h"
#include "vector_common.h"

static const std::string vct_path = "temp.vct";
static const std::string h5_path = "temp.vct.h5";

TEST(TlDenseVector_Eigen, constructer) {
    TlDenseVector_Eigen a(3);

    EXPECT_EQ(TlDenseVector_Eigen::size_type(3), a.getSize());
    EXPECT_DOUBLE_EQ(0.0, a.get(0));
    EXPECT_DOUBLE_EQ(0.0, a.get(1));
    EXPECT_DOUBLE_EQ(0.0, a.get(2));
}

// TEST(TlDenseVector_Eigen, copyConstructer) {
//   // a = {  0, 1, 2 }
//   TlDenseVector_Eigen a(3);
//   a[0] = 0.0;
//   a[1] = 1.0;
//   a[2] = 2.0;
//
//   TlDenseVector_Eigen b = a;
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(1.0, b[1]);
//   EXPECT_DOUBLE_EQ(2.0, b[2]);
// }
//
// TEST(TlDenseVector_Eigen, operatorCopy) {
//   // a = { 0, 1, 2 }
//   TlDenseVector_Eigen a(3);
//   for (int i = 0; i < 3; ++i) {
//     a[i] = static_cast<double>(i);
//   }
//
//   TlDenseVector_Eigen b;
//   b = a;
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(1.0, b[1]);
//   EXPECT_DOUBLE_EQ(2.0, b[2]);
// }
//
// TEST(TlDenseVector_Eigen, resize) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//
//   a.resize(5);
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(5), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
//   EXPECT_DOUBLE_EQ(0.0, a[3]);
//   EXPECT_DOUBLE_EQ(0.0, a[4]);
//
//   a.resize(2);
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(2), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
// }

// TEST(TlDenseVector_Eigen, operator_add) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   TlDenseVector_Eigen b = getVectorB<TlDenseVector_Eigen>();
//
//   TlDenseVector_Eigen c = a + b;
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(2.0, c[0]);
//   EXPECT_DOUBLE_EQ(5.0, c[1]);
//   EXPECT_DOUBLE_EQ(8.0, c[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_iadd) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   TlDenseVector_Eigen b = getVectorB<TlDenseVector_Eigen>();
//
//   b += a;
//
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(5.0, b[1]);
//   EXPECT_DOUBLE_EQ(8.0, b[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_sub) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   TlDenseVector_Eigen b = getVectorB<TlDenseVector_Eigen>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   TlDenseVector_Eigen c = a - b;
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(-2.0, c[0]);
//   EXPECT_DOUBLE_EQ(-3.0, c[1]);
//   EXPECT_DOUBLE_EQ(-4.0, c[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_isub) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   TlDenseVector_Eigen b = getVectorB<TlDenseVector_Eigen>();
//
//   b -= a;
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(3.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_mul_double) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//
//   TlDenseVector_Eigen b = a * 2.0;
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(2.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
//
//   TlDenseVector_Eigen c = 3.0 * a;
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(0.0, c[0]);
//   EXPECT_DOUBLE_EQ(3.0, c[1]);
//   EXPECT_DOUBLE_EQ(6.0, c[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_imul_double) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   a *= 2.0;
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(2.0, a[1]);
//   EXPECT_DOUBLE_EQ(4.0, a[2]);
// }
//
// TEST(TlDenseVector_Eigen, operator_mul_TlDenseVector_Eigen) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   TlDenseVector_Eigen b = getVectorB<TlDenseVector_Eigen>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   double c = a * b;
//   EXPECT_DOUBLE_EQ(16.0, c);
// }
//
// TEST(TlDenseVector_Eigen, getMaxAbsoluteElement) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   const double c = a.getMaxAbsoluteElement();
//   EXPECT_DOUBLE_EQ(2.0, c);
// }
//
// TEST(TlDenseVector_Eigen, push_back) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//
//   a.push_back(3.0);
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(4), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
//   EXPECT_DOUBLE_EQ(3.0, a[3]);
// }
//
// TEST(TlDenseVector_Eigen, save) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   a.save(vct_path);
// }
//
// TEST(TlDenseVector_Eigen, load) {
//   TlDenseVector_Eigen a;
//   a.load(vct_path);
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
// }
//
// #ifdef HAVE_HDF5
// TEST(TlDenseVector_Eigen, saveHdf5) {
//   TlDenseVector_Eigen a = getVectorA<TlDenseVector_Eigen>();
//   a.saveHdf5(h5_path, "vector");
// }
//
// TEST(TlDenseVector_Eigen, loadHdf5) {
//   TlDenseVector_Eigen a;
//   a.loadHdf5(h5_path, "vector");
//
//   EXPECT_EQ(TlDenseVector_Eigen::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
// }
// #endif  // HAVE_HDF5
