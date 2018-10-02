#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "TlCommunicate.h"
#include "gtest/gtest.h"
#include "tl_dense_vector_lapack.h"
#include "tl_dense_vector_scalapack.h"
#include "vector_common.h"

static const std::string vct_path = "temp.scalapack.vct";
static const std::string vct_test_scalapack_load_path =
    "test.scalapack.load.vct";
static const std::string h5_path = "temp.scalapack.vct.h5";

TEST(TlDenseVector_Scalapack, constructor) {
  TlDenseVector_Scalapack a(100);

  EXPECT_EQ(TlDenseVector_Scalapack::size_type(100), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a.get(0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2));
}

TEST(TlDenseVector_Scalapack, copyConstructor) {
  TlDenseVector_Scalapack a(100);
  a.set(0, 0.0);
  a.set(20, 1.0);
  a.set(50, 2.0);

  TlDenseVector_Scalapack b = a;
  EXPECT_EQ(TlDenseVector_Scalapack::size_type(100), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b.get(0));
  EXPECT_DOUBLE_EQ(1.0, b.get(20));
  EXPECT_DOUBLE_EQ(2.0, b.get(50));
}

TEST(TlDenseVector_Scalapack, operatorCopy) {
  TlDenseVector_Scalapack a(100);
  a.set(0, 0.0);
  a.set(20, 1.0);
  a.set(50, 2.0);

  TlDenseVector_Scalapack b;
  b = a;

  EXPECT_EQ(TlDenseVector_Scalapack::size_type(100), b.getSize());
  EXPECT_DOUBLE_EQ(0.0, b.get(0));
  EXPECT_DOUBLE_EQ(1.0, b.get(20));
  EXPECT_DOUBLE_EQ(2.0, b.get(50));
}

TEST(TlDenseVector_Scalapack, resize1) {
  TlDenseVector_Scalapack a(100);
  a.set(0, 0.0);
  a.set(20, 1.0);
  a.set(50, 2.0);

  a.resize(200);
  EXPECT_EQ(TlDenseVector_Scalapack::size_type(200), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a.get(0));
  EXPECT_DOUBLE_EQ(1.0, a.get(20));
  EXPECT_DOUBLE_EQ(2.0, a.get(50));
}

TEST(TlDenseVector_Scalapack, resize2) {
  TlDenseVector_Scalapack a(100);
  a.set(0, 0.0);
  a.set(20, 1.0);
  a.set(50, 2.0);

  a.resize(50);
  EXPECT_EQ(TlDenseVector_Scalapack::size_type(50), a.getSize());
  EXPECT_DOUBLE_EQ(0.0, a.get(0));
  EXPECT_DOUBLE_EQ(1.0, a.get(20));
}

// TEST(TlDenseVector_Scalapack, operator_add) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   TlDenseVector_Scalapack b = getVectorB<TlDenseVector_Scalapack>();
//
//   TlDenseVector_Scalapack c = a + b;
//
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(2.0, c[0]);
//   EXPECT_DOUBLE_EQ(5.0, c[1]);
//   EXPECT_DOUBLE_EQ(8.0, c[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_iadd) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   TlDenseVector_Scalapack b = getVectorB<TlDenseVector_Scalapack>();
//
//   b += a;
//
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(5.0, b[1]);
//   EXPECT_DOUBLE_EQ(8.0, b[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_sub) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   TlDenseVector_Scalapack b = getVectorB<TlDenseVector_Scalapack>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   TlDenseVector_Scalapack c = a - b;
//
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(-2.0, c[0]);
//   EXPECT_DOUBLE_EQ(-3.0, c[1]);
//   EXPECT_DOUBLE_EQ(-4.0, c[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_isub) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   TlDenseVector_Scalapack b = getVectorB<TlDenseVector_Scalapack>();
//
//   b -= a;
//   EXPECT_DOUBLE_EQ(2.0, b[0]);
//   EXPECT_DOUBLE_EQ(3.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_mul_double) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//
//   TlDenseVector_Scalapack b = a * 2.0;
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), b.getSize());
//   EXPECT_DOUBLE_EQ(0.0, b[0]);
//   EXPECT_DOUBLE_EQ(2.0, b[1]);
//   EXPECT_DOUBLE_EQ(4.0, b[2]);
//
//   TlDenseVector_Scalapack c = 3.0 * a;
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), c.getSize());
//   EXPECT_DOUBLE_EQ(0.0, c[0]);
//   EXPECT_DOUBLE_EQ(3.0, c[1]);
//   EXPECT_DOUBLE_EQ(6.0, c[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_imul_double) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   a *= 2.0;
//
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(2.0, a[1]);
//   EXPECT_DOUBLE_EQ(4.0, a[2]);
// }
//
// TEST(TlDenseVector_Scalapack, operator_mul_TlDenseVector_Scalapack) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   TlDenseVector_Scalapack b = getVectorB<TlDenseVector_Scalapack>();
//
//   // {  0, 1, 2 }
//   // {  2, 4, 6 }
//
//   double c = a * b;
//   EXPECT_DOUBLE_EQ(16.0, c);
// }
//
// TEST(TlDenseVector_Scalapack, getMaxAbsoluteElement) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   const double c = a.getMaxAbsoluteElement();
//   EXPECT_DOUBLE_EQ(2.0, c);
// }
//
// TEST(TlDenseVector_Scalapack, push_back) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//
//   a.push_back(3.0);
//
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(4), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
//   EXPECT_DOUBLE_EQ(3.0, a[3]);
// }
//
TEST(TlDenseVector_Scalapack, save) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int size = 100;
  TlDenseVector_Lapack ref = getVector<TlDenseVector_Lapack>(size);
  TlDenseVector_Scalapack v = getVector<TlDenseVector_Scalapack>(size);
  v.save(vct_path);

  rComm.barrier();
  if (rComm.isMaster()) {
    TlDenseVector_Lapack a;
    a.load(vct_path);

    EXPECT_EQ(ref.getSize(), a.getSize());
    for (int i = 0; i < ref.getSize(); ++i) {
      EXPECT_DOUBLE_EQ(ref.get(i), a.get(i));
    }
  }
}

TEST(TlDenseVector_Scalapack, load) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int size = 100;
  TlDenseVector_Lapack v = getVector<TlDenseVector_Lapack>(size);
  if (rComm.isMaster()) {
    v.save(vct_test_scalapack_load_path);
  }

  rComm.barrier();
  TlDenseVector_Scalapack a;
  a.load(vct_test_scalapack_load_path);

  rComm.barrier();
  EXPECT_EQ(v.getSize(), a.getSize());
  for (int i = 0; i < v.getSize(); ++i) {
    double vi = v.get(i);
    double ai = a.get(i);
    if (rComm.isMaster()) {
      EXPECT_DOUBLE_EQ(vi, ai);
    }
  }
}

// #ifdef HAVE_HDF5
// TEST(TlDenseVector_Scalapack, saveHdf5) {
//   TlDenseVector_Scalapack a = getVectorA<TlDenseVector_Scalapack>();
//   a.saveHdf5(h5_path, "vector");
// }
//
// TEST(TlDenseVector_Scalapack, loadHdf5) {
//   TlDenseVector_Scalapack a;
//   a.loadHdf5(h5_path, "vector");
//
//   EXPECT_EQ(TlDenseVector_Scalapack::size_type(3), a.getSize());
//   EXPECT_DOUBLE_EQ(0.0, a[0]);
//   EXPECT_DOUBLE_EQ(1.0, a[1]);
//   EXPECT_DOUBLE_EQ(2.0, a[2]);
// }
// #endif  // HAVE_HDF5
