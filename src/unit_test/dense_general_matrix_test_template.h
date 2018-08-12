#ifndef DENSE_GENERAL_MATRIX_TEST_TEMPLATE_H
#define DENSE_GENERAL_MATRIX_TEST_TEMPLATE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>
#include "gtest/gtest.h"
#include "matrix_common.h"
#include "tl_matrix_object.h"

template <typename T>
class DenseGeneralMatrixTest : public ::testing::Test {};

TYPED_TEST_CASE_P(DenseGeneralMatrixTest);

TYPED_TEST_P(DenseGeneralMatrixTest, doesConstructor) {
  TypeParam a = getMatrixA<TypeParam>();

  EXPECT_EQ(3, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesCopyConstructor) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam c = a;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesResize) {
  TypeParam a(2, 3);
  a.set(0, 1, 2.0);
  a.set(1, 2, -4.0);

  EXPECT_EQ(2, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(-4.0, a.get(1, 2));

  a.resize(4, 4);
  EXPECT_EQ(4, a.getNumOfRows());
  EXPECT_EQ(4, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(-4.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(2, 3));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 2));
  EXPECT_DOUBLE_EQ(0.0, a.get(3, 3));

  a.resize(2, 2);
  EXPECT_EQ(2, a.getNumOfRows());
  EXPECT_EQ(2, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(0.0, a.get(1, 1));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesOperatorEq) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam c;

  c = a;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesOperatorAdd) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam b = getMatrixB<TypeParam>();

  TypeParam c = a + b;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(4.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(8.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(4.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(8.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(12.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(8.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(12.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(16.0, c.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesOperatorIAdd) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam b = getMatrixB<TypeParam>();

  b += a;

  EXPECT_EQ(3, b.getNumOfRows());
  EXPECT_EQ(3, b.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, b.get(0, 0));
  EXPECT_DOUBLE_EQ(4.0, b.get(0, 1));
  EXPECT_DOUBLE_EQ(8.0, b.get(0, 2));
  EXPECT_DOUBLE_EQ(4.0, b.get(1, 0));
  EXPECT_DOUBLE_EQ(8.0, b.get(1, 1));
  EXPECT_DOUBLE_EQ(12.0, b.get(1, 2));
  EXPECT_DOUBLE_EQ(8.0, b.get(2, 0));
  EXPECT_DOUBLE_EQ(12.0, b.get(2, 1));
  EXPECT_DOUBLE_EQ(16.0, b.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesOperatorMultiMatrixAB) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam b = getMatrixB<TypeParam>();
  TypeParam c = a * b;

  EXPECT_EQ(3, c.getNumOfRows());
  EXPECT_EQ(3, c.getNumOfCols());
  EXPECT_DOUBLE_EQ(5.0, c.get(0, 0));
  EXPECT_DOUBLE_EQ(14.0, c.get(0, 1));
  EXPECT_DOUBLE_EQ(23.0, c.get(0, 2));
  EXPECT_DOUBLE_EQ(14.0, c.get(1, 0));
  EXPECT_DOUBLE_EQ(50.0, c.get(1, 1));
  EXPECT_DOUBLE_EQ(86.0, c.get(1, 2));
  EXPECT_DOUBLE_EQ(23.0, c.get(2, 0));
  EXPECT_DOUBLE_EQ(86.0, c.get(2, 1));
  EXPECT_DOUBLE_EQ(149.0, c.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesSum) {
  TypeParam a = getMatrixA<TypeParam>();

  double sum = a.sum();

  EXPECT_DOUBLE_EQ(36.0, sum);
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesTrace) {
  TypeParam a = getMatrixA<TypeParam>();

  double trace = a.trace();

  EXPECT_DOUBLE_EQ(12.0, trace);
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesMaxAbsoluteElement) {
  TypeParam a = getMatrixA<TypeParam>();

  TlMatrixObject::index_type row;
  TlMatrixObject::index_type col;
  double maxAbsoluteElement = a.getMaxAbsoluteElement(&row, &col);

  EXPECT_DOUBLE_EQ(8.0, maxAbsoluteElement);
  EXPECT_EQ(2, row);
  EXPECT_EQ(2, col);
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesTransposeInPlace) {
  TypeParam a = getMatrixA<TypeParam>();

  a.transposeInPlace();

  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(2.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(3.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(6.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(7.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesDotInPlace) {
  TypeParam a = getMatrixA<TypeParam>();
  TypeParam b = getMatrixB<TypeParam>();
  a.dotInPlace(b);

  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(12.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(3.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(16.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(35.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(12.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(35.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(64.0, a.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesSaveAndLoad) {
  static const std::string mat_save_load_path = "temp.gen.save_load.mat";

  TypeParam m = getMatrixA<TypeParam>();
  const bool isSaved = m.save(mat_save_load_path);
  EXPECT_EQ(isSaved, true);

  TypeParam a;
  const bool isLoaded = a.load(mat_save_load_path);
  EXPECT_EQ(isLoaded, true);

  EXPECT_EQ(3, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
}

TYPED_TEST_P(DenseGeneralMatrixTest, doesSaveAndLoadToHdf5) {
  static const std::string mat_h5_path = "temp.gen.save_load.h5";

#ifdef HAVE_HDF5
  TypeParam m = getMatrixA<TypeParam>();
  m.saveHdf5(mat_h5_path, "matrix_A");

  TypeParam a;
  a.loadHdf5(mat_h5_path, "matrix_A");
  EXPECT_EQ(3, a.getNumOfRows());
  EXPECT_EQ(3, a.getNumOfCols());
  EXPECT_DOUBLE_EQ(0.0, a.get(0, 0));
  EXPECT_DOUBLE_EQ(1.0, a.get(0, 1));
  EXPECT_DOUBLE_EQ(2.0, a.get(0, 2));
  EXPECT_DOUBLE_EQ(3.0, a.get(1, 0));
  EXPECT_DOUBLE_EQ(4.0, a.get(1, 1));
  EXPECT_DOUBLE_EQ(5.0, a.get(1, 2));
  EXPECT_DOUBLE_EQ(6.0, a.get(2, 0));
  EXPECT_DOUBLE_EQ(7.0, a.get(2, 1));
  EXPECT_DOUBLE_EQ(8.0, a.get(2, 2));
#else
  std::cerr << "HDF5 is not supported in this build." << std::endl;
#endif  // HAVE_HDF5
}

REGISTER_TYPED_TEST_CASE_P(DenseGeneralMatrixTest, doesConstructor,
                           doesCopyConstructor, doesResize, doesOperatorEq,
                           doesOperatorAdd, doesOperatorIAdd,
                           doesOperatorMultiMatrixAB, doesSum, doesTrace,
                           doesMaxAbsoluteElement, doesTransposeInPlace,
                           doesDotInPlace, doesSaveAndLoad,
                           doesSaveAndLoadToHdf5);

#endif  // DENSE_GENERAL_MATRIX_TEST_TEMPLATE_H
