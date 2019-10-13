#ifndef DENSE_VECTOR_TEST_TEMPLATE_H
#define DENSE_VECTOR_TEST_TEMPLATE_H

#include <iostream>
#include <string>
#include "gtest/gtest.h"
#include "vector_common.h"

template <typename T>
class DenseVectorTest : public ::testing::Test {};

TYPED_TEST_CASE_P(DenseVectorTest);

TYPED_TEST_P(DenseVectorTest, doesConstructor) {
    TypeParam a(3);

    EXPECT_EQ(3, a.getSize());
    EXPECT_DOUBLE_EQ(0.0, a.get(0));
    EXPECT_DOUBLE_EQ(0.0, a.get(1));
    EXPECT_DOUBLE_EQ(0.0, a.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesCopyConstructor) {
    TypeParam a(3);
    a.set(0, 0.0);
    a.set(1, 1.0);
    a.set(2, 2.0);

    TypeParam b = a;

    EXPECT_EQ(3, b.getSize());
    EXPECT_DOUBLE_EQ(0.0, b.get(0));
    EXPECT_DOUBLE_EQ(1.0, b.get(1));
    EXPECT_DOUBLE_EQ(2.0, b.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_eq) {
    const TypeParam a = getVectorA<TypeParam>();
    TypeParam b;

    b = a;

    EXPECT_EQ(3, b.getSize());
    EXPECT_DOUBLE_EQ(0.0, b.get(0));
    EXPECT_DOUBLE_EQ(1.0, b.get(1));
    EXPECT_DOUBLE_EQ(2.0, b.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_iadd) {
    const TypeParam a = getVectorA<TypeParam>();
    TypeParam b = getVectorB<TypeParam>();
    b += a;

    EXPECT_EQ(3, b.getSize());
    EXPECT_DOUBLE_EQ(2.0, b.get(0));
    EXPECT_DOUBLE_EQ(5.0, b.get(1));
    EXPECT_DOUBLE_EQ(8.0, b.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_isub) {
    const TypeParam a = getVectorA<TypeParam>();
    TypeParam b = getVectorB<TypeParam>();
    b -= a;

    EXPECT_EQ(3, b.getSize());
    EXPECT_DOUBLE_EQ(2.0, b.get(0));
    EXPECT_DOUBLE_EQ(3.0, b.get(1));
    EXPECT_DOUBLE_EQ(4.0, b.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_imul) {
    TypeParam a = getVectorA<TypeParam>();
    a *= 3.0;

    EXPECT_EQ(3, a.getSize());
    EXPECT_DOUBLE_EQ(0.0, a.get(0));
    EXPECT_DOUBLE_EQ(3.0, a.get(1));
    EXPECT_DOUBLE_EQ(6.0, a.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_add) {
    const TypeParam a = getVectorA<TypeParam>();
    const TypeParam b = getVectorB<TypeParam>();
    TypeParam c = a + b;

    EXPECT_EQ(3, c.getSize());
    EXPECT_DOUBLE_EQ(2.0, c.get(0));
    EXPECT_DOUBLE_EQ(5.0, c.get(1));
    EXPECT_DOUBLE_EQ(8.0, c.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesOperator_sub) {
    const TypeParam a = getVectorA<TypeParam>();
    const TypeParam b = getVectorB<TypeParam>();
    TypeParam c = a - b;

    EXPECT_EQ(3, c.getSize());
    EXPECT_DOUBLE_EQ(-2.0, c.get(0));
    EXPECT_DOUBLE_EQ(-3.0, c.get(1));
    EXPECT_DOUBLE_EQ(-4.0, c.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesSum) {
    TypeParam a = getVectorA<TypeParam>();
    const double sum = a.sum();

    EXPECT_DOUBLE_EQ(3.0, sum);
}

TYPED_TEST_P(DenseVectorTest, doesNorm) {
    const TypeParam a = getVectorA<TypeParam>();
    const double norm = a.norm();

    EXPECT_DOUBLE_EQ(2.23606797749979, norm);  // sqrt(5.0)
}

TYPED_TEST_P(DenseVectorTest, doesNorm2) {
    const TypeParam a = getVectorA<TypeParam>();
    const double norm2 = a.norm2();

    EXPECT_DOUBLE_EQ(5.0, norm2);
}

TYPED_TEST_P(DenseVectorTest, doesArgmax) {
    const TypeParam a = getVectorC<TypeParam>();
    const int argmax = a.argmax(0, a.getSize());

    EXPECT_EQ(1, argmax);  // a[1] == 7.0
}

TYPED_TEST_P(DenseVectorTest, doesSortByGreater) {
    TypeParam a = getVectorC<TypeParam>();
    a.sortByGreater();

    EXPECT_EQ(6, a.getSize());
    EXPECT_DOUBLE_EQ(7.0, a.get(0));
    EXPECT_DOUBLE_EQ(5.0, a.get(1));
    EXPECT_DOUBLE_EQ(3.0, a.get(2));
    EXPECT_DOUBLE_EQ(2.0, a.get(3));
    EXPECT_DOUBLE_EQ(1.0, a.get(4));
    EXPECT_DOUBLE_EQ(0.0, a.get(5));
}

TYPED_TEST_P(DenseVectorTest, doesDotInPlace) {
    TypeParam a = getVectorA<TypeParam>();
    const TypeParam b = getVectorB<TypeParam>();
    a.dotInPlace(b);

    EXPECT_EQ(3, a.getSize());
    EXPECT_DOUBLE_EQ(0.0, a.get(0));
    EXPECT_DOUBLE_EQ(4.0, a.get(1));
    EXPECT_DOUBLE_EQ(12.0, a.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesSaveAndLoad) {
    static const std::string vtr_save_load_path = "temp.dens.save_load.vtr";

    TypeParam v = getVectorA<TypeParam>();
    const bool isSaved = v.save(vtr_save_load_path);
    EXPECT_EQ(isSaved, true);

    TypeParam a;
    const bool isLoaded = a.load(vtr_save_load_path);
    EXPECT_EQ(isLoaded, true);

    EXPECT_EQ(3, a.getSize());
    EXPECT_DOUBLE_EQ(0.0, a.get(0));
    EXPECT_DOUBLE_EQ(1.0, a.get(1));
    EXPECT_DOUBLE_EQ(2.0, a.get(2));
}

TYPED_TEST_P(DenseVectorTest, doesTransformStdVector) {
    TypeParam a(3);
    a.set(0, 0.0);
    a.set(1, 1.0);
    a.set(2, 2.0);

    std::vector<double> b = a;

    EXPECT_EQ(3, b.size());
    EXPECT_DOUBLE_EQ(0.0, b[0]);
    EXPECT_DOUBLE_EQ(1.0, b[1]);
    EXPECT_DOUBLE_EQ(2.0, b[2]);
}

REGISTER_TYPED_TEST_CASE_P(DenseVectorTest, doesConstructor,
                           doesCopyConstructor, doesOperator_eq,
                           doesOperator_iadd, doesOperator_isub,
                           doesOperator_imul, doesOperator_add,
                           doesOperator_sub, doesSum, doesNorm, doesNorm2,
                           doesArgmax, doesSortByGreater, doesDotInPlace,
                           doesSaveAndLoad, doesTransformStdVector);

#endif  // DENSE_VECTOR_TEST_TEMPLATE_H
