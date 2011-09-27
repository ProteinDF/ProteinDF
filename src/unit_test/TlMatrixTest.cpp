#include <limits>
#include <iostream>
#include "TlMatrixTest.h"

const double TlMatrixTest::threshold = std::numeric_limits<double>::epsilon();


// 以下の要素を設定した行列を返す
// [ 0  1  2 ]
// [ 3  4  5 ]
// [ 6  7  8 ]
TlMatrix TlMatrixTest::getMatrixA(){
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

// 以下の要素を設定した行列を返す
// [ 0  3  6 ]
// [ 1  4  7 ]
// [ 2  5  8 ]
TlMatrix TlMatrixTest::getMatrixB(){
    TlMatrix b(3, 3);
    b(0, 0) = 0.0;
    b(1, 0) = 1.0;
    b(2, 0) = 2.0;
    b(0, 1) = 3.0;
    b(1, 1) = 4.0;
    b(2, 1) = 5.0;
    b(0, 2) = 6.0;
    b(1, 2) = 7.0;
    b(2, 2) = 8.0;
    
    return b;
}

// 以下の要素を設定した行列を返す
// [ 1   2  3 ]
// [ 2  -1  1 ]
// [ 4   3  2 ]
TlMatrix TlMatrixTest::getMatrixC(){
    TlMatrix b(3, 3);
    b(0, 0) = 1.0;
    b(1, 0) = 2.0;
    b(2, 0) = 3.0;
    b(0, 1) = 2.0;
    b(1, 1) =-1.0;
    b(2, 1) = 1.0;
    b(0, 2) = 4.0;
    b(1, 2) = 3.0;
    b(2, 2) = 2.0;
    
    return b;
}

void TlMatrixTest::testConstructer(){
    TlMatrix a(3, 3);
    
    CPPUNIT_ASSERT_EQUAL(3, a.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(3, a.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(2, 2), threshold);
}

void TlMatrixTest::testCopyConstructer(){
    TlMatrix a = this->getMatrixA();
    TlMatrix c(a);
    
    CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, c(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, c(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, c(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, c(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, c(2, 2), threshold);
}

// void TlMatrixTest::testConvertFromVector(){
//   TlMatrix A(2, 2);
  
//   // b =
//   // { 0 1 2 3 4 5 6 7 8 }
//   TlVector b(9);
//   for (int i = 0; i < 9; ++i){
//     b[i] = i;
//   }

//   A.convertFromVector(3, 3, b);
  
//   // column oriented
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, A(0, 0), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, A(1, 0), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, A(2, 0), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, A(0, 1), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, A(1, 1), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, A(2, 1), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, A(0, 2), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, A(1, 2), threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, A(2, 2), threshold);
// }

void TlMatrixTest::testOperatorEqual(){
  TlMatrix a = this->getMatrixA();
  TlMatrix c;

  c = a;

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, c(2, 2), threshold);
}

void TlMatrixTest::testOperatorPlus(){
  TlMatrix a = this->getMatrixA();
  TlMatrix b = this->getMatrixB();

  TlMatrix c = a + b;

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(16.0, c(2, 2), threshold);
}

void TlMatrixTest::testOperatorPlusEqual(){
  TlMatrix a = this->getMatrixA();
  TlMatrix b = this->getMatrixB();

  b += a;

  CPPUNIT_ASSERT_EQUAL(3, b.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, b.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, b(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, b(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, b(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, b(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(16.0, b(2, 2), threshold);
}

void TlMatrixTest::testSave(){
  TlMatrix a = this->getMatrixA();
  a.save("normal_matrix.a");
}

void TlMatrixTest::testLoad(){
  TlMatrix a;
  a.load("normal_matrix.a");

  CPPUNIT_ASSERT_EQUAL(3, a.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, a.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, a(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, a(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, a(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, a(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, a(2, 2), threshold);
}

void TlMatrixTest::testInverse(){
    TlMatrix a = this->getMatrixC();
    TlMatrix b = a;
    
    b.inverse();

    //a.print(std::cout);
    //b.print(std::cout);
    
    TlMatrix c = a * b;

    //c.print(std::cout);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(2, 2), threshold);
}

void TlMatrixTest::testOperatorMul_AB() {
    TlMatrix a = this->getMatrixA();
    TlMatrix b = this->getMatrixB();
    TlMatrix c  = a * b;

    CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, c(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, c(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, c(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, c(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(50.0, c(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(86.0, c(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, c(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(86.0, c(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(149.0, c(2, 2), threshold);
}

void TlMatrixTest::testOperatorMul_AX() {
    TlMatrix a = this->getMatrixA();
    TlVector x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    std::cerr << "flag1" << std::endl;
    TlVector z = a * x;
    std::cerr << "flag2" << std::endl;
    
    CPPUNIT_ASSERT_EQUAL(3, (int)z.getSize());
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, z[0], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26.0, z[1], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(44.0, z[2], threshold);
}


void TlMatrixTest::testOperatorMul_XA() {
    TlMatrix a = this->getMatrixA();
    TlVector x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    TlVector z = x * a;
    
    CPPUNIT_ASSERT_EQUAL(3, (int)z.getSize());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(24.0, z[0], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, z[1], threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(36.0, z[2], threshold);
}


void TlMatrixTest::testDot() {
    TlMatrix a = this->getMatrixA();
    TlMatrix b = this->getMatrixB();
    a.dot(b);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, a(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(16.0, a(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, a(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, a(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, a(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(64.0, a(2, 2), threshold);
}


void TlMatrixTest::testSum() {
    TlMatrix a = this->getMatrixA();
    double s = a.sum();
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(36.0, s, threshold);
}
