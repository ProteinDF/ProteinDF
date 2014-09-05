#include <limits>
#include "TlVectorTest.h"

// CPPUNIT_ASSERT( condition );
// conditionが偽(false,0)であったとき、失敗します。
//
// CPPUNIT_ASSERT_MESSAGE( message, condition );
// conditionが偽であったとき、失敗します。このときmessageを出力します。
//
// CPPUNIT_FAIL( message );
// 必ず失敗します。messageを出力します。
//
// CPPUNIT_ASSERT_EQUAL( expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// 得られた結果actualが期待する値expectedでなかったとき、すなわちexpected != actualのときに失敗します。このときmessageを出力します。
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// 得られた結果actualと期待する値expectedとの差がdeltaより大きいとき、失敗します。

const double TlVectorTest::threshold = std::numeric_limits<double>::epsilon();

// {  0, 1, 2 } を返す
TlVector TlVectorTest::getVectorA(){
  TlVector a(3);
  a[0] = 0.0;
  a[1] = 1.0;
  a[2] = 2.0;

  return a;
}

// {  2, 4, 6 } を返す
TlVector TlVectorTest::getVectorB(){
  TlVector b(3);
  b[0] = 2.0;
  b[1] = 4.0;
  b[2] = 6.0;

  return b;
}

void TlVectorTest::testConstructer(){
  TlVector a(3);

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a[2], threshold);
}

void TlVectorTest::testCopyConstructer(){
  // a = {  0, 1, 2 }
  TlVector a(3);
  a[0] = 0.0;
  a[1] = 1.0;
  a[2] = 2.0;

  TlVector b = a;
  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), b.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, b[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b[2], threshold);
}

// void TlVectorTest::testConvertFromMatrix(){
//   // a = { 0, 1, 2 }
//   TlVector a(3);
//   for (long i=0; i<3; i++){
//     a[i] = static_cast<double>(i);
//   }

//   // B = 
//   // 1 2 3
//   // 4 5 6
//   // 7 8 9
//   TlMatrix B(3, 3);
//   B(0, 0) = 1;
//   B(0, 1) = 2;
//   B(0, 2) = 3;
//   B(1, 0) = 4;
//   B(1, 1) = 5;
//   B(1, 2) = 6;
//   B(2, 0) = 7;
//   B(2, 1) = 8;
//   B(2, 2) = 9;

//   a.convertFromMatrix(B);

//   CPPUNIT_ASSERT_EQUAL((long)9, a.getSize());

//   // coulumn oriented
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[0], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, a[1], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.0, a[2], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[3], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, a[4], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, a[5], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a[6], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, a[7], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 9.0, a[8], threshold);
 
// }

// void TlVectorTest::testConvertFromMatrix_Symmetric(){
//   // a = { 0, 1, 2 }
//   TlVector a(3);
//   for (long i=0; i<3; i++){
//     a[i] = static_cast<double>(i);
//   }

//   // B = 
//   // 1 2 4
//   // 2 3 5
//   // 4 5 6
//   TlMatrix_Symmetric B(3);
//   B(0, 0) = 1;
//   B(1, 0) = 2;
//   B(1, 1) = 3;
//   B(2, 0) = 4;
//   B(2, 1) = 5;
//   B(2, 2) = 6;

//   a.convertFromMatrix(B);

//   CPPUNIT_ASSERT_EQUAL((long)6, a.getSize());

//   // coulumn oriented
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[0], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[1], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, a[2], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a[3], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, a[4], threshold);
//   CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, a[5], threshold);
  
// }

void TlVectorTest::testOperatorCopy(){
  // a = { 0, 1, 2 }
  TlVector a(3);
  for (long i=0; i<3; i++){
    a[i] = static_cast<double>(i);
  }

  TlVector b;
  b = a;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), b.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, b[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b[2], threshold);
}

void TlVectorTest::testResize(){
  TlVector a = this->getVectorA();

  a.resize(5);
  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(5), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[2], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[3], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[4], threshold);

  a.resize(2);
  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(2), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[1], threshold);
}

void TlVectorTest::testOperatorPlus(){
  TlVector a = this->getVectorA();
  TlVector b = this->getVectorB();

  TlVector c = a + b;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), c.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, c[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, c[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c[2], threshold);
}

void TlVectorTest::testOperatorPlusEqual(){
  TlVector a = this->getVectorA();
  TlVector b = this->getVectorB();
  
  b += a;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, b[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b[2], threshold);
}

void TlVectorTest::testOperatorMinus(){
  TlVector a = this->getVectorA();
  TlVector b = this->getVectorB();

  // {  0, 1, 2 }
  // {  2, 4, 6 }

  TlVector c = a - b;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), c.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -2.0, c[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -3.0, c[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -4.0, c[2], threshold);
}

void TlVectorTest::testOperatorMinusEqual(){
  TlVector a = this->getVectorA();
  TlVector b = this->getVectorB();

  b -= a;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(  2.0, b[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(  3.0, b[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(  4.0, b[2], threshold);
}

void TlVectorTest::testOperatorMultipleDouble(){
  TlVector a = this->getVectorA();

  TlVector b = a * 2.0;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), b.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, b[2], threshold);

  TlVector c = 3.0 * a;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), c.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, c[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, c[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, c[2], threshold);

}

void TlVectorTest::testOperatorMultipleEqualDouble(){
  TlVector a = this->getVectorA();

  a *= 2.0;

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, a[2], threshold);
}

void TlVectorTest::testOperatorMultipleVector(){
  TlVector a = this->getVectorA();
  TlVector b = this->getVectorB();
  
  // {  0, 1, 2 }
  // {  2, 4, 6 }

  double c = a * b;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 16.0, c, threshold);
}

void TlVectorTest::testGetMaxAbsoluteElement(){
  TlVector a = this->getVectorA();

  double c = a.getMaxAbsoluteElement();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, c, threshold);
}

void TlVectorTest::testAdd(){
  TlVector a = this->getVectorA();

  a.push_back(3.0);

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(4), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[0], 1.0E-16);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[1], 1.0E-16);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[2], 1.0E-16);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, a[3], 1.0E-16);
}

void TlVectorTest::testSave(){
  TlVector a = this->getVectorA();
  
  a.save("tlvector.a");
}

void TlVectorTest::testLoad(){
  TlVector a;
  a.load("tlvector.a");

  CPPUNIT_ASSERT_EQUAL(TlVector::size_type(3), a.getSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, a[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, a[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, a[2], threshold);
}
