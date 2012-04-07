#include <limits>
#include "TlSymmetricMatrixTest.h"

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

const double TlSymmetricMatrixTest::threshold = std::numeric_limits<double>::epsilon() * 100;

// 以下の要素を設定した行列を返す
// [ 0  -  - ]
// [ 1  2  - ]
// [ 3  4  5 ]
TlSymmetricMatrix TlSymmetricMatrixTest::getMatrixA(){
    TlSymmetricMatrix a(3);
    a(0, 0) = 0.0;
    a(1, 0) = 1.0;
    a(1, 1) = 2.0;
    a(2, 0) = 3.0;
    a(2, 1) = 4.0;
    a(2, 2) = 5.0;
    
    return a;
}

// 以下の要素を設定した行列を返す
// [ 0  -  - ]
// [ 1  3  - ]
// [ 2  4  5 ]
TlSymmetricMatrix TlSymmetricMatrixTest::getMatrixB(){
    TlSymmetricMatrix b(3);
    b(0, 0) = 0.0;
    b(1, 0) = 1.0;
    b(1, 1) = 3.0;
    b(2, 0) = 2.0;
    b(2, 1) = 4.0;
    b(2, 2) = 5.0;
    
    return b;
}

//
// [ 0.937162
// [ 0.064600 0.233206
// [ 0.880494 0.228902 1.820559
// [ 0.633540 0.053748 1.080290 0.731896
TlSymmetricMatrix TlSymmetricMatrixTest::getMatrixC(){
    TlSymmetricMatrix c(4);
    c(0, 0) = 0.937162;
    c(1, 0) = 0.064600;
    c(1, 1) = 0.233206;
    c(2, 0) = 0.880494;
    c(2, 1) = 0.228902;
    c(2, 2) = 1.820559;
    c(3, 0) = 0.633540;
    c(3, 1) = 0.053748;
    c(3, 2) = 1.080290;
    c(3, 3) = 0.731896;
    
    return c;
}



void TlSymmetricMatrixTest::testConstructer(){
    TlSymmetricMatrix a(3);
    
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

void TlSymmetricMatrixTest::testOperaterRoundBracket(){
  TlSymmetricMatrix a(3);
  
  // [ 0  -  - ]  
  // [ 1  2  - ]
  // [ 3  4  5 ]
  int t = 0;
  for (int i=0; i<3; i++){
    for (int j=0; j<=i; j++){
      a(i, j) = t;
      t++;
    }
  }

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, a(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, a(2, 2), threshold);
}

void TlSymmetricMatrixTest::testCopyConstructer(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix c(a);

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, c(2, 2), threshold);
}

void TlSymmetricMatrixTest::testConvertFromTlVector1(){
  // b =
  // { 0 1 2 3 4 5 }
  TlVector b(6);
  for (int i = 0; i < 6; ++i){
    b[i] = i;
  }

  TlSymmetricMatrix A(b, 3);
  
  // column oriented
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, A(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, A(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, A(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, A(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, A(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, A(2, 2), threshold);
}

void TlSymmetricMatrixTest::testConvertFromTlVector2(){
  TlSymmetricMatrix a = this->getMatrixA();

  TlVector v = a.getVector();

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v[0], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v[1], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, v[2], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, v[3], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, v[4], threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, v[5], threshold);
  
  TlSymmetricMatrix c(v, 3);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, c(2, 2), threshold);
}

void TlSymmetricMatrixTest::testOperatorEqual(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix c;

  c = a;

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, c(2, 2), threshold);
}

void TlSymmetricMatrixTest::testOperatorPlus(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix b = this->getMatrixB();

  TlSymmetricMatrix c = a + b;

//   0 1 3
//   1 2 4
//   3 4 5

//   0 1 2
//   1 3 4
//   2 4 5

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, c(2, 2), threshold);
}

void TlSymmetricMatrixTest::testOperatorPlusEqual(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix b = this->getMatrixB();

  b += a;

  CPPUNIT_ASSERT_EQUAL(3, b.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, b.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, b(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, b(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, b(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, b(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, b(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, b(2, 2), threshold);
}

void TlSymmetricMatrixTest::testOperatorAsterisk(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix b = this->getMatrixB();

  TlMatrix c  = a * b;
  //c.print(std::cout);

  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(3, c.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.0, c(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0, c(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.0, c(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, c(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, c(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, c(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, c(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, c(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(47.0, c(2, 2), threshold);
}

void TlSymmetricMatrixTest::testSave(){
  TlSymmetricMatrix a = this->getMatrixA();
  a.save("symmetric_matrix.a");
}

void TlSymmetricMatrixTest::testLoad(){
  TlSymmetricMatrix a;
  a.load("symmetric_matrix.a");

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, a(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, a(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, a(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, a(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, a(2, 2), threshold);
}

void TlSymmetricMatrixTest::testInverse(){
  TlSymmetricMatrix a = this->getMatrixA();
  TlSymmetricMatrix b = a;

  b.inverse();

  TlMatrix c = a * b;

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


void TlSymmetricMatrixTest::testMulti1() {
  TlSymmetricMatrix A = this->getMatrixA();
  // [ 0  -  - ]
  // [ 1  2  - ]
  // [ 3  4  5 ]
  
  TlMatrix B(3, 3);
  B(0, 0) = 0.0;
  B(0, 1) = 1.0;
  B(0, 2) = 2.0;
  B(1, 0) = 3.0;
  B(1, 1) = 4.0;
  B(1, 2) = 5.0;
  B(2, 0) = 6.0;
  B(2, 1) = 7.0;
  B(2, 2) = 8.0;

  TlMatrix C = A * B;
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.0, C(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.0, C(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29.0, C(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, C(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(37.0, C(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(44.0, C(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0, C(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(54.0, C(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(66.0, C(2, 2), threshold);
}


void TlSymmetricMatrixTest::testMulti2() {
  TlSymmetricMatrix A = this->getMatrixA();
  // [ 0  -  - ]
  // [ 1  2  - ]
  // [ 3  4  5 ]
  
  TlMatrix B(3, 3);
  B(0, 0) = 0.0;
  B(0, 1) = 1.0;
  B(0, 2) = 2.0;
  B(1, 0) = 3.0;
  B(1, 1) = 4.0;
  B(1, 2) = 5.0;
  B(2, 0) = 6.0;
  B(2, 1) = 7.0;
  B(2, 2) = 8.0;

  TlMatrix C = B * A;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.0, C(0, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, C(0, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, C(0, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.0, C(1, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(31.0, C(1, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(50.0, C(1, 2), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(31.0, C(2, 0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(52.0, C(2, 1), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(86.0, C(2, 2), threshold);
}

void TlSymmetricMatrixTest::testDot() {
    TlSymmetricMatrix A = this->getMatrixA();
    TlSymmetricMatrix B = this->getMatrixB();
    A.dot(B);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, A(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, A(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, A(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, A(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(16.0, A(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, A(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(16.0, A(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(25.0, A(2, 2), threshold);
}


void TlSymmetricMatrixTest::testSum() {
    TlSymmetricMatrix A = this->getMatrixA();
    double s = A.sum();
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(23.0, s, threshold);
}

void TlSymmetricMatrixTest::testCholeskyDecomposition()
{
    TlSymmetricMatrix A = this->getMatrixC();
    //A.print(std::cout);

    //TlMatrix L = A.choleskyFactorization();
    TlMatrix L = A.choleskyFactorization2();
    //L.print(std::cout);

    TlMatrix Lt = L;
    Lt.transpose();

    TlMatrix LL = L * Lt;
    //LL.print(std::cout);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(0, 0), LL(0, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(0, 1), LL(0, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(0, 2), LL(0, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(0, 3), LL(0, 3), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(1, 0), LL(1, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(1, 1), LL(1, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(1, 2), LL(1, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(1, 3), LL(1, 3), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(2, 0), LL(2, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(2, 1), LL(2, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(2, 2), LL(2, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(2, 3), LL(2, 3), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(3, 0), LL(3, 0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(3, 1), LL(3, 1), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(3, 2), LL(3, 2), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(A(3, 3), LL(3, 3), threshold);
}
