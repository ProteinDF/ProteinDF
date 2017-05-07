#include <limits>
#include "gtest/gtest.h"
#include "TlSymmetricMatrix.h"

static const double EPS = 1.0E-10; //std::numeric_limits<double>::epsilon();
static const std::string smat_path = "temp_s.vct";

// 以下の要素を設定した行列を返す
// [ 0  1  3 ]
// [ -  2  4 ]
// [ -  -  5 ]
TlSymmetricMatrix getSymMatrixA()
{
    TlSymmetricMatrix a(3);
    a(0, 0) = 0.0;
    a(0, 1) = 1.0;
    a(1, 1) = 2.0;
    a(0, 2) = 3.0;
    a(1, 2) = 4.0;
    a(2, 2) = 5.0;
    
    return a;
}


// 以下の要素を設定した行列を返す
// [ 0  1  2 ]
// [ -  3  4 ]
// [ -  -  5 ]
TlSymmetricMatrix getSymMatrixB()
{
    TlSymmetricMatrix b(3);
    b(0, 0) = 0.0;
    b(0, 1) = 1.0;
    b(1, 1) = 3.0;
    b(0, 2) = 2.0;
    b(1, 2) = 4.0;
    b(2, 2) = 5.0;
    
    return b;
}


//
// [ 0.937162
// [ 0.064600 0.233206
// [ 0.880494 0.228902 1.820559
// [ 0.633540 0.053748 1.080290 0.731896
TlSymmetricMatrix getSymMatrixC()
{
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


TEST(TlSymmetricMatrix, constructer)
{
    TlSymmetricMatrix a(3);
    
    ASSERT_EQ(3, a.getNumOfRows());
    ASSERT_EQ(3, a.getNumOfCols());
    EXPECT_NEAR(0.0, a(0, 0), EPS);
    EXPECT_NEAR(0.0, a(0, 1), EPS);
    EXPECT_NEAR(0.0, a(0, 2), EPS);
    EXPECT_NEAR(0.0, a(1, 0), EPS);
    EXPECT_NEAR(0.0, a(1, 1), EPS);
    EXPECT_NEAR(0.0, a(1, 2), EPS);
    EXPECT_NEAR(0.0, a(2, 0), EPS);
    EXPECT_NEAR(0.0, a(2, 1), EPS);
    EXPECT_NEAR(0.0, a(2, 2), EPS);
}


TEST(TlSymmetricMatrix, pperaterRoundBracket)
{
    TlSymmetricMatrix a(3);
    
    // [ 0  -  - ]  
    // [ 1  2  - ]
    // [ 3  4  5 ]
    int t = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            a(i, j) = t;
            ++t;
        }
    }
    
    EXPECT_NEAR(0.0, a(0, 0), EPS);
    EXPECT_NEAR(1.0, a(0, 1), EPS);
    EXPECT_NEAR(3.0, a(0, 2), EPS);
    EXPECT_NEAR(1.0, a(1, 0), EPS);
    EXPECT_NEAR(2.0, a(1, 1), EPS);
    EXPECT_NEAR(4.0, a(1, 2), EPS);
    EXPECT_NEAR(3.0, a(2, 0), EPS);
    EXPECT_NEAR(4.0, a(2, 1), EPS);
    EXPECT_NEAR(5.0, a(2, 2), EPS);
}


TEST(TlSymmetricMatrix, copyConstructer)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix c(a);
    
    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_NEAR(0.0, c(0, 0), EPS);
    EXPECT_NEAR(1.0, c(0, 1), EPS);
    EXPECT_NEAR(3.0, c(0, 2), EPS);
    EXPECT_NEAR(1.0, c(1, 0), EPS);
    EXPECT_NEAR(2.0, c(1, 1), EPS);
    EXPECT_NEAR(4.0, c(1, 2), EPS);
    EXPECT_NEAR(3.0, c(2, 0), EPS);
    EXPECT_NEAR(4.0, c(2, 1), EPS);
    EXPECT_NEAR(5.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, convertFromTlVector1)
{
    // b =
    // { 0 1 2 3 4 5 }
    TlVector b(6);
    for (int i = 0; i < 6; ++i){
        b[i] = i;
    }
    
    TlSymmetricMatrix A(b, 3);
    
    // column oriented
    EXPECT_NEAR(0.0, A(0, 0), EPS);
    EXPECT_NEAR(1.0, A(0, 1), EPS);
    EXPECT_NEAR(2.0, A(1, 1), EPS);
    EXPECT_NEAR(3.0, A(0, 2), EPS);
    EXPECT_NEAR(4.0, A(1, 2), EPS);
    EXPECT_NEAR(5.0, A(2, 2), EPS);
}


TEST(TlSymmetricMatrix, convertFromTlVector2)
{
    TlSymmetricMatrix a = getSymMatrixA();
    
    TlVector v = a.getVector();
    
    EXPECT_NEAR(0.0, v[0], EPS);
    EXPECT_NEAR(1.0, v[1], EPS);
    EXPECT_NEAR(2.0, v[2], EPS);
    EXPECT_NEAR(3.0, v[3], EPS);
    EXPECT_NEAR(4.0, v[4], EPS);
    EXPECT_NEAR(5.0, v[5], EPS);
    
    TlSymmetricMatrix c(v, 3);
    
    EXPECT_NEAR(0.0, c(0, 0), EPS);
    EXPECT_NEAR(1.0, c(0, 1), EPS);
    EXPECT_NEAR(3.0, c(0, 2), EPS);
    EXPECT_NEAR(1.0, c(1, 0), EPS);
    EXPECT_NEAR(2.0, c(1, 1), EPS);
    EXPECT_NEAR(4.0, c(1, 2), EPS);
    EXPECT_NEAR(3.0, c(2, 0), EPS);
    EXPECT_NEAR(4.0, c(2, 1), EPS);
    EXPECT_NEAR(5.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_eq)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix c;
    
    c = a;
    
    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_NEAR(0.0, c(0, 0), EPS);
    EXPECT_NEAR(1.0, c(0, 1), EPS);
    EXPECT_NEAR(3.0, c(0, 2), EPS);
    EXPECT_NEAR(1.0, c(1, 0), EPS);
    EXPECT_NEAR(2.0, c(1, 1), EPS);
    EXPECT_NEAR(4.0, c(1, 2), EPS);
    EXPECT_NEAR(3.0, c(2, 0), EPS);
    EXPECT_NEAR(4.0, c(2, 1), EPS);
    EXPECT_NEAR(5.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_add)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix b = getSymMatrixB();
    
    TlSymmetricMatrix c = a + b;
    
//   0 1 3
//   1 2 4
//   3 4 5
    
//   0 1 2
//   1 3 4
//   2 4 5

    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_NEAR( 0.0, c(0, 0), EPS);
    EXPECT_NEAR( 2.0, c(0, 1), EPS);
    EXPECT_NEAR( 5.0, c(0, 2), EPS);
    EXPECT_NEAR( 2.0, c(1, 0), EPS);
    EXPECT_NEAR( 5.0, c(1, 1), EPS);
    EXPECT_NEAR( 8.0, c(1, 2), EPS);
    EXPECT_NEAR( 5.0, c(2, 0), EPS);
    EXPECT_NEAR( 8.0, c(2, 1), EPS);
    EXPECT_NEAR(10.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_iadd)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix b = getSymMatrixB();
    
    b += a;
    
    ASSERT_EQ(3, b.getNumOfRows());
    ASSERT_EQ(3, b.getNumOfCols());
    EXPECT_NEAR( 0.0, b(0, 0), EPS);
    EXPECT_NEAR( 2.0, b(0, 1), EPS);
    EXPECT_NEAR( 5.0, b(0, 2), EPS);
    EXPECT_NEAR( 2.0, b(1, 0), EPS);
    EXPECT_NEAR( 5.0, b(1, 1), EPS);
    EXPECT_NEAR( 8.0, b(1, 2), EPS);
    EXPECT_NEAR( 5.0, b(2, 0), EPS);
    EXPECT_NEAR( 8.0, b(2, 1), EPS);
    EXPECT_NEAR(10.0, b(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_mul)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix b = getSymMatrixB();
    
    TlMatrix c  = a * b;
    //c.print(std::cout);
    
    ASSERT_EQ(3, c.getNumOfRows());
    ASSERT_EQ(3, c.getNumOfCols());
    EXPECT_NEAR( 7.0, c(0, 0), EPS);
    EXPECT_NEAR(15.0, c(0, 1), EPS);
    EXPECT_NEAR(19.0, c(0, 2), EPS);
    EXPECT_NEAR(10.0, c(1, 0), EPS);
    EXPECT_NEAR(23.0, c(1, 1), EPS);
    EXPECT_NEAR(30.0, c(1, 2), EPS);
    EXPECT_NEAR(14.0, c(2, 0), EPS);
    EXPECT_NEAR(35.0, c(2, 1), EPS);
    EXPECT_NEAR(47.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, save)
{
    TlSymmetricMatrix a = getSymMatrixA();
    a.save("symmetric_matrix.a");
}


TEST(TlSymmetricMatrix, load)
{
    TlSymmetricMatrix a;
    a.load("symmetric_matrix.a");
    
    EXPECT_NEAR(0.0, a(0, 0), EPS);
    EXPECT_NEAR(1.0, a(0, 1), EPS);
    EXPECT_NEAR(3.0, a(0, 2), EPS);
    EXPECT_NEAR(1.0, a(1, 0), EPS);
    EXPECT_NEAR(2.0, a(1, 1), EPS);
    EXPECT_NEAR(4.0, a(1, 2), EPS);
    EXPECT_NEAR(3.0, a(2, 0), EPS);
    EXPECT_NEAR(4.0, a(2, 1), EPS);
    EXPECT_NEAR(5.0, a(2, 2), EPS);
}


TEST(TlSymmetricMatrix, inverse)
{
    TlSymmetricMatrix a = getSymMatrixA();
    TlSymmetricMatrix b = a;
    
    b.inverse();
    
    TlMatrix c = a * b;
    
    EXPECT_NEAR(1.0, c(0, 0), EPS);
    EXPECT_NEAR(0.0, c(0, 1), EPS);
    EXPECT_NEAR(0.0, c(0, 2), EPS);
    EXPECT_NEAR(0.0, c(1, 0), EPS);
    EXPECT_NEAR(1.0, c(1, 1), EPS);
    EXPECT_NEAR(0.0, c(1, 2), EPS);
    EXPECT_NEAR(0.0, c(2, 0), EPS);
    EXPECT_NEAR(0.0, c(2, 1), EPS);
    EXPECT_NEAR(1.0, c(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_mul1)
{
    TlSymmetricMatrix A = getSymMatrixA();
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
    
    EXPECT_NEAR(21.0, C(0, 0), EPS);
    EXPECT_NEAR(25.0, C(0, 1), EPS);
    EXPECT_NEAR(29.0, C(0, 2), EPS);
    EXPECT_NEAR(30.0, C(1, 0), EPS);
    EXPECT_NEAR(37.0, C(1, 1), EPS);
    EXPECT_NEAR(44.0, C(1, 2), EPS);
    EXPECT_NEAR(42.0, C(2, 0), EPS);
    EXPECT_NEAR(54.0, C(2, 1), EPS);
    EXPECT_NEAR(66.0, C(2, 2), EPS);
}


TEST(TlSymmetricMatrix, operator_multi2)
{
    TlSymmetricMatrix A = getSymMatrixA();
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
    
    EXPECT_NEAR( 7.0, C(0, 0), EPS);
    EXPECT_NEAR(10.0, C(0, 1), EPS);
    EXPECT_NEAR(14.0, C(0, 2), EPS);
    EXPECT_NEAR(19.0, C(1, 0), EPS);
    EXPECT_NEAR(31.0, C(1, 1), EPS);
    EXPECT_NEAR(50.0, C(1, 2), EPS);
    EXPECT_NEAR(31.0, C(2, 0), EPS);
    EXPECT_NEAR(52.0, C(2, 1), EPS);
    EXPECT_NEAR(86.0, C(2, 2), EPS);
}


TEST(TlSymmetricMatrix, imul1)
{
    TlSymmetricMatrix A = getSymMatrixA();
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
    
    TlMatrix C = A;
    C *= B;
    
    EXPECT_NEAR(21.0, C(0, 0), EPS);
    EXPECT_NEAR(25.0, C(0, 1), EPS);
    EXPECT_NEAR(29.0, C(0, 2), EPS);
    EXPECT_NEAR(30.0, C(1, 0), EPS);
    EXPECT_NEAR(37.0, C(1, 1), EPS);
    EXPECT_NEAR(44.0, C(1, 2), EPS);
    EXPECT_NEAR(42.0, C(2, 0), EPS);
    EXPECT_NEAR(54.0, C(2, 1), EPS);
    EXPECT_NEAR(66.0, C(2, 2), EPS);
}


TEST(TlSymmetricMatrix, testMultiEqual2)
{
    TlSymmetricMatrix A = getSymMatrixA();
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
    
    TlMatrix C = B;
    C *= A;

    EXPECT_NEAR( 7.0, C(0, 0), EPS);
    EXPECT_NEAR(10.0, C(0, 1), EPS);
    EXPECT_NEAR(14.0, C(0, 2), EPS);
    EXPECT_NEAR(19.0, C(1, 0), EPS);
    EXPECT_NEAR(31.0, C(1, 1), EPS);
    EXPECT_NEAR(50.0, C(1, 2), EPS);
    EXPECT_NEAR(31.0, C(2, 0), EPS);
    EXPECT_NEAR(52.0, C(2, 1), EPS);
    EXPECT_NEAR(86.0, C(2, 2), EPS);
}


TEST(TlSymmetricMatrix, dot)
{
    TlSymmetricMatrix A = getSymMatrixA();
    TlSymmetricMatrix B = getSymMatrixB();
    A.dot(B);
    
    EXPECT_NEAR( 0.0, A(0, 0), EPS);
    EXPECT_NEAR( 1.0, A(0, 1), EPS);
    EXPECT_NEAR( 6.0, A(0, 2), EPS);
    EXPECT_NEAR( 1.0, A(1, 0), EPS);
    EXPECT_NEAR( 6.0, A(1, 1), EPS);
    EXPECT_NEAR(16.0, A(1, 2), EPS);
    EXPECT_NEAR( 6.0, A(2, 0), EPS);
    EXPECT_NEAR(16.0, A(2, 1), EPS);
    EXPECT_NEAR(25.0, A(2, 2), EPS);
}


TEST(TlSymmetricMatrixTest, sum)
{
    TlSymmetricMatrix A = getSymMatrixA();
    double s = A.sum();
    
    EXPECT_NEAR(23.0, s, EPS);
}


TEST(TlSymmetricMatrix, choleskyDecomposition)
{
    TlSymmetricMatrix A = getSymMatrixC();
    //A.print(std::cout);

    //TlMatrix L = A.choleskyFactorization();
    TlMatrix L = A.choleskyFactorization2(1.0E-16);
    //L.print(std::cout);

    TlMatrix Lt = L;
    Lt.transpose();

    TlMatrix LL = L * Lt;
    //LL.print(std::cout);

    EXPECT_NEAR(A(0, 0), LL(0, 0), EPS);
    EXPECT_NEAR(A(0, 1), LL(0, 1), EPS);
    EXPECT_NEAR(A(0, 2), LL(0, 2), EPS);
    EXPECT_NEAR(A(0, 3), LL(0, 3), EPS);
    EXPECT_NEAR(A(1, 0), LL(1, 0), EPS);
    EXPECT_NEAR(A(1, 1), LL(1, 1), EPS);
    EXPECT_NEAR(A(1, 2), LL(1, 2), EPS);
    EXPECT_NEAR(A(1, 3), LL(1, 3), EPS);
    EXPECT_NEAR(A(2, 0), LL(2, 0), EPS);
    EXPECT_NEAR(A(2, 1), LL(2, 1), EPS);
    EXPECT_NEAR(A(2, 2), LL(2, 2), EPS);
    EXPECT_NEAR(A(2, 3), LL(2, 3), EPS);
    EXPECT_NEAR(A(3, 0), LL(3, 0), EPS);
    EXPECT_NEAR(A(3, 1), LL(3, 1), EPS);
    EXPECT_NEAR(A(3, 2), LL(3, 2), EPS);
    EXPECT_NEAR(A(3, 3), LL(3, 3), EPS);
}
