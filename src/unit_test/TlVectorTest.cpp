#include <limits>
#include "gtest/gtest.h"
#include "TlVector.h"

static const double EPS = 1.0E-10; // std::numeric_limits<double>::epsilon();
static const std::string vct_path = "temp.vct";

// {  0, 1, 2 } を返す
TlVector getVectorA()
{ 
    TlVector a(3);
    a[0] = 0.0;
    a[1] = 1.0;
    a[2] = 2.0;
    
    return a;
}


// {  2, 4, 6 } を返す
TlVector getVectorB()
{
    TlVector b(3);
    b[0] = 2.0;
    b[1] = 4.0;
    b[2] = 6.0;
    
    return b;
}


TEST(TlVector, constructer)
{
    TlVector a(3);
    
    EXPECT_EQ(TlVector::size_type(3), a.getSize());
    EXPECT_NEAR(0.0, a[0], EPS);
    EXPECT_NEAR(0.0, a[1], EPS);
    EXPECT_NEAR(0.0, a[2], EPS);
}


TEST(TlVector, copyConstructer)
{
    // a = {  0, 1, 2 }
    TlVector a(3);
    a[0] = 0.0;
    a[1] = 1.0;
    a[2] = 2.0;

    TlVector b = a;
    EXPECT_EQ(TlVector::size_type(3), b.getSize());
    EXPECT_NEAR( 0.0, b[0], EPS);
    EXPECT_NEAR( 1.0, b[1], EPS);
    EXPECT_NEAR( 2.0, b[2], EPS);
}


TEST(TlVector, operatorCopy)
{
    // a = { 0, 1, 2 }
    TlVector a(3);
    for (int i = 0; i < 3; ++i) {
        a[i] = static_cast<double>(i);
    }
    
    TlVector b;
    b = a;

    EXPECT_EQ(TlVector::size_type(3), b.getSize());
    EXPECT_NEAR( 0.0, b[0], EPS);
    EXPECT_NEAR( 1.0, b[1], EPS);
    EXPECT_NEAR( 2.0, b[2], EPS);
}


TEST(TlVector, resize)
{
    TlVector a = getVectorA();
    
    a.resize(5);
    EXPECT_EQ(TlVector::size_type(5), a.getSize());
    EXPECT_NEAR( 0.0, a[0], EPS);
    EXPECT_NEAR( 1.0, a[1], EPS);
    EXPECT_NEAR( 2.0, a[2], EPS);
    EXPECT_NEAR( 0.0, a[3], EPS);
    EXPECT_NEAR( 0.0, a[4], EPS);
    
    a.resize(2);
    EXPECT_EQ(TlVector::size_type(2), a.getSize());
    EXPECT_NEAR( 0.0, a[0], EPS);
    EXPECT_NEAR( 1.0, a[1], EPS);
}


TEST(TlVector, operator_add)
{
    TlVector a = getVectorA();
    TlVector b = getVectorB();
    
    TlVector c = a + b;
    
    EXPECT_EQ(TlVector::size_type(3), c.getSize());
    EXPECT_NEAR( 2.0, c[0], EPS);
    EXPECT_NEAR( 5.0, c[1], EPS);
    EXPECT_NEAR( 8.0, c[2], EPS);
}


TEST(TlVector, operator_iadd)
{
    TlVector a = getVectorA();
    TlVector b = getVectorB();
    
    b += a;
    
    EXPECT_NEAR( 2.0, b[0], EPS);
    EXPECT_NEAR( 5.0, b[1], EPS);
    EXPECT_NEAR( 8.0, b[2], EPS);
}


TEST(TlVector, operator_sub)
{
    TlVector a = getVectorA();
    TlVector b = getVectorB();
    
    // {  0, 1, 2 }
    // {  2, 4, 6 }
    
    TlVector c = a - b;
    
    EXPECT_EQ(TlVector::size_type(3), c.getSize());
    EXPECT_NEAR( -2.0, c[0], EPS);
    EXPECT_NEAR( -3.0, c[1], EPS);
    EXPECT_NEAR( -4.0, c[2], EPS);
}


TEST(TlVector, operator_isub)
{
    TlVector a = getVectorA();
    TlVector b = getVectorB();
    
    b -= a;
    
    EXPECT_NEAR(  2.0, b[0], EPS);
    EXPECT_NEAR(  3.0, b[1], EPS);
    EXPECT_NEAR(  4.0, b[2], EPS);
}


TEST(TlVector, operator_mul_double)
{
    TlVector a = getVectorA();
    
    TlVector b = a * 2.0;
    
    EXPECT_EQ(TlVector::size_type(3), b.getSize());
    EXPECT_NEAR( 0.0, b[0], EPS);
    EXPECT_NEAR( 2.0, b[1], EPS);
    EXPECT_NEAR( 4.0, b[2], EPS);
    
    TlVector c = 3.0 * a;
    
    EXPECT_EQ(TlVector::size_type(3), c.getSize());
    EXPECT_NEAR( 0.0, c[0], EPS);
    EXPECT_NEAR( 3.0, c[1], EPS);
    EXPECT_NEAR( 6.0, c[2], EPS);
}


TEST(TlVector, operator_imul_double)
{
    TlVector a = getVectorA();
    
    a *= 2.0;
    
    EXPECT_EQ(TlVector::size_type(3), a.getSize());
    EXPECT_NEAR( 0.0, a[0], EPS);
    EXPECT_NEAR( 2.0, a[1], EPS);
    EXPECT_NEAR( 4.0, a[2], EPS);
}


TEST(TlVector, operator_mul_TlVector)
{
    TlVector a = getVectorA();
    TlVector b = getVectorB();
    
    // {  0, 1, 2 }
    // {  2, 4, 6 }
    
    double c = a * b;
    EXPECT_NEAR( 16.0, c, EPS);
}


TEST(TlVector, getMaxAbsoluteElement)
{
    TlVector a = getVectorA();
    
    double c = a.getMaxAbsoluteElement();
    EXPECT_NEAR( 2.0, c, EPS);
}


TEST(TlVector, push_back)
{
  TlVector a = getVectorA();

  a.push_back(3.0);

  EXPECT_EQ(TlVector::size_type(4), a.getSize());
  EXPECT_NEAR( 0.0, a[0], EPS);
  EXPECT_NEAR( 1.0, a[1], EPS);
  EXPECT_NEAR( 2.0, a[2], EPS);
  EXPECT_NEAR( 3.0, a[3], EPS);
}


TEST(TlVector, save)
{
    TlVector a = getVectorA();
    
    a.save(vct_path);
}


TEST(TlVector, load)
{
    TlVector a;
    a.load(vct_path);
    
    EXPECT_EQ(TlVector::size_type(3), a.getSize());
    EXPECT_NEAR( 0.0, a[0], EPS);
    EXPECT_NEAR( 1.0, a[1], EPS);
    EXPECT_NEAR( 2.0, a[2], EPS);
}

