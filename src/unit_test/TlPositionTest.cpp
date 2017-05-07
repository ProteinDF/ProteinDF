#include <vector>
#include <string>
#include <limits>
#include "gtest/gtest.h"
#include "TlPosition.h"

static const double EPS = 1.0E-10; // std::numeric_limits<double>::epsilon();

TEST(TlPosition, constructer)
{
    TlPosition a;
    EXPECT_NEAR(0.0, a.x(), EPS);
    EXPECT_NEAR(0.0, a.y(), EPS);
    EXPECT_NEAR(0.0, a.z(), EPS);
    
    TlPosition b(1.0, 2.0, 3.0);
    EXPECT_NEAR(1.0, b.x(), EPS);
    EXPECT_NEAR(2.0, b.y(), EPS);
    EXPECT_NEAR(3.0, b.z(), EPS);
}


TEST(TlPosition, copyConstructer)
{
    TlPosition a(1.0, 2.0, 3.0);
    TlPosition b(a);
    
    EXPECT_NEAR(1.0, b.x(), EPS);
    EXPECT_NEAR(2.0, b.y(), EPS);
    EXPECT_NEAR(3.0, b.z(), EPS);
}


TEST(TlPositionTest, pperator_eq)
{
    TlPosition a(1.0, 2.0, 3.0);
    TlPosition b(2.0, 3.0, 4.0);
    
    b = a;
    EXPECT_NEAR(1.0, b.x(), EPS);
    EXPECT_NEAR(2.0, b.y(), EPS);
    EXPECT_NEAR(3.0, b.z(), EPS);
}


TEST(TlPosition, squareDistanceFrom)
{
    TlPosition a(1.0, 2.0, 3.0);
    TlPosition b(4.0, 5.0, 6.0);
    
    double squareDistance = a.squareDistanceFrom(b);
    EXPECT_NEAR(27.0, squareDistance, EPS);
}


