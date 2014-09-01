#ifndef TLCOLVECTORMATRIXTEST_H
#define TLCOLVECTORMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlColVectorMatrix.h"

class TlColVectorMatrixTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlColVectorMatrixTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testConstructer2);
    CPPUNIT_TEST(testResize);
    CPPUNIT_TEST(testContents);
    CPPUNIT_TEST(testSaveLoad);
    CPPUNIT_TEST(testToTlMatrix);
    CPPUNIT_TEST_SUITE_END();

public:
    void testConstructer();
    void testConstructer2();
    void testResize();
    void testContents();
    void testSaveLoad();
    void testToTlMatrix();

public:
    TlColVectorMatrixTest() {
    }
    
    void setUp() {
    }
    
    void tearDown() {
    }

private:
    static const double threshold;
};

#endif // TLCOLVECTORMATRIXTEST_H
