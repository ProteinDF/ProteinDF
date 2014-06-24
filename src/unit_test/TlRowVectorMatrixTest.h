#ifndef TLROWVECTORMATRIXTEST_H
#define TLROWVECTORMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlRowVectorMatrix.h"

class TlRowVectorMatrixTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlRowVectorMatrixTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testConstructer2);
    CPPUNIT_TEST(testResize);
    CPPUNIT_TEST(testContents);
    CPPUNIT_TEST(testSaveLoad);
    // CPPUNIT_TEST(testToTlMatrix);
    CPPUNIT_TEST_SUITE_END();

public:
    void testConstructer();
    void testConstructer2();
    void testResize();
    void testContents();
    void testSaveLoad();
    void testToTlMatrix();

public:
    TlRowVectorMatrixTest() {
    }
    
    void setUp() {
    }
    
    void tearDown() {
    }

private:
    static const double threshold;
};

#endif // TLROWVECTORMATRIXTEST_H
