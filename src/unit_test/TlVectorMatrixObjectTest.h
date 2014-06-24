#ifndef TLVECTORMATRIXOBJECTTEST_H
#define TLVECTORMATRIXOBJECTTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlVectorMatrixObject.h"

class TlVectorMatrixObjectTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlVectorMatrixObjectTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testConstructer2);
    CPPUNIT_TEST(testResize);
    CPPUNIT_TEST(testContents);
    CPPUNIT_TEST_SUITE_END();

public:
    void testConstructer();
    void testConstructer2();
    void testResize();
    void testContents();

public:
    TlVectorMatrixObjectTest() {
    }
    
    void setUp() {
    }
    
    void tearDown() {
    }

private:
    static const double threshold;
};

#endif // TLVECTORMATRIXOBJECTTEST_H
