#ifndef TLFILEMATRIXTEST_H
#define TLFILEMATRIXTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlFileMatrix.h"

class TlFileMatrixTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(TlFileMatrixTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testGet);
    CPPUNIT_TEST(testSet);
    CPPUNIT_TEST(testAdd);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void testConstructer();
    void testGet();
    void testSet();
    void testAdd();
    
public:
    TlFileMatrixTest() {
    }
    
    void setUp() {
    }
    
    void tearDown() {
    }
    
private:
    static const double threshold;
};

#endif // TLFILEMATRIXTEST_H

