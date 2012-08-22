#ifndef TLSPARSEVECTORTEST_H
#define TLSPARSEVECTORTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "TlSparseVector.h"

class TlSparseVectorTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlSparseVectorTest);
    CPPUNIT_TEST(testConstructer);
    CPPUNIT_TEST(testAdd_double);
    CPPUNIT_TEST(testAdd_TlSparseVector);
    CPPUNIT_TEST_SUITE_END();

public:
    void testConstructer();
    void testAdd_double();
    void testAdd_TlSparseVector();

public:
    TlSparseVectorTest() {
    }

    void setUp() {
    }

    void tearDown() {
    }

private:
    static const double threshold;
};

#endif // TLSPARSEVECTORTEST_H
