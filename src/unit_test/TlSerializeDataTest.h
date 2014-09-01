#ifndef TLSERIALIZEOBJECTTEST_H
#define TLSERIALIZEOBJECTTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "TlSerializeData.h"

class TlSerializeDataTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TlSerializeDataTest);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testCopyConstructor);
    CPPUNIT_TEST(testMapAccess);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void testConstructor();
    void testCopyConstructor();
    void testMapAccess();
    
public:
    TlSerializeDataTest(){
    }

    void setUp(){
    }
    
    void tearDown(){
    }
};

#endif // TLSERIALIZEOBJECTTEST_H

