#include <limits>
#include "TlFileMatrixTest.h"
#include "TlFile.h"
#include "TlMatrix.h"

const double TlFileMatrixTest::threshold = std::numeric_limits<double>::epsilon();

void TlFileMatrixTest::testConstructer()
{
    TlMatrix m(100, 200);
    m( 3, 17) =  51.0;
    m( 0, 28) = - 1.0;
    m.save("matrix.bin");

    TlFileMatrix fm("matrix.bin");
}

void TlFileMatrixTest::testGet() {
    TlMatrix m(100, 200);
    m( 3, 17) =  51.0;
    m( 0, 28) = - 1.0;
    m.save("matrix.bin");

    TlFileMatrix fm("matrix.bin");

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fm.get( 0,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fm.get( 0, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fm.get(99,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fm.get(99, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(51.0, fm.get( 3, 17), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, fm.get( 0, 28), threshold);
}

void TlFileMatrixTest::testSet() {
    const std::string filename = "file_matrix.bin";
    TlFile::remove(filename);
    {
        TlFileMatrix fm(filename, 100, 200);
        
        fm.set( 3, 17, 51.0);
        fm.set( 0, 28, -1.0);
    }
    
    TlMatrix m;
    m.load(filename);
    
    CPPUNIT_ASSERT_EQUAL(100, m.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(200, m.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m( 0,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m( 0, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m(99,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m(99, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(51.0, m( 3, 17), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, m( 0, 28), threshold);
}

void TlFileMatrixTest::testAdd() {
    {
        TlMatrix m(100, 200);
        m( 3, 17) =  51.0;
        m( 0, 28) = - 1.0;
        m.save("matrix.bin");
    }

    {
        TlFileMatrix fsm("matrix.bin");
        fsm.add(3, 17, 12.3);
    }
    
    TlMatrix m;
    m.load("matrix.bin");
    CPPUNIT_ASSERT_EQUAL(100, m.getNumOfRows());
    CPPUNIT_ASSERT_EQUAL(200, m.getNumOfCols());
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m( 0,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m( 0, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m(99,  0), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, m(99, 99), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(63.3, m( 3, 17), threshold);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, m( 0, 28), threshold);
}


