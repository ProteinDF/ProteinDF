#include <limits>
#include "TlFileSymmetricMatrixTest.h"
#include "TlFile.h"
#include "TlSymmetricMatrix.h"

const double TlFileSymmetricMatrixTest::threshold = std::numeric_limits<double>::epsilon();

void TlFileSymmetricMatrixTest::testConstructer(){
  TlSymmetricMatrix sm(100);
  sm( 3, 17) =  51.0;
  sm( 0, 28) = - 1.0;
  sm.save("sym_matrix.bin");

  TlFileSymmetricMatrix fsm("sym_matrix.bin");
}

void TlFileSymmetricMatrixTest::testGet() {
  TlSymmetricMatrix sm(100);
  sm( 3, 17) =  51.0;
  sm( 0, 28) = - 1.0;
  sm.save("sym_matrix.bin");

  TlFileSymmetricMatrix fsm("sym_matrix.bin");

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fsm.get( 0,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fsm.get( 0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fsm.get(99,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fsm.get(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(51.0, fsm.get( 3, 17), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, fsm.get( 0, 28), threshold);
}

void TlFileSymmetricMatrixTest::testSet() {
  const std::string filename = "file_sym_matrix.bin";
  TlFile::remove(filename);
  {
    TlFileSymmetricMatrix fsm(filename, 100);
    
    fsm.set( 3, 17, 51.0);
    fsm.set( 0, 28, -1.0);
  }

  TlSymmetricMatrix sm;
  sm.load(filename);

  CPPUNIT_ASSERT_EQUAL(100, sm.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(100, sm.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm( 0,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm( 0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm(99,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(51.0, sm( 3, 17), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, sm( 0, 28), threshold);

}

void TlFileSymmetricMatrixTest::testAdd() {
  {
    TlSymmetricMatrix sm(100);
    sm( 3, 17) =  51.0;
    sm( 0, 28) = - 1.0;
    sm.save("sym_matrix.bin");
  }

  {
    TlFileSymmetricMatrix fsm("sym_matrix.bin");
    fsm.add(3, 17, 12.3);
  }

  TlSymmetricMatrix sm;
  sm.load("sym_matrix.bin");
  CPPUNIT_ASSERT_EQUAL(100, sm.getNumOfRows());
  CPPUNIT_ASSERT_EQUAL(100, sm.getNumOfCols());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm( 0,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm( 0, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm(99,  0), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, sm(99, 99), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(63.3, sm( 3, 17), threshold);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, sm( 0, 28), threshold);

}


