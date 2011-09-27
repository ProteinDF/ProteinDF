#include <vector>
#include <string>
#include "TlStringTokenizerTest.h"

// CPPUNIT_ASSERT( condition );
// condition$B$,56(B(false,0)$B$G$"$C$?$H$-!"<:GT$7$^$9!#(B
//
// CPPUNIT_ASSERT_MESSAGE( message, condition );
// condition$B$,56$G$"$C$?$H$-!"<:GT$7$^$9!#$3$N$H$-(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_FAIL( message );
// $BI,$:<:GT$7$^$9!#(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_ASSERT_EQUAL( expected, actual );
// $BF@$i$l$?7k2L(Bactual$B$,4|BT$9$kCM(Bexpected$B$G$J$+$C$?$H$-!"$9$J$o$A(Bexpected != actual$B$N$H$-$K<:GT$7$^$9!#(B
//
// CPPUNIT_ASSERT_EQUAL_MESSAGE( message, expected, actual );
// $BF@$i$l$?7k2L(Bactual$B$,4|BT$9$kCM(Bexpected$B$G$J$+$C$?$H$-!"$9$J$o$A(Bexpected != actual$B$N$H$-$K<:GT$7$^$9!#$3$N$H$-(Bmessage$B$r=PNO$7$^$9!#(B
//
// CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, delta );
// $BF@$i$l$?7k2L(Bactual$B$H4|BT$9$kCM(Bexpected$B$H$N:9$,(Bdelta$B$h$jBg$-$$$H$-!"<:GT$7$^$9!#(B

void TlStringTokenizerTest::testConstructer(){
  TlStringTokenizer st("This is a test string.");
}

void TlStringTokenizerTest::testCountTokens(){
  TlStringTokenizer st("This is a test string.");

  int count = st.countTokens();

  CPPUNIT_ASSERT_EQUAL(5, count);
}

void TlStringTokenizerTest::testGetTokens(){
  TlStringTokenizer st("This is a test string.");
  
  std::vector<std::string> tokens;
  while (st.hasMoreTokens() == true){
    std::string s = st.nextToken();
    tokens.push_back(s);
  }

  CPPUNIT_ASSERT_EQUAL(size_t(5), tokens.size());
  CPPUNIT_ASSERT_EQUAL(std::string("This"),    tokens[0]);
  CPPUNIT_ASSERT_EQUAL(std::string("is"),      tokens[1]);
  CPPUNIT_ASSERT_EQUAL(std::string("a"),       tokens[2]);
  CPPUNIT_ASSERT_EQUAL(std::string("test"),    tokens[3]);
  CPPUNIT_ASSERT_EQUAL(std::string("string."), tokens[4]);
}


