#include <vector>
#include <string>
#include "TlStlUtilsTest.h"

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

void TlStlUtilsTest::testCache() {
    TlCache<int, std::string> cache(10);
    CPPUNIT_ASSERT_EQUAL(std::size_t(10), cache.getMaxItems());

    cache.setMaxItems(8);
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getMaxItems());

    cache.set(1, "one");
    cache.set(2, "two");
    cache.set(3, "three");
    cache.set(4, "four");
    cache.set(5, "five");

    CPPUNIT_ASSERT_EQUAL(std::size_t(5), cache.getNumOfItems());
    CPPUNIT_ASSERT_EQUAL(std::string("three"), *(cache.get(3)));

    cache.set(6, "six");
    cache.set(7, "seven");
    cache.set(8, "eight");
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getNumOfItems());

    cache.set(9, "nine");
    CPPUNIT_ASSERT_EQUAL(std::size_t(8), cache.getNumOfItems());

    cache.set(10, "ten");
    CPPUNIT_ASSERT_EQUAL(std::string("five"), *(cache.get(5)));

    // empty
    CPPUNIT_ASSERT_EQUAL(std::string(""), *(cache.get(1)));
}

