#include <vector>
#include <string>
#include "gtest/gtest.h"
#include "TlStlUtils.h"

TEST(TlStlUtils, cache)
{
    TlCache<int, std::string> cache(10);
    EXPECT_EQ(std::size_t(10), cache.getMaxItems());

    cache.setMaxItems(8);
    EXPECT_EQ(std::size_t(8), cache.getMaxItems());

    cache.set(1, "one");
    cache.set(2, "two");
    cache.set(3, "three");
    cache.set(4, "four");
    cache.set(5, "five");

    EXPECT_EQ(std::size_t(5), cache.getNumOfItems());
    EXPECT_EQ(std::string("three"), *(cache.get(3)));

    cache.set(6, "six");
    cache.set(7, "seven");
    cache.set(8, "eight");
    EXPECT_EQ(std::size_t(8), cache.getNumOfItems());

    cache.set(9, "nine");
    EXPECT_EQ(std::size_t(8), cache.getNumOfItems());

    cache.set(10, "ten");
    EXPECT_EQ(std::string("five"), *(cache.get(5)));

    // empty
    EXPECT_EQ(std::string(""), *(cache.get(1)));
}

