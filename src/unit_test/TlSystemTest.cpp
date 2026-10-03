#include <cstdlib>
#include <string>

#include "TlSystem.h"
#include "gtest/gtest.h"

TEST(TlSystem, getEnv_unset) {
    ::unsetenv("TEST_TLSYSTEM_UNSET_VAR_12345");
    EXPECT_EQ("", TlSystem::getEnv("TEST_TLSYSTEM_UNSET_VAR_12345"));
}

TEST(TlSystem, getEnv_set) {
    ::setenv("TEST_TLSYSTEM_SET_VAR_12345", "hello_world", 1);
    EXPECT_EQ("hello_world", TlSystem::getEnv("TEST_TLSYSTEM_SET_VAR_12345"));
    ::unsetenv("TEST_TLSYSTEM_SET_VAR_12345");
}
