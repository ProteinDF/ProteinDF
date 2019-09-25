#include "digit.h"
#include <cstdio>

int digit(int c) {
    char dummy[256];
    int a = sprintf(dummy, "%d", c);

    return a;
}
