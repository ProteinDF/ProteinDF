#!/bin/sh
set -x

rm -rf autom4te.cache

libtoolize -c -f &&
aclocal -I macros &&
autoheader &&
automake --foreign --add-missing --copy &&
autoconf

