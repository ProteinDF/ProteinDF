#!/bin/sh
#set -x

rm -rf autom4te.cache

OS=`uname -s`

if [ $OS = "Darwin" ]; then
    LIBTOOLIZE=glibtoolize
else
    LIBTOOLIZE=libtoolize
fi

${LIBTOOLIZE} -c -f &&
aclocal -I macros &&
autoheader &&
automake --foreign --add-missing --copy &&
autoconf

