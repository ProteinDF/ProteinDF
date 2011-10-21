#!/bin/sh

BASEDIR=..
CP_CMD="cp"

$CP_CMD $BASEDIR/include/*.h .
$CP_CMD $BASEDIR/src/pdflib/*.cpp .
$CP_CMD $BASEDIR/src/pdflib/*.cxx .
$CP_CMD $BASEDIR/src/pdf/*.h .
$CP_CMD $BASEDIR/src/pdf/*.c .
$CP_CMD $BASEDIR/src/pdf/*.cpp .
$CP_CMD $BASEDIR/src/pdf/*.cxx .

OBJ=`ls *.c* | sed "s/\.\(cpp\|cxx\)/\.o/g"`


