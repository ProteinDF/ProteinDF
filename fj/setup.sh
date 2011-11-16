#!/bin/sh

BASEDIR=..
CP="cp -f"

$CP $BASEDIR/include/*.h .
$CP $BASEDIR/src/pdflib/*.cpp .
$CP $BASEDIR/src/pdflib/*.cxx .
$CP $BASEDIR/src/pdf/*.h .
$CP $BASEDIR/src/pdf/*.c .
$CP $BASEDIR/src/pdf/*.cpp .
$CP $BASEDIR/src/pdf/*.cxx .

$CP $BASEDIR/src/performace_test/*.cpp .
$CP $BASEDIR/src/performace_test/*.h .

OBJ=`ls *.c* | sed "s/\.\(cpp\|cxx\)/\.o/g"`


