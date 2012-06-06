#!/bin/sh

BASEDIR=../
CP="ln -sf"

$CP $BASEDIR/include/*.h .
$CP $BASEDIR/src/pdflib/*.cpp .
$CP $BASEDIR/src/pdflib/*.cxx .
$CP $BASEDIR/src/pdf/*.h .
$CP $BASEDIR/src/pdf/*.c .
$CP $BASEDIR/src/pdf/*.cpp .
$CP $BASEDIR/src/pdf/*.cxx .
$CP $BASEDIR/src/tools/*.cpp .
$CP $BASEDIR/src/tools/*.cxx .

$CP $BASEDIR/src/performace_test/*.cpp .
$CP $BASEDIR/src/performace_test/*.h .

$CP fj_config.h config.h
$CP fj_pdflib-int.h pdflib-int.h

#OBJ=`ls *.c* | sed "s/\.\(cpp\|cxx\)/\.o/g"`


