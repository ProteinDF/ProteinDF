#!/bin/sh

echo ">>>> start pdf serial"
echo "PDF_HOME is $PDF_HOME"
PDF_PWD=`pwd`
echo "pwd = $PDF_PWD"

cmd=${PDF_HOME}/bin/PDF.x
eval ${cmd}
echo "<<<< end pdf serial"


