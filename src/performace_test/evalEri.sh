#!/bin/sh


CMD="./evalEri"
PDF_PARAM="./1FUL.param"
P_MATRIX="./1FUL_P.mtx"
K_MATRIX="./1FUL_K.check.mtx"

/usr/bin/time ${CMD} -v -d ${P_MATRIX} -p ${PDF_PARAM} -s ${K_MATRIX} 

echo 
echo ">>>> opreport"
opreport ${CMD} 

echo
echo ">>>> opreport -c"
opreport -c ${CMD} 

echo
echo ">>>> opreport -dg"
opreport -dg ${CMD} 

echo 
echo ">>>> opreport -l"
opreport -l ${CMD} 

echo
echo ">>>> opannotate -s"
opannotate -s ${CMD} 


