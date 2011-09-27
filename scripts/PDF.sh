#!/bin/sh

PDF_HOME=${HOME}/local
export PDF_HOME

OMP_NUM_THREADS=1
export OMP_NUM_THREADS


cmd=${PDF_HOME}/bin/PDF.x
eval ${cmd}


