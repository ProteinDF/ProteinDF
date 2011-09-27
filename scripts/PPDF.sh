#!/bin/sh

PDF_HOME=${HOME}/local
export PDF_HOME

MPI_NUM_PROCS=4
OMP_NUM_THREADS=1
export OMP_NUM_THREADS

cmd="mpiexec -np ${MPI_NUM_PROCS} ${PDF_HOME}/bin/PPDF.x"
eval ${cmd}



