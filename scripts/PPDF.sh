#!/bin/sh

cmd="mpiexec -np ${MPI_NUM_PROCS} ${PDF_HOME}/bin/PPDF.x"
eval ${cmd}



