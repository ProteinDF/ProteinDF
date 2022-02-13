===================
Execution Procedure
===================

Program execution procedure
===========================

Follow the procedure below to execute ProteinDF from the command line:


Preparation
-----------

Create directories (default: ``fl_Work``) for outputting intermediate data,
under the execution directory of ProteinDF (where the input file is located).

.. note::

  These directories can be created with the ``pdf-setup`` command.

.. note::

  The data written in these directories will be extremely large.
  It is recommended to create the directories in a high-speed disk storage
  with large capacity.


Executing the program (serial version)
---------------------------------------

To execute the serial version of ProteinDF, use the following command:

.. code-block:: bash

  % $PDF_HOME/bin/PDF.x

When the computation starts,
total energy at each SCF calculation is sequentially displayed
in the standard output.
In addition, the series of the calculation result data is output
in the log file (``fl_Out_Std``).
Intermediate data during all-electron calculation is also output
in the log file.


Executing the program (parallel version)
----------------------------------------

To execute the parallel version of ProteinDF, use the following command:

.. code-block:: bash

  % mpiexec -n N $PDF_HOME/bin/PPDF.x

Here, specify ``N``, the number of processors for parallel computation.

.. note::

  Execution procedure of the MPI program varies depending
  on the computing system environment.
  For details, refer to the system manuals.

When the computation starts,
total energy at each SCF calculation is sequentially displayed
in the standard output.
In addition, the series of the calculation result data is
output in a text file (default file name: ``fl_Out_Std``), as in the serial mode.


========
Run Type
========

Overview
========

ProteinDF has several run types to efficiently compute
a large object with limited computing resources.

======================== ================= ============
run type                 parallel method   matrix
======================== ================= ============
serial                   OpenMP only       replica
------------------------ ----------------- ------------
replica_static           MPI/OpenMP hybrid replica
------------------------ ----------------- ------------
replica_dynamic          MPI/OpenMP hybrid replica
------------------------ ----------------- ------------
distributed              MPI/OpenMP hybrid distributed
======================== ================= ============


serial
======

* Calculations are performed by only one process. No inter-process communication is performed.
* Thread parallel computing by OpenMP is possible.
* LAPACK is used for matrix operations.
* The size of the computable system depends on the amount of memory available to the process.


replica_static
==============

* Performs parallel computation using inter-process communication (MPI). Within each process, OpenMP parallel computation is performed.
* The process replicates and maintains all elements of the matrix required for the calculation in each MPI process.
* If the matrix cannot be allocated with the specified memory size, the matrix is kept in the disk area.
* Divide and conquer method is used to distribute the tasks.
* Use LAPACK for matrix operations.
* The ``linear_algebra_package`` keyword is ``LAPACK``.

.. note::

  The divide-and-conquer method is effective when the number of processes is small, since all processes participate in the arithmetic operations. However It has the disadvantage that the load cannot be distributed evenly.


replica_dynamic
===============

* Performs parallel computation using inter-process communication (MPI). Within each process, OpenMP parallel computation is performed.
* The process replicates and maintains all elements of the matrix required for the calculation in each MPI process.
* If the matrix cannot be allocated with the specified memory size, the matrix will be kept in the disk area.
* The leader/follower method is used for task distribution.
* Use LAPACK for matrix operations.
* ``linear_algebra_package`` keyword is ``LAPACK``.
* The ``parallel_processing_type`` keyword must be ``MS``.


.. note::

  In the leader/follower method, the master process concentrates on distributing tasks and does not perform operations.
  Since the load can be distributed evenly, it is effective when the number of processes is large.



distributed
===========

* It performs parallel computation by inter-process communication (MPI). OpenMP parallel computing is performed within each process.
* The global matrix is divided and maintained in each MPI process.
* If the matrix cannot be allocated in the specified memory size, the matrix is held in the disk area.
* ScaLAPACK is used for matrix operations.
* Specify ``ScaLAPACK`` for the ``linear_algebra_package`` keyword.
