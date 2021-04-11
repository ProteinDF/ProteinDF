===================
Execution Procedure
===================

Program execution procedure
===========================

Follow the procedure below to execute ProteinDF from the command line:


Preparation
-----------

Create directories (default: fl_Work) for outputting intermediate data, 
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
in the log file (fl_Out_Std). 
Intermediate data during all-electron calculation is also output 
in the log file.


Executing the program (parallel version)
----------------------------------------

To execute the parallel version of ProteinDF, use the following command:

.. code-block:: bash

  % mpiexec -n N $PDF_HOME/bin/PDF.x

Here, specify ``N``, the number of processors for parallel computation.

.. note::

  Execution procedure of the MPI program varies depending 
  on the computing system environment. 
  For details, refer to the system manuals.

When the computation starts, 
total energy at each SCF calculation is sequentially displayed 
in the standard output. 
In addition, the series of the calculation result data is 
output in a text file (default file name: fl_Out_Std), as in the serial mode.


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

* 1プロセスのみで計算を行います。プロセス間通信は行いません。
* OpenMPによるスレッド並列計算が可能です。
* 行列演算はLAPACKを利用します。
* 計算可能な系のサイズは、プロセスが使用できるメモリ容量に依存します。


replica_static
==============

* プロセス間通信(MPI)による並列計算を行います。各プロセス内ではOpenMP並列計算を行います。
* プロセスは計算に必要な行列の全要素を各MPIプロセスで複製、保持します。
* 指定されたメモリサイズでは行列が確保できない場合は、ディスク領域に行列を保持します。
* タスクの分散は分割統治法を採用します。
* 行列演算はLAPACKを利用します。
* ``linear_algebra_package`` キーワードに ``LAPACK`` を指定します。

.. note::

  分割統治法では、すべてのプロセスが演算処理に参加するので、
  プロセス数が少ない場合に有効です。
  負荷が均等に分散できない欠点があります。

.. note::

  memory_size キーワードでプロセスが利用できるメモリ量を指定します。

.. warning::

  メモリ不足によってディスクを利用した場合は、パフォーマンスが劣化する場合があります。



replica_dynamic
===============

* プロセス間通信(MPI)による並列計算を行います。各プロセス内ではOpenMP並列計算を行います。
* プロセスは計算に必要な行列の全要素を各MPIプロセスで複製、保持します。
* 指定されたメモリサイズでは行列が確保できない場合は、ディスク領域に行列を保持します。
* タスクの分散はマスター-スレーブ法を採用します。
* 行列演算はLAPACKを利用します。
* ``linear_algebra_package`` キーワードに ``LAPACK`` を指定します。
* ``parallel_processing_type`` キーワードに ``MS`` を指定します。

.. note::

  マスター-スレーブ法では、マスタープロセスはタスクの分散に専念し、演算を行いません。
  負荷が均等に分散できるので、プロセス数が多い場合に有効です。


.. note::

  memory_size キーワードでプロセスが利用できるメモリ量を指定します。

.. warning::

  メモリ不足によってディスクを利用した場合は、パフォーマンスが劣化する場合があります。



distributed
===========

* プロセス間通信(MPI)による並列計算を行います。各プロセス内ではOpenMP並列計算を行います。
* 大域行列を各MPIプロセスに分割して保持します。
* 指定されたメモリサイズでは行列が確保できない場合は、ディスク領域に行列を保持します。
* 行列演算はScaLAPACKを利用します。
* ``linear_algebra_package`` キーワードに ``ScaLAPACK`` を指定します。
