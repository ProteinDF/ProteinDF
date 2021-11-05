*****************
Quick Start Guide
*****************

This chapter describes how to perform a single-point calculation of amino acid (glycine) using a single process.

.. index:: fl_Work

Creating calculation directories
================================

ProteinDF writes a large number of large intermediate files with matrices/vectors on disk storage. 
Before starting actual calculations, create the following calculation directories:

* fl_Work


.. note::

  Exec ``${PDF_HOME}/bin/pdf-setup`` command, and prepare the environment for the calculation.


.. warning::

  The program may abnormally terminate if the intermediate files cannot be written properly. 
  Before executing MPI parallel computation, make sure that the files can be properly written from all nodes.


Create input file
==================

.. index:: fl_Userinput

Create a text file with the following content, and save it under the name ``fl_Userinput``:


.. literalinclude:: fl_Userinput


.. note::

  Input files can be changed with a startup option of ProteinDF.


Exec ProteinDF
===============

Configure the environment variable PDF_HOME appropriately, and execute ProteinDF (in serial mode).

.. code-block:: bash

  % ${PDF_HOME}/bin/PDF.x


If the program terminates properly, the system returns to the command prompt.


Results
==========

The calculation results are output in files. The following shows an example of the output:

.. note::

  The output file location can be changed with a startup option of ProteinDF.


The beginning of the output file shows the version of ProteinDF and the number of parallel processes 
(MPI process count, the number of OpenMP threads). 
Make sure that the calculation was performed as intended.

.. code-block:: none
   
   [0:2012/**/07 17:17:02:INFO] **************************************
   [0:2012/**/07 17:17:02:INFO] ProteinDF version 20xx.x:xxxx (serial)
   [0:2012/**/07 17:17:02:INFO] 
   [0:2012/**/07 17:17:02:INFO]  OpenMP threads: 12
   [0:2012/**/07 17:17:02:INFO] 

The calculation is performed according to the procedure described in ``step_control``. 
The output date is indicated at left in the log.

.. code-block:: none
   
   ===============================================
    >>>> INTEGRAL
   ===============================================
    >>>> Hpq
   
    ...
   
   ===============================================
    >>>> GUESS
   ===============================================
   
    ...
   
   ===============================================
    >>>> SCF
   ===============================================


Energy data are output in the ``Total Energy`` block.

.. code-block:: none
   
   ------------------------------------------------
    >>>> Total Energy
   ------------------------------------------------
    Ts+Vn          =        -745.7264230071891689
    E_J[Rho, Rho~] =         629.4907739256434525
    E_J[Rho~,Rho~] =        -311.1519936178789294
    E_xc(pure)     =         -27.8222685895087842
    E_K            =          -6.8830892822011522
    E_nuclei       =         179.8930288412503558
    TE             =        -282.1999717298841688
   ------------------------------------------------

Other information, such as calculation size or cut off data, is output as needed.

When the calculation completes properly, the CPU time and elapsed time are output.

.. code-block:: none
   
   ************************************************
    ProteinDF successful completion
    CPU_TIME:         3454 sec
    ELAPS_TIME:        542 sec
   ************************************************
