
**************
 Input files
**************

ProteinDF reads calculation parameters from ASCII text input files. 
The user can use any preferred editor to create and edit the input files.

.. index:: fl_Userinput

The default input file name is ``fl_Userinput``. 
If a file named ``fl_Userinput`` exists in the current directory, 
the program reads the file as an input.

The input file for ProteinDF consists of the following sections:

* MAIN

* MOLECULE

.. note::

  Although the previous versions of ProteinDF categorized the keywords in groups, 
  the current version does not require such grouping. 
  Specify keywords at any arbitrary location."

.. warning::
   
   The input files only accept the default locale (LANG=C). 
   If a Japanese Kanji-code (e.g. UTF-8, EUC, Shift-JIS) is used in input files, 
   the program may not obtain an accurate result. 
   In particular, pay attention to the use of a blank (space) character.

.. warning::
   
   In input files, be sure to use the line feed code in accordance with the user's system. 
   Otherwise, the program may not obtain an accurate result. 
   Most of the UNIX systems use the LF (0x0A) code. 
   Be sure to take an appropriate measure when transferring the files 
   created in a Windows system (line feed code: CR(0x0D)+LF) 
   to UNIX systems via FTP or SFTP.


Syntax
======

Specify the keywords in the following format:

.. code-block:: none
                  
   keyword = value

The keywords are case-insensitive.

The value can be of the following three types:

* A single value

A value which does not contain a blank (space, tab, line feed, etc.) 
can be directly specified as follows:

.. code-block:: none
   
   max-iteration = 100


* A single-line text

A text which does not contain a line feed can be specified 
between brackets (``[ ]``), as follows:

.. code-block:: none
   
   method/nsp/occlevel = [ 1 - 20 ]


* Multiple-line text

Specify a text which contains a line feed between braces (``{``, ``}end``) 
as follows:

.. code-block:: none
   
   geometry/cartesian/input = {
                // molecule1
                N        -1.888000        0.035000       -0.211000
                H        -1.766000        0.945000        0.189000
                H        -1.817000        0.099000       -1.205000
                C        -0.758000       -0.730000        0.287000
                H        -0.893000       -0.915000        1.372000
                H        -0.720000       -1.725000       -0.200000
                C         0.529000        0.065000        0.064000
                O         0.520000        1.294000        0.114000
                O         1.742000       -0.451000       -0.186000
                H         1.692000       -1.400000       -0.203000
   }

.. note::
   
   Any text following double-slash ``//`` to the end of the line is 
   treated as a comment. 
   Any line beginning with a hash ``#`` is also ignored as a comment.

.. note::

  If identical keywords are specified, the latter one overwrites the earlier.


MAIN section
==============

Specify the parameters necessary for calculation. 
For details, see the parameter list in Appendix. 
The following describes the most frequently used parameters:

.. index:: step_control

step_control
^^^^^^^^^^^^

Specifies calculation schemes. 
The ProteinDF performs calculations according to the scheme order specified here.

Values
""""""""""

* create

Analyzes the parameters. No calculation is performed.

* integral

Executes the pre-processing for SCF loops, such as core Hamiltonian, overlap integrals, and grid generation. 

* guess

Generates initial guess.

* SCF

Executes SCF loops.

* force

Calculates derivatives of energy with respect to nuclear coordinates.


Example
"""""""

.. code-block:: none
                
   step_control = [create integral guess scf]


.. index:: scf_start_guess

scf_start_guess
^^^^^^^^^^^^^^^

Specifies initial guess for SCF calculations.

Values
""""""""""

* huckel

  Obtains initial guess with the Hückel method.

* harris

  Obtains initial guess using the Harris functional from the electron density of atoms previously prepared.

* core

  Obtains initial guess from the wave function obtained by diagonalizing the core Hamiltonian.

* rho

  Merges the approximated electron density of each atom previously prepared, and generates an approximated electron density of the model molecule.

* file_rho

  Obtains an approximated electron density using an auxiliary basis expansion coefficient file ``guess_rho`` created by the user, and uses the value as initial guess.

* lcao

  Generates initial guess from the user-created LCAO matrix file (``guess.lcao.rks``) and occupation number file (``guess.occ.rks``).

* density_matrix

  Uses the user-created electron density file as initial guess.

MOLECULE section
==================

Specify the following keywords:

.. index:: geometry/cartesian/unit

geometry/cartesian/unit
^^^^^^^^^^^^^^^^^^^^^^^

Specifies the unit of length used for input coordinates.

Value
""""""""""

angstrom, au


Exaple
""""""

.. code-block:: none
                
   geometry/cartesian/unit = angstrom


.. index:: geometry/cartesian/input

geometry/cartesian/input
^^^^^^^^^^^^^^^^^^^^^^^^

Specifies the target atomic species to be calculated and their XYZ coordinates. 
In each line, specify atomic species, X, Y, and Y coordinates sequentially, 
separating each value with a blank (space or tab) character. 
Specify atomic species with the atomic symbols.


Example
"""""""

.. code-block:: none
   
   geometry/cartesian/input        = {
                N        -1.888000        0.035000       -0.211000
                H        -1.766000        0.945000        0.189000
                H        -1.817000        0.099000       -1.205000
                C        -0.758000       -0.730000        0.287000
                H        -0.893000       -0.915000        1.372000
                H        -0.720000       -1.725000       -0.200000
                C         0.529000        0.065000        0.064000
                O         0.520000        1.294000        0.114000
                O         1.742000       -0.451000       -0.186000
                H         1.692000       -1.400000       -0.203000
   }end

.. note::
   
   The user can add a label by attaching ``@`` after atomic symbols. 
   This function is useful when assigning basis sets to the same element 
   in separate groups.

.. note::
   
   A dummy atom can be specified with ``X``. 
   In that case, specify the electric charge of the dummy atom 
   in the fifth column.


.. index:: basis-set/orbital

basis-set/orbital
^^^^^^^^^^^^^^^^^

Specifies basis sets to all atomic species used for calculation. 
Describe the name of the basis set for each atomic species. 
The names of the assigned basis sets must be previously specified 
in the basis2 file. See Appendix for the basis2 file.


Example
"""""""

.. code-block:: none
   
   basis-set/orbital = {
                H = "O-HYDROGEN (41) DZVP"
                O = "O-OXYGEN (621/41) by FS"
                C = "O-CARBON (621/41) by FS"
                N = "O-NITROGEN (621/41) by FS"
   }end


.. note::
   
   The user can add a label by attaching ``@`` after atomic symbols.


.. index:: basis-set/density-auxiliary

basis-set/density-auxiliary
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specifies the auxiliary basis sets used for Coulomb term calculation. 
Use this keyword when calculating the term in the RI_J method. 
The specification procedure is the same as that for basis sets.


.. index:: basis-set/exchange-auxiliary

basis-set/exchange-auxiliary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specifies the auxiliary basis sets used for exchange-correlation term calculation. 
Use this keyword when calculating the term in the RI method 
(i.e. when the user attached ``~`` to the end of the specified exchange-correlation functional.)
The specification procedure is the same as that for basis sets.


basis-set/gridfree
^^^^^^^^^^^^^^^^^^

Specifies the auxiliary basis sets used for grid-free method.