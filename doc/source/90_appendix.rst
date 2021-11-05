********
Appendix
********

startup options of ProteinDF
============================

* -r

Resume (restart) the calculation from where it can be restarted.
The calculation conditions are read from the parameter file (pdfparam.mpac).


* -i input_file

Specify an input file to start the calculation.
If not specified, fl_Userinput on the current directory will be used.


* -o output_file

Specifies the output file for the result.
If not specified, fl_Out_Std will be output to the current directory.


* -d

Debug output.
Note that a large number of messages are output to the file and standard output.
Normally, you do not need to select this option.


Input parameter list
====================

.. include:: pdfkwd.rst


Input file format
========================


basis2 file format
==========================

The basis2 file stores the basis function information.
The following is a concrete example.

.. literalinclude:: basis2_sample.txt



* Format

** line 1(Title)

Describes the name of the basis function.
A name starting with "O-" indicates a basis function, and a name starting with "A-" indicates an auxiliary basis function.
In the case of basis functions, a single CGTO block.
In the case of auxiliary basis functions, it consists of two CGTO blocks, one for the Coulomb term and the other for the auxiliary basis functions for the exchange-correlation term.


** line 2(the number of CGTOs)

Specifies the form and number of orbits of a contracted Gaussian basis function (CGTO).
It is represented by three integer values separated by whitespace.
The three integer values, separated by spaces, represent the number of s-type, p-type, and d-type orbits, in that order.
After this line, the CGTO block continues.

If it is written as "3 2 1,"
There are three contracted Gaussian basis functions of type s, two contracted Gaussian basis functions of type p, and one contracted Gaussian basis function of type d.
Therefore, after this line, six CGTO blocks will be described.


** line 3(beginning of CGTO block)

The beginning of a CGTO block indicates the number of primitive Gaussian basis functions (PGTO) that the CGTO block contains.


** line 4(PGTO)

The exponents and coefficients of the primitive basis functions (PGTO) are described, separated by whitespace.
For the auxiliary basis functions, the coefficient (1.0) can be omitted.

