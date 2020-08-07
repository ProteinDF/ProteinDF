.. -*- coding: utf-8; -*-

************************
Ground State Calculation
************************

Calculation method
==================

ProteinDF allows the user to perform closed-shell, open-shell, 
and open-shell restricted calculations according 
to the electron configuration of the object. 
The calculation method can be specified using the keyword ``method``.

Specifying electron count/occupied state
========================================

restricted Kohn-Sham method
---------------------------

This method performs the calculation of a closed-shell molecule (singlet) 
where an electron pair occupies a single molecular orbital. 
To perform the calculation, specify ``method = rks``.

Further, specify the electron configuration using the keyword ``method/rks/occlevel``.
For example, when ten paired electrons occupy five low energy molecular orbitals, 
specify as follows:

.. code-block:: none
   
   method/rks/occlevel = [1 - 5]

Here the lowest energy orbital is No. 1, and as the number increases, 
a higher energy orbital is assigned.


Specify electron count by using the keyword ``rks/electrons``. 
In the case of ten electrons, specify as follows:

.. code-block:: none
   
   rks/electrons = 10


.. note::
   
   For closed-shell calculations, the ``rks/electrons`` value must be an even number.

.. note::

   The user can also specify the electron configuration as in ``method/rks/occlevel = [1-4, 6]``. 
   In this case, electrons occupy the first to fourth orbitals, and the sixth orbital.


unrestricted Kohn-Sham method
-------------------------------------------

This method performs the calculation of an open-shell molecule, 
where the electron configurations of alpha spin and beta spin differ. 
To perform the calculation, specify ``method=uks``.

Specify electron counts for alpha and beta spins using ``uks/alpha_electrons`` and 
``uks/beta_electrons``, respectively. 
Similarly, specify the electron configuration for alpha and beta spins 
with ``method/uks/alpha-spin-occlevel`` and ``method/uks/beta-spin-occlevel``, respectively. 

.. warning::
   
   Specify a ``uks/alpha_electrons`` value equal to or larger than the ``uks/beta_electrons`` value.


restricted open shell Kohn-Sham method
----------------------------------------------------

This method performs the calculation by categorizing closed- and open-shell electron configurations. 
To perform the calculation, specify ``method=roks``.

Specify electron counts with ``method/roks/closed_shell_electrons`` and 
``method/roks/open_shell_electrons``, 
and electron configuration with ``method/roks/closed_shell_occlevel`` 
and ``method/roks/open_shell_occlevel``.


Initial guess
=============

Use the keyword ``guess`` to specify the initial guess for SFC loops.

core
----

Generates initial guess using core Hamiltonian. 
Specify ``guess=core``

Hückel
------

enerates initial guess with the Hückel method. 
Specify ``guess=huckel``.


Harris functional
-----------------

Generates initial guess with the Harris functional. 
Specify ``guess=harris``. 
This function is not available for some atomic species.


Approximated electron density
-----------------------------

Generates initial guess with an approximated electron density. 
Specify ``guess=rho``. 
A reliable result may not be obtained unless the RI method is applied.


LCAO coefficient matrix
-----------------------

Generates initial guess with the LCAO coefficient matrix. 
Specify ``guess=lcao``. When using this function, 
prepare the LCAO coefficient matrix beforehand.

.. note::
   
   In the current version of the program, 
   it is necessary to prepare LCAO text files and OCC text files in the calculation directory. 
   This specification may change in the future.


Density matrix
--------------

Uses a density matrix as the initial guess. 
Specify ``guess=density``. 
Prepare the density matrix beforehand.

.. note::
   
   In the current version of the program, 
   it is necessary to put the 0th density matrix file in the work directory (fl_Work). 
   This specification may change in the future.


oulomb term calculation
=======================

Selecting calculation engine
----------------------------

Four-center two electron integrals required in Coulomb term calculation is
a rate-determining process. 
Several calculation engines are implemented on ProteinDF for the calculation. 
Use the keyword ``J_engine`` for selection.

conventional
^^^^^^^^^^^^

Calculates four-center two electron integrals at each SCF iteration 
to obtain the Coulomb term.


RI_J
^^^^

Calculates three-center integrals at each SCF iteration based on the RI method 
to obtain the Coulomb term. 
The calculation accuracy depends on auxiliary basis sets. 

Cholesky decomposition
^^^^^^^^^^^^^^^^^^^^^^

Based on the Cholesky decomposition method, 
obtains Cholesky vectors for four-center two-electron integrals before SCF loops. 
The Coulomb term is obtained during each SCF iteration through density matrix operation. 
High-speed computation is allowed since no molecular integral is 
executed during SCF calculations, 
but a large amount of memory and disk is consumed. 
Specify ``J_engine=CD`` to select this engine.
This is the default engine of ProteinDF.


Fock exchange term calculation
==============================

Selecting calculation engine
----------------------------

The Fock exchange term calculation is also rate-determining 
since it requires four-center two electron integrals. 
Use the keyword ``k_engine`` to select a calculation engine.

conventional
^^^^^^^^^^^^

Calculates four-center two electron integrals at each SCF iteration 
to obtain the Fock exchange term. 
This is the default engine of ProteinDF. 
Specify ``K_engine=conventional`` to select this engine.


Cholesky decomposition
^^^^^^^^^^^^^^^^^^^^^^

Obtains the Fock exchange term using the Cholesky decomposition method, 
as in the Coulomb term calculation. 
This engine uses the Cholesky vectors obtained through the Cholesky decomposition 
for the Coulomb term calculation. 
Specify ``K_engine=CD`` to select this engine.


Hybrid functional method and Hartree-Fock method
------------------------------------------------

The user can perform a hybrid functional calculation or Hartree-Fock calculation 
by specifying the following value in the parameter ``xc-potential``:

* HF

  Performs electron state calculations by the Hartree-Fock method.

* B3LYP

  Performs hybrid functional calculations with the Becke 3-parameter.


Exchange-correlation term calculation
=====================================

In ProteinDF, 
the user can use numerical integral calculation or analytical calculation (grid-free method) 
to obtain the exchange-correlation term of the Kohn-Sham matrix, 
as well as the exchange-correlation energy. 
The default is the numerical integrals.


Selecting the grid
------------------

Specify the numerical grid with the parameter ``xc-potential/grid-type``. 
The default is the SG-1 grid. Refer to Appendix for details.


Functionals available for numerical integral method
---------------------------------------------------

Specify functionals with ``xc_potential``. 
The available functionals are as follows:

* SVWN~
* SVWN
* BLYP
* B3LYP
* HFB

.. note::
   
   For the exchange-correlation functional followed by a tilde ``~``, 
   ProteinDF obtains the exchange-correlation term 
   with an approximated electron density based on the RI method.


Grid free method
----------------

Calculates the exchange-correlation term with a grid-free method. 
For details, see the keyword ``grid_free`` in Appendix.


Level shift calculation
=======================

This method allows shifting the energy level of a particular orbital. 
For details, see the keyword ``level_shift`` in Appendix.


Convergence acceleration techniques
===================================

ProteinDF provides several convergence algorisms to achieve a stable 
and efficient convergence during SFC loops. 
Use the keyword ``scf_acceleration`` for selection.


damping method
--------------

The physical quantity used in the last iteration is mixed to the current in a certain ratio. 
When Y (n) represents the physical quantity obtained at the nth SCF iteration, 
the updated amount X (n) can be obtained as follows:

.. math::
   
   X^{\left(n\right)}\leftarrow aX^{\left(n-1\right)}+\left(1-a\right)Y^{\left(n\right)} 

   \left(0<a<1\right)


Here, specify the mixing ratio (a) and target physical quantity 
with ``scf_acceleration/damping/damping_factor`` 
and ``scf_acceleration/damping/damping_type``, respectitvely.


* Example

.. code-block:: none
   
   scf_acceleration/damping/damping_factor = 0.85
   scf_acceleration/damping/damping_type = density_matrix


Anderson's method
-----------------

Employs the quadratic convergence method developed by Anderson. 
The equations when using the physical quantities at the past two points are as follows:

.. math::
   
   X^{\left(n\right)}=u^{\left(n-1\right)}+b^{\left(n-1\right)}\left(v^{\left(n\right)}-u^{\left(n-1\right)}\right)

   u^{\left(n-1\right)}=X^{\left(n-1\right)}+\theta^{\left(n-1\right)}\left(X^{\left(n-2\right)}-X^{\left(n-1\right)}\right)

   v^{\left(n\right)}=Y^{\left(n\right)}+\theta^{\left(n-1\right)}\left(Y^{\left(n-1\right)}-Y^{\left(n\right)}\right)

   \theta^{\left(n-1\right)}=\frac{\left(r^{\left(n-1\right)},r^{\left(n-1\right)}-r^{\left(n-2\right)}\right)}{\left(r^{\left(n-1\right)}-r^{\left(n-2\right)},r^{\left(n-1\right)}-r^{\left(n-2\right)}\right)}

   r^{\left(n-1\right)}=Y^{\left(n\right)}-X^{\left(n-1\right)}

   \left(u,v\right)=\sum_{i}u_{i}v_{i}w_{i}


Here, specify the b(n-1) with ``scf-acceleration/anderson/damping-factor``. 

The damping method is applied before the Andarson's method is started. 
Specify the SCF iteration number starting the Anderson's method 
with ``scf_acceleration/anderson/start_number``. 


DIIS method
-----------

Employs the Direct Inversion of the Iterative Subspace (DIIS) method by Pulay. 
The DIIS method assumes that a new physical quantity X (n) can be obtained 
by the linear combination of X (n-i) in the past.

.. math::
   
   \displaystyle{X^{\left(n\right)}=\sum_{i=i_{0}}^{M}c_{i}X^{\left(n-i\right)}}

   \left(i_{0} \ge 1,\ i_{0}<M \le n-1\right)
   

Here, specify the number of references M with ``scf-acceleration/diis/number-of-diis``.

The damping method is applied before the DIIS method is started. 
Specify the SCF iteration number starting the DIIS method 
with ``scf-acceleration/diis/start-number``.
