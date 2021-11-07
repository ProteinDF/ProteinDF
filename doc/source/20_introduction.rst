.. -*- coding: utf-8; -*-

********
Overview
********

ProteinDF is an application program focusing on precisely executing/analyzing all-electron canonical molecular orbital (CMO) calculations of proteins.


Features
========

The ProteinDF features are as follows:

* Density functional theory (DFT) program using Gaussian-type basis sets

  * Allows Hartree-Fock (HF) calculations, pure DFT calculations, and hybrid DFT calculations.

* Implemented using the object-oriented program language C++.

* Ground-state all-electron CMO calculation for large molecules including proteins.

  * High-speed molecular integrals with the Resolution of Identity (RI) method.
  * Precise and fast computation with Cholesky decomposition.
  * Large-scale computation by distributing large-scale matrices.
  * Hybrid parallel computation with MPI/OpenMP.

* QM Analysis of Large molecule

  * Analytic calculation of energy gradient at each atomic coordinate.
  * Mulliken population analysis.
  * Calculation of MOs, electron density, and electrostatic potential.


Calculation method
==================

To handle the entire structure of a protein and to predict its functions quantitatively,
it is appropriate to use a DFT method which solves Kohn-Sham equations,
where electron correlation can be incorporated effectively.
The ProteinDF is a Gaussian-basis DFT calculation program based on the Kohn-Sham-Roothaan (KSR) method,
a standard MO calculation method in chemistry.
For details on DFT calculations and large-scale computation,
please refer to the other materials.
This section briefly describes the ProteinDF overview.


To carry out large-scale canonical calculation,
the ProteinDF employs RI technique and Cholesky decomposition to speed up the calculation.
Restricted, Unrestricted, and Restricted Open-Shell Kohn-Sham (KS) methods are available.
The hybrid functional B3LYP is set to the default exchange-correlation potential,
and a local density functional SVWN and generalized gradient approximation (GGA) BLYP are also available.
The SVWN can also be used for calculation of exchange-correlation potential
and also for RI calculations.
Energy gradient calculations can be achieved by these methods.
To use the obtained energy gradient for geometry optimization,
quantum molecular dynamics calculation, and physical quantity calculation,
the user is required to work ProteinDF with the relevant external programs.
The basis sets and auxiliary basis sets are prepared in the basis2 file in the form optimized for Gaussian basis DFT calculations.
The user can also add desired basis sets to the file.
The default set for ProteinDF is Valence Double (VD) equivalent.
For details, refer to the content of the basis2 file.
The user can operate ProteinDF by specifying keywords as in the other MO calculation programs.
Default values are already set to most keywords, some of which were mentioned above.
See Appendix for details.


Computation size
================

The most distinctive feature of ProteinDF is its computation size.
The user can execute all-electron calculation of proteins with several dozens of residues in the current PC system without any difficulties. Among all of our efforts to make the program possible, the biggest contributing factor was our coding conventions for memory management. Utilizing dynamic memory management, ProteinDF calculates all matrices, with as many dimensions as basis sets, by using only a limited number of memory units. With this scheme, we can reversely calculate the maximum computation size from the memory size of the computing machine. When the number of basis sets is Norb, we can roughly estimate the memory size for each matrix to be 8 x Norb x Norb bytes. Note, however, that the OS or transferring processes also use some memory spaces. In particular, the OS may execute dynamic memory deletion on its own timing (even when the timing is expressly declared in the program). Note also that ProteinDF outputs intermediate files to disk storage in order to take full advantage of the memory. During computation, be sure to secure disk space sufficient for at least 200 matrices, each with the estimated computation size. ProteinDF saves all intermediate outputs for safe computation.

Others
======

ProteinDF is optimized for all-electron calculation of proteins made up of peptide chains. However, the automatic calculation methods and GUI are not well supported for more complicated calculations, such as those involving hetero molecules. To perform such calculations, thoroughly familiarize yourself with the program as well as applied methods, and modify the program manually. Refer to separate reference manuals for details. Finally, the following sections explain the points which, based on our experience of all-electron calculation of peptide chain proteins, seem most important:

Empirical rule on protein size
------------------------------

When performing DFT calculation of proteins, one of the most important first considerations is estimating the size of computation. In the Roothaan method, the eigenvalue problem of the Kohn-Sham equation is replaced with that of the matrix equation. This means that the size of the matrix can indicate the scale of computation. The following empirical rule traditionally defines the relationship between protein molecular mass and the total number of amino acid residues:

(Protein molecular mass) = 110 x (Number of amino acid residues)

This relational expression is well matched to the actual cases, due to the fact that proteins are made of peptide chains, and that the composition ratios of atoms which make up amino acids, such as hydrogen, carbon, nitrogen, oxygen, and sulfur, are averaged in large molecules. By expanding this rule, we can derive a proportional relationship on a protein, among the number of the amino acid residues (Nres), number of the atoms (Natom), number of the electrons (Nele), and number of the orbitals (Norb) (i.e. matrix dimensions). For example, in the case of the ProteinDF default set (VD equivalent), the following relationship is approximately established:

N\ :sub:`res` : N\ :sub:`atom` : N\ :sub:`ele` : N\ :sub:`orb` = 1 : 20 : 70 : 100

Here, it is empirically known that the half of Natom consists of hydrogen. This relational expression is useful to estimate computation size. Although Norb is dependent on the size of basis set, we can easily derive a similar proportional relation for different basis sets. Obtain a good estimation of computation size before actual computation, by deriving the relational expression according to the basis set.


Distortion in protein structure
-------------------------------

In general, we can obtain protein conformations from the Protein Data Bank (PDB). The PDB collects three-dimensional coordinate data of protein atoms experimentally determined with X-ray structural analysis, neutron scattering method, multi-dimensional NMR, and other methods. In the PDB, however, there are a number of data structurally distorted due to the characteristics of the experiment itself or by later data tuning. This distortion may cause adverse effects especially on DFT calculations. Meanwhile, it is not yet practical under the current computer resources to optimize the entire structure of proteins using only the DFT calculations based on all-electron CMO method. To optimize the structure, therefore, we recommend using the MM method, semi-empirical MO method such as MOZYME, QM/MM method, ONIOM method, and FMO method, etc. To check the distortion in protein structure, refer to the PROCHECK, etc.

Protein surface properties
--------------------------

Isolated systems have been used as the standard method for simulating chemical molecules. In our experience, however, the DFT cannot solve water-soluble proteins containing a number of dissociable groups on their surface in a vacuum environment. This fact may prove the accuracy of DFT calculations, but also proves troublesome for the person performing the calculation. Elucidating the surface property of proteins features among the forefront of current research subjects. In the surface property case, the most ideal means will involve appropriately arranging a number of solvent molecules (i.e. water molecules or buffer ions) around the protein. However, handling all molecules quantum-mechanically will significantly increase the computation size. Although there is an attraction in attempting that kind of calculation, an alternative means should be applied in practical computation, such as arranging classically-handled water molecules or counter ions in droplets around the proteins.
