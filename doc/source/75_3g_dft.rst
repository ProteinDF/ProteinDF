***********************************************************************
The third-generation density functional calculation method
***********************************************************************


Outline
***********************************************************************

The third-generation density functional theory is an indispensable method
for achieving canonical molecular orbitals of large molecules in ProteinDF.
Please refer to the paper [ref:Hirano]_ for details. Here is an overview.

In general, the computational complexity of a canonical molecular orbital calculation is :math:`N`, where N is the total number of basis functions.
Formally, it is expressed as :math:`O(N^{4})`.
This calculation size dependence is due to the four-center two-electron repulsive integral required to calculate the Coulomb and exchange terms.


By using the cutoff method to ignore small integrals in advance,
the computational complexity can be reduced to :math:`O(N^{2})` or even :math:`O(N)` by using empirical parameters.
However, no :math:`O(N)` method has yet been found that can stably achieve canonical molecular orbital calculations for any large molecule,
including proteins with a variety of properties.


In the third-generation density functional theory, the computational complexity is reduced by mathematically exact cutoff using the Cholesky decomposition.
The third-generation density functional method uses the Cholesky decomposition and mathematically exact cutoff to reduce the computational complexity and to achieve efficient parallel computation even on the current mainstream distributed parallel computers.


The third-generation density functional theory can be said to be an efficient storage of the molecular integrals
in the File method in the form of Cholesky vectors.
Nevertheless, the data size of Cholesky vectors for large molecules is still large,
which is a drawback of the third-generation density functional theory.



Calculation of Coulomb and exchange terms in the third-generation density functional calculation method
*******************************************************************************************************


Consider the supermatrix :math:`\pmb{V}` with the 4-center integral :math:`\langle pq | rs \rangle` as its matrix element.
Since this is a positive definite symmetric matrix, it can be Cholesky decomposed.


.. math::

    V^{ERI}_{pq,rs} = \langle pq | rs \rangle \approx \sum_{K=1}^{M}{L^{ERI}_{K,pq} L^{ERI}_{K,rs}}


Using the obtained Cholesky vector :math:`\pmb{L}`, the Coulomb and Fock exchange terms can be obtained as follows


.. math::

    J_{pq} = \sum_{rs}{P_{rs} \langle pq | rs \rangle} \approx \sum_{rs}\sum_{I}{L_{I,pq} L_{I,rs} P_{rs}}


.. math::

    K_{pq} &= - \frac{1}{2} \sum_{rs}{P_{rs} \langle pr | qs \rangle} \approx - \frac{1}{2} \sum_{I}\sum_{i}{X_{I,pi} X_{I,qi}}

    X_{I,pi} &= \sum_{r} {L^{ERI}_{I,pr} Q_{ri}}

    P_{rs} &= \sum_{i}{Q_{ri} Q_{si}}


The computational precision of the Cholesky decomposition, :math:`\delta`, can then be obtained as follows
The decomposition can be mathematically exact with any computational precision determined by the user, and the rank of the matrix can be reduced accordingly.
In other words, the size of the Cholesky vector :math:`\pmb{L}` can be reduced according to the desired computational precision.


.. math::

    \left| \langle pq | rs \rangle - \sum_{K=1}^{M}{L_{K,pq} L_{K,rs}} \right| \leq \delta


You can also reduce the size of the Cholesky vector :math:`\pmb{L}` by sorting out the diagonal components of the supermatrix :math:`\pmb{V}` that have values larger than the threshold :math:`\tau`.
We can also reduce the size of the Cholesky vector :math:`\pmb{L}` by selecting those whose diagonal components are larger than the threshold :math:`\pmb{V}`.
This is equivalent to the commonly used cutoff method with Schwarz's inequality.


.. math::

    \langle pq | pq \rangle \ge \tau


Implementation of the third-generation density functional calculation method in the ProteinDF
*********************************************************************************************


In recent years, computers have greatly improved their performance through the use of caches.
The direction of data access for Cholesky vectors differs between the generation of Cholesky vectors
and the calculation of Coulomb and exchange terms using Cholesky vectors.
When generating Cholesky vectors, the data is accessed in the line direction,
while when calculating Coulomb and exchange terms, the data is accessed in the column direction.
To make the best use of the computer cache, ProteinDF switches the data access direction before and during SCF iterations.


When switching the data access direction for Cholesky vectors and when calculating Coulomb and exchange terms, file I/O is a major performance issue.
ProteinDF allows you to specify the storage area of Cholesky vectors with the `local_temp_path` keyword.
By specifying a device that provides fast file I/O, such as SSD, you can achieve a high-performance calculation.


Reference
***********************************************************************

.. [ref:Hirano] T. Hirano, F. Sato, "A third-generation density-functional-theory-based method for calculating canonical molecular orbitals of large molecules", Phys. Chem. Chem. Phys., 16(28), 14496-14503 (2014); https://doi.org/10.1039/C3CP55514C

