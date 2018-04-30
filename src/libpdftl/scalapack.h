// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SCALAPACK_H
#define SCALAPACK_H

#ifdef F77_WITH_NO_UNDERSCORE
#define descinit_ descinit
#define infog2l_ infog2l
#define numroc_ numroc
#define pdelget_ pdelget
#define pdelset_ pdelset
#define pdgeadd_ pdgeadd
#define pdgemm_ pdgemm
#define pdsyev_ pdsyev
#define pdsyevd_ pdsyevd
#define pdpotrf_ pdpotrf
#define pdpotri_ pdpotri
#define pdsymm_ pdsymm
#define pdetrf_ pdetrf
#define pdetri_ pdetri
#define pdtran_ pdtran
#define pdgemv_ pdgemv
#define pdsymv_ pdsymv
#define pddot_ pddot
#endif  // F77_WITH_NO_UNDERSCORE

extern "C" {
// blacs_pinfo Reports a unique identifier (rank) for the calling process,
// and the total number of processes available in the BLACS environment.
// blacs_pinfo is available before other initialization or setup steps are
// taken.
//
// These routines take the following arguments:
//
// mypnum    Integer. (output)
//           An integer between 0 and (NPROCS-1) uniquely representing the
//           current process.
//
// nprocs    Integer. (output)
//           The total number of processes available to BLACS.
void Cblacs_pinfo(int* mypnum, int* nprocs);

// blacs_setup is relevent only for PVM based systems. In the SCSL
// implementation based on MPI, it is functionally equivilent to
// blacs_pinfo
void Cblacs_setup(int* mypnum, int* nprocs);

// blacs_get reports information about the BLACS state. The request may be
// tied to a particular context or handle, or may request system wide
// information, depending upon the nature of the request. blacs_get can be
// used to return the system default context needed to construct a first
// BLACS grid and its associated context. In this implementation the
// default context corresponds to the MPI default context, of
// MPI_COMM_WORLD.
//
// The routines take the following arguments:
//
// context   Integer. (input)
//           The context/handle which is the target of the enquiry. For
//           context-independent requests, this is unused.
//
// request   Integer. (input)
//           Which BLACS parameters should be returned.
//           request = 0 : Report the default system context.
//           request = 1 : Report the BLACS message ID range.
//           request = 2 : Report the debug level, BLACS was  compiled with.
//           request = 10 : Report the system context from which the
//                          specified context was constructed from.
//           request = 11 : Report the number of rings, multiring topology
//                          is presently using.
//           request = 12 : Report the number of branches, the general
void Cblacs_get(int context, int request, int* value);

// blacs_gridinit maps available processes onto a BLACS grid. The BLACS
// grid created is identified by a context handle that is subsequently
// used to identify that particular process grid among the many that may
// be generated. A BLACS grid must be created before using other BLACS
// calls, except where noted. blacs_gridinit is a collective or globally
// blocking operation amongst the participating processes.
//
// The routines take the following arguments:
//
// context   Integer. (input/output)
//           On input the handle of the system context to use (a default
//           system context can be obtained from blacs_get. On return, the
//           handle to the newly created grid.
//
// order     character*1. (input)
//           Specifies the nature of the mapping of processes onto the new
//           BLACS grid. A value that is not explicitly specified as
//           below, will default to row major ordering.
//           order = ’R’ : Use row major ordering
//           order = ’C’ : Use column major ordering
//
// np_row    Integer (input)
//           Specifies the number of rows in the process grid.
//
// np_col    Integer (input)
//           Specifies the number of columns in the process grid.
int Cblacs_gridinit(int* context, const char* order, int np_row, int np_col);

// blacs_gridinfo reports information about a BLACS context. If the
// specified context is invalid, all quantities are returned as -1.
//
// The routine takes the following arguments:
//
// context   Integer. (input)
//           The context/handle which is the target of the enquiry.
//
// np_row    Integer. (output)
//           The number of rows in the process grid.
//
// np_col    Integer. (output)
//           The number of columns in the process grid.
//
// my_row    Integer. (output)
//           The calling process row coordinate in the process grid.
//
// my_col    Integer. (output)
//           The calling process column coordinate in the process grid.
void Cblacs_gridinfo(const int context, int* np_row, int* np_col, int* my_row,
                     int* my_col);

// blacs_gridmap maps available processes onto a BLACS grid. The BLACS
// grid created is identified by a context handle that is subsequently
// used to identify that particular process grid among the many that may
// be generated. A BLACS grid must be created before using other BLACS
// calls, except where noted. blacs_gridmap is a collective or globally
// blocking operation amongst the participating processes. The 2D usermap
// array specifies which process corresponds to that location in the new
// 2D grid.
//
// The routines take the following arguments:
//
// context   Integer. (input/output)
//           On input the handle of the system context to use (a default
//           system context can be obtained from blacs_get. On return, the
//           handle to the newly created grid.
//
// usermap   Integer array of dimention (ld_usermap,np_col). (input)
//           This array contains process ranks that will be mapped onto
//           the corresponding location in the new process grid. The
//           process ranks can be obtained from e.g. blacs_pnum
//
// ld_usermap
//           Integer (input)
//           Specifies the leading dimention of the usermap array.
//
// np_row    Integer (input)
//           Specifies the number of rows in the process grid.
int Cblacs_gridmap(int* context, int* usermap, int ld_usermap, int np_row,
                   int np_col);

// blacs_gridexit frees resources associated with the specified context.
// After the call to blacs_gridexit the context is invalid. The numeric
// value of the context may be recycled if subsequent process grids are
// created with blacs_gridinit or blacs_gridmap.
//
// The routine takes the following arguments:
//
// context   Integer. (input)
//           The handle of the context to be freed.
void Cblacs_gridexit(int context);

// blacs_exit indicates that all work using the BLACS has been completed.
// The call will initiate destruction of the BLACS internal structures and
// free all memory possible.
//
// The routines take the following arguments:
//
// continue  Integer. (input)
//           If continue is set to a non zero value, the application will
//           be able to continue using the message passing subsystems
//           after the call. The application is then responsible for
//           shutting down MPI appropriately and in an orderly fashion.
//           Otherwise, BLACS will call MPI_FINALIZE() internally and
//           further communication using MPI will not be possible (If the
//           application has incomplete MPI operations [not instigated by
//           BLACS] at the time that this call is made then the program
//           state will be undefined).
void Cblacs_exit(int ncontinue);

// DESCINIT initializes the descriptor vector with the 8 input arguments
// M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
// *
// *  Arguments
// *  =========
// *
// *  DESC    (output) INTEGER array of dimension DLEN_.
// *          The array descriptor of a distributed matrix to be set.
// *
// *  M       (global input) INTEGER
// *          The number of rows in the distributed matrix. M >= 0.
// *
// *  N       (global input) INTEGER
// *          The number of columns in the distributed matrix. N >= 0.
// *
// *  MB      (global input) INTEGER
// *          The blocking factor used to distribute the rows of the
// *          matrix. MB >= 1.
// *
// *  NB      (global input) INTEGER
// *          The blocking factor used to distribute the columns of the
// *          matrix. NB >= 1.
// *
// *  IRSRC   (global input) INTEGER
// *          The process row over which the first row of the matrix is
// *          distributed. 0 <= IRSRC < NPROW.
// *
// *  ICSRC   (global input) INTEGER
// *          The process column over which the first column of the
// *          matrix is distributed. 0 <= ICSRC < NPCOL.
// *
// *  ICTXT   (global input) INTEGER
// *          The BLACS context handle, indicating the global context of
// *          the operation on the matrix. The context itself is global.
// *
// *  LLD     (local input)  INTEGER
// *          The leading dimension of the local array storing the local
// *          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
// *
// *  INFO    (output) INTEGER
// *          = 0: successful exit
// *          < 0: if INFO = -i, the i-th argument had an illegal value
// *
void descinit_(int* desc, const int* m, const int* n, const int* mb,
               const int* nb, const int* irsrc, const int* icsrc,
               const int* ictxt, const int* lld, int* info);

// *  INFOG2L computes the starting local indexes LRINDX, LCINDX corres-
// *  ponding to the distributed submatrix starting globally at the entry
// *  pointed by GRINDX, GCINDX. This routine returns the coordinates in
// *  the grid of the process owning the matrix entry of global indexes
// *  GRINDX, GCINDX, namely RSRC and CSRC.
// *
// *  Arguments
// *  =========
// *
// *  GRINDX    (global input) INTEGER
// *            The global row starting index of the submatrix.
// *
// *  GCINDX    (global input) INTEGER
// *            The global column starting index of the submatrix.
// *
// *  DESC      (input) INTEGER array of dimension DLEN_.
// *            The array descriptor for the underlying distributed matrix.
// *
// *  NPROW     (global input) INTEGER
// *            The total number of process rows over which the distributed
// *            matrix is distributed.
// *
// *  NPCOL     (global input) INTEGER
// *            The total number of process columns over which the
// *            distributed matrix is distributed.
// *
// *  MYROW     (local input) INTEGER
// *            The row coordinate of the process calling this routine.
// *
// *  MYCOL     (local input) INTEGER
// *            The column coordinate of the process calling this routine.
// *
// *  LRINDX    (local output) INTEGER
// *            The local rows starting index of the submatrix.
// *
// *  LCINDX    (local output) INTEGER
// *            The local columns starting index of the submatrix.
// *
// *  RSRC      (global output) INTEGER
// *            The row coordinate of the process that possesses the first
// *            row and column of the submatrix.
// *  CSRC      (global output) INTEGER
// *            The column coordinate of the process that possesses the
// *            first row and column of the submatrix.
void infog2l_(const int* grindx, const int* gcindx, const int* desc,
              const int* nprow, const int* npcol, const int* myrow,
              const int* mycol, int* lrindx, int* lcindx, int* rsrc, int* csrc);

// *  NUMROC computes the NUMber of Rows Or Columns of a distributed
// *  matrix owned by the process indicated by IPROC.
// *
// *  Arguments
// *  =========
// *
// *  N         (global input) INTEGER
// *            The number of rows/columns in distributed matrix.
// *
// *  NB        (global input) INTEGER
// *            Block size, size of the blocks the distributed matrix is
// *            split into.
// *
// *  IPROC     (local input) INTEGER
// *            The coordinate of the process whose local array row or
// *            column is to be determined.
// *
// *  ISRCPROC  (global input) INTEGER
// *            The coordinate of the process that possesses the first
// *            row or column of the distributed matrix.
// *
// *  NPROCS    (global input) INTEGER
// *            The total number processes over which the matrix is
// *            distributed.
int numroc_(const int* n, const int* nb, const int* iproc, const int* isrcproc,
            const int* nprocs);

// *  PDELGET sets alpha to the distributed matrix entry A( IA, JA ).
// *  The value of alpha is set according to the scope.
// *
// *  Arguments
// *  =========
// *
// *  SCOPE   (global input) CHARACTER*1
// *          The BLACS scope in which alpha is updated.
// *          If SCOPE = 'R', alpha is updated only in the process row
// *                          containing A( IA, JA ),
// *          If SCOPE = 'C', alpha is updated only in the process column
// *                          containing A( IA, JA ),
// *          If SCOPE = 'A', alpha is updated in all the processes of the
// *                          grid,
// *          otherwise alpha is updated only in the process containing
// *           A( IA, JA ).
// *
// *  TOP     (global input) CHARACTER*1
// *          The topology to be used if broadcast is needed.
// *
// *  ALPHA   (global output) DOUBLE PRECISION, the scalar alpha.
// *
// *  A       (local input) DOUBLE PRECISION pointer into the local memory
// *          to an array of dimension (LLD_A,*) containing the local
// *          pieces of the distributed matrix A.
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
void pdelget_(const char* SCOPE, const char* TOP, double* ALPHA,
              const double* A, const int* IA, const int* JA, const int* DESCA);

// *  PDELSET sets the distributed matrix entry A( IA, JA ) to ALPHA.
// *
// *  Arguments
// *  =========
// *
// *  A       (local output) DOUBLE PRECISION pointer into the local memory
// *          to an array of dimension (LLD_A,*) containing the local
// *          pieces of the distributed matrix A.
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  ALPHA   (local input) DOUBLE PRECISION
// *          The scalar alpha.
// *
void pdelset_(double* a, const int* ia, const int* ja, const int* desca,
              const double* alpha);

// *  PDGEADD  adds a matrix to another
// *
// *     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
// *
// *  where
// *
// *     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
// *
// *     op( X ) = X   or   op( X ) = X'.
// *
// *  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+N-1)   if TRANS = 'N',
// *                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'T',
// *                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'C'.
// *
// *  Alpha  and  beta  are scalars, sub( C ) and op( sub( A ) ) are m by n
// *  submatrices.
// *
// *  Notes
// *  =====
// *
// *  A description  vector  is associated with each 2D block-cyclicly dis-
// *  tributed matrix.  This  vector  stores  the  information  required to
// *  establish the  mapping  between a  matrix entry and its corresponding
// *  process and memory location.
// *
// *  In  the  following  comments,   the character _  should  be  read  as
// *  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
// *  block cyclicly distributed matrix.  Its description vector is DESC_A:
// *
// *  NOTATION         STORED IN       EXPLANATION
// *  ---------------- --------------- ------------------------------------
// *  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
// *  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
// *                                   the NPROW x NPCOL BLACS process grid
// *                                   A  is  distributed over. The context
// *                                   itself  is  global,  but  the handle
// *                                   (the integer value) may vary.
// *  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
// *                                   ted matrix A, M_A >= 0.
// *  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
// *                                   buted matrix A, N_A >= 0.
// *  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
// *                                   block of the matrix A, IMB_A > 0.
// *  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
// *                                   left   block   of   the  matrix   A,
// *                                   INB_A > 0.
// *  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
// *                                   bute the last  M_A-IMB_A  rows of A,
// *                                   MB_A > 0.
// *  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
// *                                   bute the last  N_A-INB_A  columns of
// *                                   A, NB_A > 0.
// *  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
// *                                   row of the matrix  A is distributed,
// *                                   NPROW > RSRC_A >= 0.
// *  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
// *                                   first column of  A  is  distributed.
// *                                   NPCOL > CSRC_A >= 0.
// *  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
// *                                   array  storing  the  local blocks of
// *                                   the distributed matrix A,
// *                                   IF( Lc( 1, N_A ) > 0 )
// *                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
// *                                   ELSE
// *                                      LLD_A >= 1.
// *
// *  Let K be the number of  rows of a matrix A starting at the global in-
// *  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
// *  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
// *  receive if these K rows were distributed over NPROW processes.  If  K
// *  is the number of columns of a matrix  A  starting at the global index
// *  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
// *  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
// *  these K columns were distributed over NPCOL processes.
// *
// *  The values of Lr() and Lc() may be determined via a call to the func-
// *  tion PB_Cnumroc:
// *  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A A )A spe   On   On(lfies C's
void pdgeadd_(const char* TRANS, const int* M, const int* N,
              const double* ALPHA, const double* A, const int* IA,
              const int* JA, const int* DESCA, const double* BETA, double* C,
              const int* IC, const int* JC, const int* DESCC);

// *  PDGEMM  performs one of the matrix-matrix operations
// *
// *     sub( C ) := alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C ),
// *
// *  where sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
// *
// *        op( X )  is one of
// *        op( X ) = X   or   op( X ) = X',
// *
// *  thus  op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+K-1)  if TRANSA = 'N',
// *                               A(IA:IA+K-1,JA:JA+M-1)' if TRANSA = 'T',
// *                               A(IA:IA+K-1,JA:JA+M-1)' if TRANSA = 'C',
// *
// *        op( sub( B ) ) denotes B(IB:IB+K-1,JB:JB+N-1)  if TRANSB = 'N',
// *                               B(IB:IB+N-1,JB:JB+K-1)' if TRANSB = 'T',
// *                               B(IB:IB+N-1,JB:JB+K-1)' if TRANSB = 'C',
// *
// *  alpha and beta are scalars, and sub( A ), sub( B ) and sub( C ) are
// *  distributed matrices, with op( sub( A ) ) an M-by-K distributed
// *  matrix, op( sub( B ) ) a K-by-N distributed matrix and  sub( C ) an
// *  M-by-N distributed matrix.
// *
// *  TRANSA  (global input) pointer to CHARACTER
// *          The form of op( A ) to be used in the matrix multiplication
// *          as follows:
// *
// *                TRANSA = 'N' or 'n',    op( A ) = A,
// *
// *                TRANSA = 'T' or 't',    op( A ) = A',
// *                TRANSA = 'C' or 'c',    op( A ) = A'.
// *
// *  TRANSB  (global input) pointer to CHARACTER
// *          The form of op( B ) to be used in the matrix multiplication
// *          as follows:
// *
// *                TRANSB = 'N' or 'n',    op( B ) = B,
// *                TRANSB = 'T' or 't',    op( B ) = B',
// *                TRANSB = 'C' or 'c',    op( B ) = B'.
// *
// *  M       (global input) pointer to INTEGER
// *          The number of rows of the distributed matrices op( sub( A ) )
// *          and sub( C ).  M >= 0.
// *
// *  N       (global input) pointer to INTEGER.
// *          The number of columns of the distributed matrices
// *           op( sub( B ) ) and sub( C ). N >= 0.
// *
// *  K       (global input) pointer to INTEGER.
// *          The number of columns of the distributed matrix
// *          op( sub( A ) ) and the number of rows of the distributed
// *          matrix op( B ). K >= 0.
// *
// *  ALPHA   (global input) pointer to DOUBLE PRECISION
// *          On entry, ALPHA specifies the scalar alpha.
// *
// *  A       (local input) DOUBLE PRECISION pointer into the local memory
// *          to an array of dimension (LLD_A, KLa), where KLa is
// *          LOCc(JA+K-1) when  TRANSA = 'N' or 'n',  and is LOCc(JA+M-1)
// *          otherwise.  Before entry, this array must contain the local
// *          pieces of the distributed matrix sub( A ).
// *
// *  IA      (global input) pointer to INTEGER
// *          The global row index of the submatrix of the distributed
// *          matrix A to operate on.
// *
// *  JA      (global input) pointer to INTEGER
// *          The global column index of the submatrix of the distributed
// *          matrix A to operate on.
// *
// *  DESCA   (global and local input) INTEGER array of dimension 8.
// *          The array descriptor of the distributed matrix A.
// *
// *  B       (local input) DOUBLE PRECISION pointer into the local memory
// *          to an array of dimension (LLD_B, KLb), where KLb is
// *          LOCc(JB+N-1) when  TRANSB = 'N' or 'n', and is LOCc(JB+K-1)
// *          otherwise. Before entry this array must contain the local
// *          pieces of the distributed matrix sub( B ).
// *
// *  IB      (global input) pointer to INTEGER
// *          The global row index of the submatrix of the distributed
// *          matrix B to operate on.
// *
// *  JB      (global input) pointer to INTEGER
// *          The global column index of the submatrix of the distributed
// *          matrix B to operate on.
// *
// *  DESCB   (global and local input) INTEGER array of dimension 8.
// *          The array descriptor of the distributed matrix B.
// *
// *  BETA    (global input) DOUBLE PRECISION
// *          On entry,  BETA  specifies the scalar  beta.  When  BETA  is
// *          supplied as zero then sub( C ) need not be set on input.
// *
// *  C       (local input/local output) DOUBLE PRECISION pointer into the
// *          local memory to an array of dimension (LLD_C, LOCc(JC+N-1)).
// *          Before entry, this array must contain the local pieces of the
// *          distributed matrix sub( C ). On exit, the distributed matrix
// *          sub( C ) is overwritten by the M-by-N distributed matrix
// *          alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C ).
// *
// *  IC      (global input) pointer to INTEGER
// *          The global row index of the submatrix of the distributed
// *          matrix C to operate on.
// *  JC      (global input) pointer to INTEGER
// *          The global column index of the submatrix of the distributed
// *          matrix C to operate on.
// *
// *  DESCC   (global and local input) INTEGER array of dimension 8.
// *          The array descriptor of the distributed matrix C.
void pdgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N,
             const int* K, const double* ALPHA, const double* A, const int* IA,
             const int* JA, const int* DESCA, const double* B, const int* IB,
             const int* JB, const int* DESCB, const double* BETA, double* C,
             const int* IC, const int* JC, const int* DESCC);

// *  PDSYEV computes all eigenvalues and, optionally, eigenvectors
// *  of a real symmetric matrix A by calling the recommended sequence
// *  of ScaLAPACK routines.
// *
// *  In its present form, PDSYEV assumes a homogeneous system and makes
// *  no checks for consistency of the eigenvalues or eigenvectors across
// *  the different processes.  Because of this, it is possible that a
// *  heterogeneous system may return incorrect results withiout any error
// *  messages.
// *
// *  Arguments
// *  =========
// *
// *     NP = the number of rows local to a given process.
// *     NQ = the number of columns local to a given process.
// *
// *  JOBZ    (global input) CHARACTER*1
// *          Specifies whether or not to compute the eigenvectors:
// *          = 'N':  Compute eigenvalues only.
// *          = 'V':  Compute eigenvalues and eigenvectors.
// *
// *  UPLO    (global input) CHARACTER*1
// *          Specifies whether the upper or lower triangular part of the
// *          symmetric matrix A is stored:
// *          = 'U':  Upper triangular
// *          = 'L':  Lower triangular
// *
// *  N       (global input) INTEGER
// *          The number of rows and columns of the matrix A.  N >= 0.
// *
// *  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
// *          global dimension (N, N), local dimension ( LLD_A,
// *          LOCc(JA+N-1) )
// *
// *          On entry, the symmetric matrix A.  If UPLO = 'U', only the
// *          upper triangular part of A is used to define the elements of
// *          the symmetric matrix.  If UPLO = 'L', only the lower
// *          triangular part of A is used to define the elements of the
// *          symmetric matrix.
// *
// *          On exit, the lower triangle (if UPLO='L') or the upper
// *          triangle (if UPLO='U') of A, including the diagonal, is
// *          destroyed.
// *
// *  IA      (global input) INTEGER
// *          A's global row index, which points to the beginning of the
// *          submatrix which is to be operated on.
// *
// *  JA      (global input) INTEGER
// *          A's global column index, which points to the beginning of
// *          the submatrix which is to be operated on.
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *          If DESCA( CTXT_ ) is incorrect, PDSYEV cannot guarantee
// *          correct error reporting.
// *
// *  W       (global output) DOUBLE PRECISION array, dimension (N)
// *          On normal exit, the first M entries contain the selected
// *          eigenvalues in ascending order.
// *
// *  Z       (local output) DOUBLE PRECISION array,
// *          global dimension (N, N), local dimension ( LLD_Z, LOCc(JZ+N-1) )
// *          If JOBZ = 'V', then on normal exit the first M columns of Z
// *          contain the orthonormal eigenvectors of the matrix
// *          corresponding to the selected eigenvalues.
// *          If JOBZ = 'N', then Z is not referenced.
// *
// *  IZ      (global input) INTEGER
// *          Z's global row index, which points to the beginning of the
// *          submatrix which is to be operated on.
// *
// *  JZ      (global input) INTEGER
// *          Z's global column index, which points to the beginning of
// *          the submatrix which is to be operated on.
// *
// *  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix Z.
// *          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
// *
// *  WORK    (local workspace/output) DOUBLE PRECISION array,
// *          dimension (LWORK)
// *          Version 1.0:  on output, WORK(1) returns the workspace
// *          needed to guarantee completion.
// *          If the input parameters are incorrect, WORK(1) may also be
// *          incorrect.
// *
// *          If JOBZ='N' WORK(1) = minimal=optimal amount of workspace
// *          If JOBZ='V' WORK(1) = minimal workspace required to
// *             generate all the eigenvectors.
// *
// *
// *  LWORK   (local input) INTEGER
// *          See below for definitions of variables used to define LWORK.
// *          If no eigenvectors are requested (JOBZ = 'N') then
// *             LWORK >= 5*N + SIZESYTRD + 1
// *          where
// *             SIZESYTRD = The workspace requirement for PDSYTRD
// *                         and is MAX( NB * ( NP +1 ), 3 * NB )
// *          If eigenvectors are requested (JOBZ = 'V' ) then
// *             the amount of workspace required to guarantee that all
// *             eigenvectors are computed is:
// *
// *             QRMEM = 2*N-2
// *             LWMIN = 5*N + N*LDC + MAX( SIZEMQRLEFT, QRMEM ) + 1
// *
// *          Variable definitions:
// *             NB = DESCA( MB_ ) = DESCA( NB_ ) =
// *                  DESCZ( MB_ ) = DESCZ( NB_ )
// *             NN = MAX( N, NB, 2 )
// *             DESCA( RSRC_ ) = DESCA( RSRC_ ) = DESCZ( RSRC_ ) =
// *                              DESCZ( CSRC_ ) = 0
// *             NP = NUMROC( NN, NB, 0, 0, NPROW )
// *             NQ = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
// *             NRC = NUMROC( N, NB, MYPROWC, 0, NPROCS)
// *             LDC = MAX( 1, NRC )
// *             SIZEMQRLEFT = The workspace requirement for PDORMTR
// *                           when it's SIDE argument is 'L'.
// *
// *          With MYPROWC defined when a new context is created as:
// *             CALL BLACS_GET( DESCA( CTXT_ ), 0, CONTEXTC )
// *             CALL BLACS_GRIDINIT( CONTEXTC, 'R', NPROCS, 1 )
// *             CALL BLACS_GRIDINFO( CONTEXTC, NPROWC, NPCOLC, MYPROWC,
// *                                  MYPCOLC )
// *
// *          If LWORK = -1, the LWORK is global input and a workspace
// *          query is assumed; the routine only calculates the minimum
// *          size for the WORK array.  The required workspace is returned
// *          as the first element of WORK and no error message is issued
// *          by PXERBLA.
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegal value, then
// *                INFO = -i.
// *          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
// *                converge in DSTEQR2 after a total of 30*N iterations.
// *                If INFO = N+1, then PDSYEV has detected heterogeneity
// *                by finding that eigenvalues were not identical across
// *                the process grid.  In this case, the accuracy of
// *                the results from PDSYEV cannot be guaranteed.
void pdsyev_(const char* jobz, const char* uplo, const int* n, double* a,
             const int* ia, const int* ja, const int* desca, double* w,
             double* z, const int* iz, const int* jz, const int* descz,
             double* work, const int* lwork, int* info);

// *  PDSYEVD computes  all the eigenvalues and eigenvectors
// *  of a real symmetric matrix A by calling the recommended sequence
// *  of ScaLAPACK routines.
// *
// *  In its present form, PDSYEVD assumes a homogeneous system and makes
// *  no checks for consistency of the eigenvalues or eigenvectors across
// *  the different processes.  Because of this, it is possible that a
// *  heterogeneous system may return incorrect results without any error
// *  messages.
// *
// *  Arguments
// *  =========
// *
// *     NP = the number of rows local to a given process.
// *     NQ = the number of columns local to a given process.
// *
// *  JOBZ    (input) CHARACTER*1
// *          = 'N':  Compute eigenvalues only;     (NOT IMPLEMENTED YET)
// *          = 'V':  Compute eigenvalues and eigenvectors.
// *
// *  UPLO    (global input) CHARACTER*1
// *          Specifies whether the upper or lower triangular part of the
// *          symmetric matrix A is stored:
// *          = 'U':  Upper triangular
// *          = 'L':  Lower triangular
// *
// *  N       (global input) INTEGER
// *          The number of rows and columns to be operated on, i.e. the
// *          order of the distributed submatrix sub( A ). N >= 0.
// *
// *  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
// *          global dimension (N, N), local dimension ( LLD_A,
// *          LOCc(JA+N-1) )
// *          On entry, the symmetric matrix A.  If UPLO = 'U', only the
// *          upper triangular part of A is used to define the elements of
// *          the symmetric matrix.  If UPLO = 'L', only the lower
// *          triangular part of A is used to define the elements of the
// *          symmetric matrix.
// *          On exit, the lower triangle (if UPLO='L') or the upper
// *          triangle (if UPLO='U') of A, including the diagonal, is
// *          destroyed.
// *
// *  IA      (global input) INTEGER
// *          A's global row index, which points to the beginning of the
// *          submatrix which is to be operated on.
// *
// *  JA      (global input) INTEGER
// *          A's global column index, which points to the beginning of
// *          the submatrix which is to be operated on.
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  W       (global output) DOUBLE PRECISION array, dimension (N)
// *          If INFO=0, the eigenvalues in ascending order.
// *
// *  Z       (local output) DOUBLE PRECISION array,
// *          global dimension (N, N),
// *          local dimension ( LLD_Z, LOCc(JZ+N-1) )
// *          Z contains the orthonormal eigenvectors
// *          of the symmetric matrix A.
// *
// *  IZ      (global input) INTEGER
// *          Z's global row index, which points to the beginning of the
// *          submatrix which is to be operated on.
// *
// *  JZ      (global input) INTEGER
// *          Z's global column index, which points to the beginning of
// *          the submatrix which is to be operated on.
// *
// *  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix Z.
// *          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
// *
// *  WORK    (local workspace/output) DOUBLE PRECISION array,
// *          dimension (LWORK)
// *          On output, WORK(1) returns the workspace required.
// *
// *  LWORK   (local input) INTEGER
// *          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
// *          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
// *          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
// *          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
// *
// *          If LWORK = -1, the LWORK is global input and a workspace
// *          query is assumed; the routine only calculates the minimum
// *          size for the WORK array.  The required workspace is returned
// *          as the first element of WORK and no error message is issued
// *          by PXERBLA.
// *
// *  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
//*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
// *
// *  LIWORK  (input) INTEGER
// *          The dimension of the array IWORK.
// *          LIWORK = 7*N + 8*NPCOL + 2
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegal value, then
// *                INFO = -i.
// *          > 0:  The algorithm failed to compute the INFO/(N+1) th
// *                eigenvalue while working on the submatrix lying in
// *                global rows and columns mod(INFO,N+1).
void pdsyevd_(const char* jobz, const char* uplo, const int* n, double* a,
              const int* ia, const int* ja, const int* desca, double* w,
              double* z, const int* iz, const int* jz, const int* descz,
              double* work, const int* lwork, int* iwork, const int* liwork,
              int* info);

// *  PDPOTRF computes the Cholesky factorization of an N-by-N real
// *  symmetric positive definite distributed matrix sub( A ) denoting
// *  A(IA:IA+N-1, JA:JA+N-1).
// *
// *  The factorization has the form
// *
// *            sub( A ) = U' * U ,  if UPLO = 'U', or
// *
// *            sub( A ) = L  * L',  if UPLO = 'L',
// *
// *  where U is an upper triangular matrix and L is lower triangular.
// *
// *  Arguments
// *  =========
// *
// *  UPLO    (global input) CHARACTER
// *          = 'U':  Upper triangle of sub( A ) is stored;
// *          = 'L':  Lower triangle of sub( A ) is stored.
// *
// *  N       (global input) INTEGER
// *          The number of rows and columns to be operated on, i.e. the
// *          order of the distributed submatrix sub( A ). N >= 0.
// *
// *  A       (local input/local output) DOUBLE PRECISION pointer into the
// *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
// *          On entry, this array contains the local pieces of the
// *          N-by-N symmetric distributed matrix sub( A ) to be factored.
// *          If UPLO = 'U', the leading N-by-N upper triangular part of
// *          sub( A ) contains the upper triangular part of the matrix,
// *          and its strictly lower triangular part is not referenced.
// *          If UPLO = 'L', the leading N-by-N lower triangular part of
// *          sub( A ) contains the lower triangular part of the distribu-
// *          ted matrix, and its strictly upper triangular part is not
// *          referenced. On exit, if UPLO = 'U', the upper triangular
// *          part of the distributed matrix contains the Cholesky factor
// *          U, if UPLO = 'L', the lower triangular part of the distribu-
// *          ted matrix contains the Cholesky factor L.
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegal value, then
// *                INFO = -i.
// *          > 0:  If INFO = K, the leading minor of order K,
// *                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
// *                the factorization could not be completed.
void pdpotrf_(const char* UPLO, const int* N, double* A, const int* IA,
              const int* JA, const int* DESCA, int* INFO);

// *  PDPOTRI computes the inverse of a real symmetric positive definite
// *  distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the
// *  Cholesky factorization sub( A ) = U**T*U or L*L**T computed by
// *  PDPOTRF.
// *  Arguments
// *  =========
// *
// *  UPLO    (global input) CHARACTER*1
// *          = 'U':  Upper triangle of sub( A ) is stored;
// *
// *          = 'L':  Lower triangle of sub( A ) is stored.
// *
// *  N       (global input) INTEGER
// *          The number of rows and columns to be operated on, i.e. the
// *          order of the distributed submatrix sub( A ). N >= 0.
// *
// *  A       (local input/local output) DOUBLE PRECISION pointer into the
// *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
// *          On entry, the local pieces of the triangular factor U or L
// *          from the Cholesky factorization of the distributed matrix
// *          sub( A ) = U**T*U or  L*L**T, as computed by PDPOTRF.
// *          On exit, the local pieces of the upper or lower triangle of
// *          the (symmetric) inverse of sub( A ), overwriting the input
// *          factor U or L.
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegal value, then
// *                INFO = -i.
// *          > 0:  If INFO = i, the (i,i) element of the factor U or L is
// *                zero, and the inverse could not be computed.
void pdpotri_(const char* UPLO, const int* N, double* A, const int* IA,
              const int* JA, const int* DESCA, int* INFO);

// PvSYMM performs one of the distributed matrix-matrix operations
//
//   * sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ), or
//   * sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ),
//
//   where sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
//
//   sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
//   A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R',
//   sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
//
//   Alpha and beta are scalars, sub( A ) is a symmetric distributed matrix and
//   sub( B ) and sub( C ) are M-by-N distributed matrices.
//
// Arguments
//
// SIDE
// (global input) CHARACTER
//   On entry, SIDE specifies whether the symmetric distributed matrix sub( A )
//   appears on the left or right in the operation as follows:
//
//   * SIDE = 'L' sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ),
//   * SIDE = 'R' sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ),
//
// UPLO
//   (global input) CHARACTER
//   On entry, UPLO specifies whether the upper or lower triangular part of the
//   symmetric distributed matrix sub( A )
// is to be referenced as follows:
//
//         * UPLO = 'U' Only the upper triangular part of the symmetric
//         distributed matrix is to be referenced.
//         * UPLO = 'L' Only the lower triangular part of the symmetric
//         distributed matrix is to be referenced.
//
// M
//   (global input) INTEGER
//   The number of rows to be operated on i.e., the number of rows of the
//   distributed submatrix sub( C ). M >= 0.
// N
//   (global input) INTEGER
//   The number of columns to be operated on i.e the number of columns of the
//   distributed submatrix sub( C ). N >= 0.
// ALPHA
//   (global input) REAL/COMPLEX
//   On entry, ALPHA specifies the scalar alpha.
// A
//   (local input) array of dimension (LLD_A, LOCq(JA+NA-1))
//   Before entry this array contains the local pieces of the symmetric
//   distributed matrix sub( A ),
// such that when UPLO = 'U', the NA-by-NA upper triangular part of the
// distributed matrix sub( A ) must contain the upper triangular part of the
// symmetric distributed matrix and the strictly lower triangular part of sub( A
// ) is not referenced, and when UPLO = 'L', the NA-by-NA lower triangular part
// of the distributed matrix sub( A ) must contain the lower triangular part of
// the symmetric distributed matrix and the strictly lower triangular part of
// sub( A ) is not referenced.
// IA
//   (global input) INTEGER
//     The global row index of the submatrix of the distributed matrix A to
//     operate on.
// JA
//   (global input) INTEGER
//     The global column index of the submatrix of the distributed matrix A to
//     operate on.
// DESCA
//   (global and local input) INTEGER array of dimension 8
//     The array descriptor of the distributed matrix A.
// B
//   (local input) array of dimension (LLD_B, LOCq(JB+N-1))
//   Before entry, this array contains the local pieces of the distributed
//   matrix sub( B ).
// IB
//   (global input) INTEGER
//     The global row index of the submatrix of the distributed matrix B to
//     operate on.
// JB
//   (global input) INTEGER
//     The global column index of the submatrix of the distributed matrix B to
//     operate on.
// DESCB
//   (global and local input) INTEGER array of dimension 8
//     The array descriptor of the distributed matrix B.
// BETA
//   (global input) REAL/COMPLEX
//   On entry, BETA specifies the scalar beta. When BETA is supplied as zero
//   then sub( C ) need not be set on input.
// C
//   (local input/local output) array of dimension (LLD_C, LOCq(JC+N-1))
//   Before entry, this array must contain the local pieces of the distributed
//   matrix sub( C ).
// On exit, the distributed matrix sub( C ) is overwritten by the M-by-N updated
// distributed matrix.
//  IC
//  (global input) INTEGER
//     The global row index of the submatrix of the distributed matrix C to
//     operate on.
// JC
//  (global input) INTEGER
//     The global column index of the submatrix of the distributed matrix C to
//     operate on.
// DESCC
//  (global and local input) INTEGER array of dimension 8
//     The array descriptor of the distributed matrix C.
void pdsymm_(const char* SIDE, const char* UPLO, const int* M, const int* N,
             const double* ALPHA, const double* A, const int* IA, const int* JA,
             const int* DESCA, const double* B, const int* IB, const int* JB,
             const int* DESCB, const double* BETA, double* C, const int* IC,
             const int* JC, const int* DESCC);

// *  PDGETRF computes an l value, then
// *                INFO = -i.
// *          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
// *                The factorization has been completed, but the factor U
// *                is exactly singular, and division by zero will occur if
// *                it is used to solve a system of equations.
// *
// *  M       (global input) INTEGER
// *          The number of rows to be operated on, i.e. the number of rows
// *          of the distributed submatrix sub( A ). M >= 0.
// *
// *  N       (global input) INTEGER
// *          The number of columns to be operated on, i.e. the number of
// *          columns of the distributed submatrix sub( A ). N >= 0.
// *
// *  A       (local input/local output) DOUBLE PRECISION pointer into the
// *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
// *          On entry, this array contains the local pieces of the M-by-N
// *          distributed matrix sub( A ) to be factored. On exit, this
// *          array contains the local pieces of the factors L and U from
// *          the factorization sub( A ) = P*L*U; the unit diagonal ele-
// *          ments of L are not stored.
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
// *          This array contains the pivoting information.
// *          IPIV(i) -> The global row local row i was swapped with.
// *          This array is tied to the distributed matrix A.
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegaAL           ICEIL
void pdgetrf_(const int* M, const int* N, double* A, const int* IA,
              const int* JA, const int* DESCA, int* IPIV, int* INFO);

// *  N       (global input) INTEGER
// *          The number of rows and columns to be operated on, i.e. the
// *          order of the distributed submatrix sub( A ). N >= 0.
// *
// *  A       (local input/local output) DOUBLE PRECISION pointer into the
// *          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
// *          On entry, the local pieces of the L and U obtained by the
// *          factorization sub( A ) = P*L*U computed by PDGETRF. On
// *          exit, if INFO = 0, sub( A ) contains the inverse of the
// *          original distributed matrix sub( A ).
// *
// *  IA      (global input) INTEGER
// *          The row index in the global array A indicating the first
// *          row of sub( A ).
// *
// *  JA      (global input) INTEGER
// *          The column index in the global array A indicating the
// *          first column of sub( A ).
// *
// *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
// *          The array descriptor for the distributed matrix A.
// *
// *  IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A
// *          keeps track of the pivoting information. IPIV(i) is the
// *          global row index the local row i was swapped with.  This
// *          array is tied to the distributed matrix A.
// *
// *  WORK    (local workspace/local output) DOUBLE PRECISION array,
// *                                                     dimension (LWORK)
// *          On exit, WORK(1) returns the minimal and optimal LWORK.
// *
// *  LWORK   (local or global input) INTEGER
// *          The dimension of the array WORK.
// *          LWORK is local input and must be at least
// *          LWORK = LOCr(N+MOD(IA-1,MB_A))*NB_A. WORK is used to keep a
// *          copy of at most an entire column block of sub( A ).
// *
// *          If LWORK = -1, then LWORK is global input and a workspace
// *          query is assumed; the routine only calculates the minimum
// *          and optimal size for all work arrays. Each of these
// *          values is returned in the first entry of the corresponding
// *          work array, and no error message is issued by PXERBLA.
// *
// *  IWORK   (local workspace/local output) INTEGER array,
// *                                                    dimension (LIWORK)
// *          On exit, IWORK(1) returns the minimal and optimal LIWORK.
// *
// *  LIWORK  (local or global input) INTEGER
// *          The dimension of the array IWORK used as workspace for
// *          physically transposing the pivots.
// *          LIWORK is local input and must be at least
// *          if NPROW == NPCOL then
// *            LIWORK = LOCc( N_A + MOD(JA-1, NB_A) ) + NB_A,
// *          else
// *            LIWORK =  LOCc( N_A + MOD(JA-1, NB_A) ) +
// *                      MAX( CEIL(CEIL(LOCr(M_A)/MB_A)/(LCM/NPROW)),
// *                           NB_A )
// *              where LCM is the least common multiple of process
// *              rows and columns (NPROW and NPCOL).
// *          end if
// *
// *          If LIWORK = -1, then LIWORK is global input and a workspace
// *          query is assumed; the routine only calculates the minimum
// *          and optimal size for all work arrays. Each of these
// *          values is returned in the first entry of the corresponding
// *          work array, and no error message is issued by PXERBLA.
// *
// *  INFO    (global output) INTEGER
// *          = 0:  successful exit
// *          < 0:  If the i-th argument is an array and the j-entry had
// *                an illegal value, then INFO = -(i*100+j), if the i-th
// *                argument is a scalar and had an illegal value, then
// *                INFO = -i.
// *          > 0:  If INFO = K, U(IA+K-1,IA+K-1) is exactly zero; the
// *                matrix is singular and its inverse could not be
// *                computed.
void pdgetri_(const int* N, double* A, const int* IA, const int* JA,
              const int* DESCA, const int* IPIV, double* WORK, const int* LWORK,
              int* IWORK, const int* LIWORK, int* INFO);

// NAME
//             PDTRAN  -  Transpose  a matrix
//
// PDTRAN performs the following matrix computation:
// C<--beta*C+alpha*op(A)
//
// SYNOPSIS
// SUBROUTINE  PDTRAN(  M,  N,  ALPHA,  A, IA, JA, DESCA, BETA, C, IC, JC,
//             DESCC )
//
//   INTEGER         IA,IC,JA,JC,M,N
//
//   DOUBLE          ALPHA,BETA
//
//   INTEGER         DESCA,DESCC
//
//   DOUBLE          A ( * ), C ( * )
//
// PURPOSE
//        PDTRAN  transposes a matrix
//   sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
//   where:
// sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1)
//   sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1)
//   Thus,  op(  sub(  A  )  )  denotes  A(IA:IA+N-1,JA:JA+M-1)'.  Beta is a
//        scalar, sub( C ) is an m by n submatrix, and sub( A ) is an n by m
//        sub- matrix.
void pdtran_(const int* M, const int* N, const double* ALPHA, const double* A,
             const int* IA, const int* JA, const int* DESCA, const double* BETA,
             double* C, const int* IC, const int* JC, const int* DESCC);

// PvGEMV
// SUBROUTINE PvGEMV( TRANS, M, N, ALPHA, A, IA, JA, DESCA,
// X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY )
//
// Purpose
//
// PvGEMV performs one of the distributed matrix-vector operations
//
// * sub( Y ) := alpha*sub( A ) * sub( X ) + beta*sub( Y ), or
//   * sub( Y ) := alpha*sub( A )' * sub( X ) + beta*sub( Y ),
//
// where sub( A ) denotes A(IA:IA+M-1,JA:JA+N-1),
//
//       sub( X ) denotes if TRANS = 'N',
//                      X(IX:IX,JX:JX+N-1), if INCX = M_X,
//                      X(IX:IX+N-1,JX:JX), if INCX = 1 and INCX <> M_X,
//                    else
//                      X(IX:IX,JX:JX+M-1), if INCX = M_X,
//                      X(IX:IX+M-1,JX:JX), if INCX = 1 and INCX <> M_X,
//                    end if
//       sub( Y ) denotes if trans = 'N',
//                      Y(IY:IY,JY:JY+M-1), if INCY = M_Y,
//                      Y(IY:IY+M-1,JY:JY), if INCY = 1 and INCY <> M_Y,
//                    else
//                      Y(IY:IY,JY:JY+N-1), if INCY = M_Y,
//                      Y(IY:IY+N-1,JY:JY), if INCY = 1 and INCY <> M_Y,
//                    end if
//
// alpha and beta are scalars, and sub( X ) and sub( Y ) are distributed vectors
// and sub( A ) is a M-by-N distributed submatrix.
void pdgemv_(const char* trans, const int* M, const int* N, const double* alpha,
             const double* A, const int* IA, const int* JA, const int* DESCA,
             const double* X, const int* IX, const int* JX, const int* DESCX,
             const int* INCX, const double* beta, double* Y, const int* IY,
             const int* JY, const int* DESCY, const int* INCY);

void pdsymv_(const char* UPLO, const int* N, const double* ALPHA,
             const double* A, const int* IA, const int* JA, const int* DESCA,
             const double* X, const int* IX, const int* JX, const int* DESCX,
             const int* INCX, const double* BETA, double* Y, const int* IY,
             const int* JY, const int* DESCY, const int* INCY);

void pddot_(const int* N, double* DOT, const double* X, const int* IX,
            const int* JX, const int* DESCX, const int* INCX, const double* Y,
            const int* IY, const int* JY, const int* DESCY, const int* INCY);
}

#endif  // SCALAPACH_H
