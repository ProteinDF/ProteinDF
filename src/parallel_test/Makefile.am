TESTS = parallelTest \
	TlDistributeMatrixTest \
	TlDistributeSymmetricMatrixTest \
	TlDistributeVectorTest \
	commtest \
	comm_mpi_test \
	comm_msgpack_test \
	scabench
INCLUDES = -I$(top_srcdir)/src/pdf -I$(includedir)

if PARALLEL
PPDF_TESTS = $(TESTS)
else
PPDF_TESTS = 
endif
noinst_PROGRAMS = $(PPDF_TESTS)

parallelTest_SOURCES = main_ParallelTest.cpp
parallelTest_LDADD   = ../pdflib/libppdftl.la ../pdflib/libpdftl.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

TlDistributeMatrixTest_SOURCES = main_TlDistributeMatrixTest.cpp
TlDistributeMatrixTest_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

TlDistributeSymmetricMatrixTest_SOURCES = main_TlDistributeSymmetricMatrixTest.cpp
TlDistributeSymmetricMatrixTest_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

TlDistributeVectorTest_SOURCES = main_TlDistributeVectorTest.cpp
TlDistributeVectorTest_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

commtest_SOURCES = main_commtest.cpp
commtest_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

comm_mpi_test_SOURCES = main_mpi_test.cpp
comm_mpi_test_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

comm_msgpack_test_SOURCES = main_comm_msgpack_test.cpp
comm_msgpack_test_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

scabench_SOURCES = main_scalapack_benchmark.cpp
scabench_LDADD = ../pdflib/libppdftl.la ../pdflib/libpdftl.la ../pdf/libpdf.la ../pdf/libppdf.la \
 @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@  

