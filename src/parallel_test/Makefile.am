SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@
PARALLEL_LIBS = ../pdf/libppdf.la ../pdf/libpdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@ @FLIBS@

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
parallelTest_LDADD   = $(PARALLEL_LIBS)

TlDistributeMatrixTest_SOURCES = main_TlDistributeMatrixTest.cpp
TlDistributeMatrixTest_LDADD = $(PARALLEL_LIBS)

TlDistributeSymmetricMatrixTest_SOURCES = main_TlDistributeSymmetricMatrixTest.cpp
TlDistributeSymmetricMatrixTest_LDADD = $(PARALLEL_LIBS)

TlDistributeVectorTest_SOURCES = main_TlDistributeVectorTest.cpp
TlDistributeVectorTest_LDADD = $(PARALLEL_LIBS)

commtest_SOURCES = main_commtest.cpp
commtest_LDADD = $(PARALLEL_LIBS)

comm_mpi_test_SOURCES = main_mpi_test.cpp
comm_mpi_test_LDADD = $(PARALLEL_LIBS)

comm_msgpack_test_SOURCES = main_comm_msgpack_test.cpp
comm_msgpack_test_LDADD = $(PARALLEL_LIBS)

scabench_SOURCES = main_scalapack_benchmark.cpp
scabench_LDADD = $(PARALLEL_LIBS)
