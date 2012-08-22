SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@
PARALLEL_LIBS = ../pdf/libppdf.la ../pdf/libpdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 		\
	testERI		\
	testDensFit	\
	testXC

PARALLEL_PROGS = 	\
	testERIp	\
	testDensFitp	\
	testXCp


testDensFit_SOURCES = main_testDensFit.cpp
testDensFit_LDADD = $(SERIAL_LIBS)

testERI_SOURCES = main_testERI.cpp
testERI_LDADD = $(SERIAL_LIBS)

testXC_SOURCES = testXC_template.h main_testXC.cpp
testXC_LDADD = $(SERIAL_LIBS)

testDensFitp_SOURCES = main_testDensFit_parallel.cpp
testDensFitp_LDADD = $(PARALLEL_LIBS)

testERIp_SOURCES = main_testERI_parallel.cpp
testERIp_LDADD = $(PARALLEL_LIBS)

testXCp_SOURCES = testXC_template.h main_testXC_parallel.cpp
testXCp_LDADD = $(PARALLEL_LIBS)


if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
