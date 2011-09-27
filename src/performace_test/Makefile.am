SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@  
PARALLEL_LIBS = ../pdf/libppdf.la ../pdf/libpdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 		\
	ptERI		\
	ptDensFit	\
	ptXC

PARALLEL_PROGS = 	\
	pptERI		\
	pptDensFit


ptERI_SOURCES = main_ptERI.cpp
ptERI_LDADD = $(SERIAL_LIBS)

pptERI_SOURCES = main_ptERI_parallel.cpp
pptERI_LDADD = $(PARALLEL_LIBS)


ptDensFit_SOURCES = main_ptDensFit.cpp
ptDensFit_LDADD = $(SERIAL_LIBS)


pptDensFit_SOURCES = main_ptDensFit_parallel.cpp
pptDensFit_LDADD = $(PARALLEL_LIBS)

ptXC_SOURCES = main_ptXC.cpp
ptXC_LDADD = $(SERIAL_LIBS)


if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
