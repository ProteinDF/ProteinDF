TESTS = PdfTests
INCLUDES = -I$(top_srcdir)/src/pdf -I$(includedir)

if CPPUNIT
PDF_TESTS = PdfTests
else
PDF_TESTS =
endif
noinst_PROGRAMS = $(PDF_TESTS)

PDF_OBJS = \
	Fl_UserinputX.o \
	DfFunctional.o	\
	DfFunctional_Slater.o	\
	DfFunctional_Becke88.o	\
	DfFunctional_VWN.o	\
	DfFunctional_LYP.o	\
	DfFunctional_VWN3.o	\
	DfFunctional_Becke88_ExceptLDA.o \
	DfFunctional_B3LYP.o	\
	DfFunctional_PW91X.o

PdfTests_SOURCES = \
	digit.cpp							\
	PdfTests.cpp \
	TlUtilsTest.h			TlUtilsTest.cpp			\
	TlPositionTest.h		TlPositionTest.cpp		\
	TlParameterTest.h		TlParameterTest.cpp		\
	TlStringTokenizerTest.h		TlStringTokenizerTest.cpp	\
	TlMatrixTest.h			TlMatrixTest.cpp		\
	TlSymmetricMatrixTest.h		TlSymmetricMatrixTest.cpp	\
	TlSparseMatrixTest.h		TlSparseMatrixTest.cpp		\
	TlSparseSymmetricMatrixTest.h	TlSparseSymmetricMatrixTest.cpp	\
	TlVectorTest.h			TlVectorTest.cpp		\
	TlFileMatrixTest.h		TlFileMatrixTest.cpp		\
	TlFileSymmetricMatrixTest.h	TlFileSymmetricMatrixTest.cpp	\
	TlPartialSymmetricMatrixTest.h	TlPartialSymmetricMatrixTest.cpp\
	TlMemManagerTest.h		TlMemManagerTest.cpp		\
	TlMsgPackTest.h			TlMsgPackTest.cpp		\
	TlSerializeDataTest.h		TlSerializeDataTest.cpp		\
	Fl_UserinputXTest.h		Fl_UserinputXTest.cpp		\
	DfFunctional_SlaterTest.h	DfFunctional_SlaterTest.cpp	\
	DfFunctional_Becke88Test.h	DfFunctional_Becke88Test.cpp	\
	DfFunctional_VWNTest.h		DfFunctional_VWNTest.cpp	\
	DfFunctional_LYPTest.h		DfFunctional_LYPTest.cpp	\
	DfFunctional_B3LYPTest.h	DfFunctional_B3LYPTest.cpp	\
	DfFunctional_PW91XTest.h	DfFunctional_PW91XTest.cpp

PdfTests_LDFLAGS = -L$(libdir)
PdfTests_LDADD   = ../pdf/libpdf.la ../pdflib/libpdftl.la -lcppunit -ldl @LAPACK_LIBS@ @BLAS_LIBS@ 



