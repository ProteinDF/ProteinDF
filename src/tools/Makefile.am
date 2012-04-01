#	RLMO					
#RLMO_SOURCES = \
#	main_RLMO.cpp
#RLMO_LDADD = $(SERIAL_LIBS)

SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@ @FCLIBS@
PARALLEL_LIBS = ../pdf/libpdf.la ../pdf/libppdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FCLIBS@ @MPILIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 				\
	pdfkwd 					\
	pdfvtr2txt 				\
	pdfmat2txt 				\
	pdfdiffmat 				\
	pdfdiffvtr 				\
	topPdfMatrix				\
	catPdfMatrix addPdfMatrix	\
	combineDensityMatrix 		\
	pdfTop					\
	pdfOrbInfo 				\
	pdfMatrixCheck 			\
	mullikenPop 			\
	charge 					\
	mmap_bench 				\
	mmap_bench_matrix_calc 	\
	localize 				\
	merge 					\
	makeCrd					\
	makeMoFld				\
	makeEspFld makeDensFld 	\
	makeDensEspFld 			\
	transdim				\
	diagonal				\
	cholesky				\
	inverse					\
	multiple				\
	transpose				\
	getMaxValue				\
	rotateLCAO				\
	checkERI


PARALLEL_PROGS = \
	pdiagonal \
	plocalize \
	mergeLocalMatrix \
	divideGlobalMatrix 

pdfkwd_SOURCES = pdfkwd.cpp
pdfkwd_LDADD = $(SERIAL_LIBS)

pdfvtr2txt_SOURCES = pdfvtr2txt.cpp
pdfvtr2txt_LDADD = $(SERIAL_LIBS)

pdfmat2txt_SOURCES = pdfmat2txt.cpp
pdfmat2txt_LDADD = $(SERIAL_LIBS)

pdfdiffmat_SOURCES = pdfdiffmat.cpp
pdfdiffmat_LDADD = $(SERIAL_LIBS)

pdfdiffvtr_SOURCES = pdfdiffvtr.cpp
pdfdiffvtr_LDADD = $(SERIAL_LIBS)


topPdfMatrix_SOURCES = main_topPdfMatrix.cpp
topPdfMatrix_LDADD = $(SERIAL_LIBS)


catPdfMatrix_SOURCES = catPdfMatrix.cpp
catPdfMatrix_LDADD = $(SERIAL_LIBS)


addPdfMatrix_SOURCES = addPdfMatrix.cpp
addPdfMatrix_LDADD = $(SERIAL_LIBS)


combineDensityMatrix_SOURCES = main_combineDensityMatrix.cpp
combineDensityMatrix_LDADD = $(SERIAL_LIBS)

pdfTop_SOURCES = main_top.cpp
pdfTop_LDADD = $(SERIAL_LIBS)

pdfOrbInfo_SOURCES = pdfOrbInfo.cpp
pdfOrbInfo_LDADD = $(SERIAL_LIBS)


pdfMatrixCheck_SOURCES = pdfMatrixCheck.cpp
pdfMatrixCheck_LDADD = $(SERIAL_LIBS)


mullikenPop_SOURCES = main_mullikenPop.cpp
mullikenPop_LDADD = $(SERIAL_LIBS)


charge_SOURCES = charge.cpp
charge_LDADD = $(SERIAL_LIBS)


mmap_bench_SOURCES = mmap_bench.cpp
mmap_bench_LDADD = $(SERIAL_LIBS)


mmap_bench_matrix_calc_SOURCES = mmap_bench_matrix_calc.cpp
mmap_bench_matrix_calc_LDADD = $(SERIAL_LIBS)


localize_SOURCES = \
	DfLocalize.h	DfLocalize.cpp 	\
	main_localize.cpp
localize_LDADD = $(SERIAL_LIBS)


plocalize_SOURCES = \
	DfLocalize.h		DfLocalize.cpp	\
	DfLocalize_Parallel.h	DfLocalize_Parallel.cpp \
	main_plocalize.cpp
plocalize_LDADD = $(PARALLEL_LIBS)


merge_SOURCES = merge.cpp
merge_LDADD = $(SERIAL_LIBS)


makeCrd_SOURCES = \
	makeField_common.h	makeField_common.cpp	\
	main_makeCrd.cpp
makeCrd_LDADD = $(SERIAL_LIBS)


makeMoFld_SOURCES = \
	makeField_common.h	makeField_common.cpp	\
	TlMoField.h		TlMoField.cpp		\
	makeMoFld.cpp
makeMoFld_LDADD = $(SERIAL_LIBS)


makeEspFld_SOURCES = \
	makeField_common.h	makeField_common.cpp	\
	TlEspField.h		TlEspField.cpp		\
	makeEspFld.cpp
makeEspFld_LDADD = $(SERIAL_LIBS)


makeDensFld_SOURCES = \
	makeField_common.h	makeField_common.cpp	\
	TlDensField.h		TlDensField.cpp		\
	makeDensFld.cpp
makeDensFld_LDADD = $(SERIAL_LIBS)


makeDensEspFld_SOURCES = \
	makeField_common.h	makeField_common.cpp	\
	TlDensField.h		TlDensField.cpp		\
	TlEspField.h		TlEspField.cpp		\
	makeDensEspFld.cpp
makeDensEspFld_LDADD = $(SERIAL_LIBS)


transdim_SOURCES = \
	main_transdim.cpp
transdim_LDADD = $(SERIAL_LIBS)


inverse_SOURCES = \
	main_inverse.cpp
inverse_LDADD = $(SERIAL_LIBS)


diagonal_SOURCES = \
	main_diagonal.cpp
diagonal_LDADD = $(SERIAL_LIBS)


cholesky_SOURCES = \
	main_choleskydecompose.cpp
cholesky_LDADD = $(SERIAL_LIBS)


multiple_SOURCES = \
	main_multiple.cpp
multiple_LDADD = $(SERIAL_LIBS)


transpose_SOURCES = \
	main_transpose.cpp
transpose_LDADD = $(SERIAL_LIBS)


pdiagonal_SOURCES = main_pdiagonal.cpp
pdiagonal_LDADD = $(PARALLEL_LIBS)


getMaxValue_SOURCES = main_getMaxValue.cpp
getMaxValue_LDADD = $(SERIAL_LIBS)

rotateLCAO_SOURCES = \
	TlRotateLCAO.h		TlRotateLCAO.cpp	\
	main_rotateLCAO.cpp
rotateLCAO_LDADD = $(SERIAL_LIBS)

mergeLocalMatrix_SOURCES = main_mergeLocalMatrix.cpp
mergeLocalMatrix_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)

divideGlobalMatrix_SOURCES = main_divideGlobalMatrix.cpp
divideGlobalMatrix_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)

checkERI_SOURCES = main_checkERI.cpp
checkERI_LDADD = $(SERIAL_LIBS)

if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
