#	RLMO					
#RLMO_SOURCES = \
#	main_RLMO.cpp
#RLMO_LDADD = $(SERIAL_LIBS)

SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@
PARALLEL_LIBS = ../pdf/libpdf.la ../pdf/libppdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 				\
	pdfdiffmat 				\
	pdfdiffvtr 				\
	topPdfMatrix				\
	catPdfMatrix addPdfMatrix	\
	combineDensityMatrix 		\
	pdfMatrixCheck 			\
	mmap_bench 				\
	mmap_bench_matrix_calc 	\
	localize 				\
	merge 					\
	transdim				\
	diagonal				\
	cholesky				\
	inverse					\
	multiple				\
	transpose				\
	getMaxValue				\
	rotateLCAO				\
	\
	pdf-component			\
	pdf-eri					\
	pdf-info-charge 		\
	pdf-info-orb			\
	pdf-kwd					\
	pdf-mkfld-mo			\
	pdf-mkfld-dens			\
	pdf-mkfld-esp			\
	pdf-pop-mulliken		\
	pdf-show-mat			\
	pdf-show-vtr


PARALLEL_PROGS = \
	pdiagonal 				\
	plocalize 				\
	mergeLocalMatrix 		\
	divideGlobalMatrix 		\
	pdfdivmat

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

pdfMatrixCheck_SOURCES = pdfMatrixCheck.cpp
pdfMatrixCheck_LDADD = $(SERIAL_LIBS)


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

pdfdivmat_SOURCES = main_pdfdivmat.cpp
pdfdivmat_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)


pdf_component_SOURCES = main_component.cpp
pdf_component_LDADD = $(SERIAL_LIBS)

pdf_eri_SOURCES = main_eri.cpp
pdf_eri_LDADD = $(SERIAL_LIBS)

pdf_info_charge_SOURCES = main_info_charge.cpp
pdf_info_charge_LDADD = $(SERIAL_LIBS)

pdf_info_orb_SOURCES = main_info_orb.cpp
pdf_info_orb_LDADD = $(SERIAL_LIBS)

pdf_kwd_SOURCES = main_kwd.cpp
pdf_kwd_LDADD = $(SERIAL_LIBS)

pdf_mkfld_mo_SOURCES = \
	mkfld_common.h	mkfld_common.cpp	\
	TlMoField.h		TlMoField.cpp		\
	main_mkfld_mo.cpp
pdf_mkfld_mo_LDADD = $(SERIAL_LIBS)

pdf_mkfld_dens_SOURCES = \
	mkfld_common.h	mkfld_common.cpp	\
	TlDensField.h	TlDensField.cpp		\
	main_mkfld_dens.cpp
pdf_mkfld_dens_LDADD = $(SERIAL_LIBS)

pdf_mkfld_esp_SOURCES = \
	mkfld_common.h	mkfld_common.cpp	\
	TlEspField.h	TlEspField.cpp		\
	main_mkfld_esp.cpp
pdf_mkfld_esp_LDADD = $(SERIAL_LIBS)

pdf_pop_mulliken_SOURCES = main_pop_mulliken.cpp
pdf_pop_mulliken_LDADD = $(SERIAL_LIBS)

pdf_show_mat_SOURCES = main_show_mat.cpp
pdf_show_mat_LDADD = $(SERIAL_LIBS)

pdf_show_vtr_SOURCES = main_show_vtr.cpp
pdf_show_vtr_LDADD = $(SERIAL_LIBS)



if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
