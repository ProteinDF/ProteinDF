#	RLMO					
#RLMO_SOURCES = \
#	main_RLMO.cpp
#RLMO_LDADD = $(SERIAL_LIBS)

SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@
PARALLEL_LIBS = ../pdf/libpdf.la ../pdf/libppdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@ @FLIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 				\
	pdfdiffmat 				\
	pdfdiffvtr 				\
	topPdfMatrix			\
	catPdfMatrix 			\
	combineDensityMatrix 	\
	pdfMatrixCheck 			\
	mmap_bench 				\
	mmap_bench_matrix_calc 	\
	merge 					\
	transdim				\
	diagonal				\
	cholesky				\
	inverse					\
	getMaxValue				\
	rotateLCAO				\
	\
	pdf-component			\
	pdf-eri					\
	pdf-info-charge 		\
	pdf-info-mat			\
	pdf-info-orb			\
	pdf-init-param			\
	pdf-kwd					\
	pdf-lo					\
	pdf-mat-add				\
	pdf-mat-dot				\
	pdf-mat-ext-block		\
	pdf-mat-mul				\
	pdf-mat-n2s				\
	pdf-mat-resize			\
	pdf-mat-select			\
	pdf-mat-transpose		\
	pdf-mkfld-mo			\
	pdf-mkfld-dens			\
	pdf-mkfld-esp			\
	pdf-pop-mulliken		\
	pdf-show-mat			\
	pdf-show-vtr			\
	pdf-reconf-vectormatrix


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


combineDensityMatrix_SOURCES = main_combineDensityMatrix.cpp
combineDensityMatrix_LDADD = $(SERIAL_LIBS)

pdfMatrixCheck_SOURCES = pdfMatrixCheck.cpp
pdfMatrixCheck_LDADD = $(SERIAL_LIBS)


mmap_bench_SOURCES = mmap_bench.cpp
mmap_bench_LDADD = $(SERIAL_LIBS)


mmap_bench_matrix_calc_SOURCES = mmap_bench_matrix_calc.cpp
mmap_bench_matrix_calc_LDADD = $(SERIAL_LIBS)


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

pdf_info_mat_SOURCES = main_info_mat.cpp
pdf_info_mat_LDADD = $(SERIAL_LIBS)

pdf_info_orb_SOURCES = main_info_orb.cpp
pdf_info_orb_LDADD = $(SERIAL_LIBS)

pdf_init_param_SOURCES = main_init_param.cpp
pdf_init_param_LDADD = $(SERIAL_LIBS)

pdf_kwd_SOURCES = main_kwd.cpp
pdf_kwd_LDADD = $(SERIAL_LIBS)

pdf_lo_SOURCES = \
	DfLocalize.h	DfLocalize.cpp 	\
	main_localize.cpp
pdf_lo_LDADD = $(SERIAL_LIBS)

pdf_mat_add_SOURCES = main_mat_add.cpp
pdf_mat_add_LDADD = $(SERIAL_LIBS)

pdf_mat_dot_SOURCES = main_mat_dot.cpp
pdf_mat_dot_LDADD = $(SERIAL_LIBS)

pdf_mat_ext_block_SOURCES = main_mat_ext_block.cpp
pdf_mat_ext_block_LDADD = $(SERIAL_LIBS)

pdf_mat_mul_SOURCES = main_mat_mul.cpp
pdf_mat_mul_LDADD = $(SERIAL_LIBS)

pdf_mat_n2s_SOURCES = main_mat_n2s.cpp
pdf_mat_n2s_LDADD = $(SERIAL_LIBS)

pdf_mat_resize_SOURCES = main_mat_resize.cpp
pdf_mat_resize_LDADD = $(SERIAL_LIBS)

pdf_mat_select_SOURCES = main_mat_select.cpp
pdf_mat_select_LDADD = $(SERIAL_LIBS)

pdf_mat_transpose_SOURCES = main_mat_transpose.cpp
pdf_mat_transpose_LDADD = $(SERIAL_LIBS)

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

pdf_reconf_vectormatrix_SOURCES = main_reconf_vectormatrix.cpp
pdf_reconf_vectormatrix_LDADD = $(SERIAL_LIBS)

if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
