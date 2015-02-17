SERIAL_LIBS = ../pdf/libpdf.la ../pdflib/libpdftl.la @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@
PARALLEL_LIBS = ../pdf/libpdf.la ../pdf/libppdf.la ../pdflib/libppdftl.la ../pdflib/libpdftl.la @SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@ @FLIBS@

AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/src/pdf

SERIAL_PROGS = 				\
	pdf-avs_merge 			\
	pdf-bench-mmap			\
	pdf-bench-mmap-mat		\
	pdf-component			\
	pdf-eri					\
	pdf-info-charge 		\
	pdf-info-grid			\
	pdf-info-orb			\
	pdf-init-param			\
	pdf-kwd					\
	pdf-lo					\
	pdf-mat-add				\
	pdf-mat-cholesky		\
	pdf-mat-dot				\
	pdf-mat-diagonal		\
	pdf-mat-diff			\
	pdf-mat-ext				\
	pdf-mat-format			\
	pdf-mat-info			\
	pdf-mat-inverse			\
	pdf-mat-mul				\
	pdf-mat-symmetrize		\
	pdf-mat-resize			\
	pdf-mat-select			\
	pdf-mat-show			\
	pdf-mat-sum				\
	pdf-mat-top				\
	pdf-i2pq-show			\
	pdf-mat-transpose		\
	pdf-mkfld-mo			\
	pdf-mkfld-dens			\
	pdf-mkfld-esp			\
	pdf-pop-mulliken		\
	pdf-rotate-lcao			\
	pdf-trans-densmat		\
	pdf-reconf-vectormatrix	\
	pdf-vtr-diff			\
	pdf-vtr-show			


PARALLEL_PROGS = \
	pdiagonal 				\
	plocalize 				\
	mergeLocalMatrix 		\
	divideGlobalMatrix 		\
	pdfdivmat


plocalize_SOURCES = \
	DfLocalize.h		DfLocalize.cpp	\
	DfLocalize_Parallel.h	DfLocalize_Parallel.cpp \
	main_plocalize.cpp
plocalize_LDADD = $(PARALLEL_LIBS)


pdiagonal_SOURCES = main_pdiagonal.cpp
pdiagonal_LDADD = $(PARALLEL_LIBS)


mergeLocalMatrix_SOURCES = main_mergeLocalMatrix.cpp
mergeLocalMatrix_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)

divideGlobalMatrix_SOURCES = main_divideGlobalMatrix.cpp
divideGlobalMatrix_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)

pdfdivmat_SOURCES = main_pdfdivmat.cpp
pdfdivmat_LDADD = $(SERIAL_LIBS) $(PARALLEL_LIBS)

pdf_avs_merge_SOURCES = main_avs_merge.cpp
pdf_avs_merge_LDADD = $(SERIAL_LIBS)

pdf_bench_mmap_SOURCES = main_bench_mmap.cpp
pdf_bench_mmap_LDADD = $(SERIAL_LIBS)

pdf_bench_mmap_mat_SOURCES = main_bench_mmap_mat.cpp
pdf_bench_mmap_mat_LDADD = $(SERIAL_LIBS)

pdf_component_SOURCES = main_component.cpp
pdf_component_LDADD = $(SERIAL_LIBS)

pdf_eri_SOURCES = main_eri.cpp
pdf_eri_LDADD = $(SERIAL_LIBS)

pdf_info_charge_SOURCES = main_info_charge.cpp
pdf_info_charge_LDADD = $(SERIAL_LIBS)

pdf_info_grid_SOURCES = main_info_grid.cpp
pdf_info_grid_LDADD = $(SERIAL_LIBS)

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

pdf_mat_cholesky_SOURCES = main_mat_cholesky.cpp
pdf_mat_cholesky_LDADD = $(SERIAL_LIBS)

pdf_mat_dot_SOURCES = main_mat_dot.cpp
pdf_mat_dot_LDADD = $(SERIAL_LIBS)

pdf_mat_diagonal_SOURCES = main_diagonal.cpp
pdf_mat_diagonal_LDADD = $(SERIAL_LIBS)

pdf_mat_diff_SOURCES = main_mat_diff.cpp
pdf_mat_diff_LDADD = $(SERIAL_LIBS)

pdf_mat_ext_SOURCES = main_mat_ext.cpp
pdf_mat_ext_LDADD = $(SERIAL_LIBS)

pdf_mat_format_SOURCES = main_mat_format.cpp
pdf_mat_format_LDADD = $(SERIAL_LIBS)

pdf_mat_info_SOURCES = main_mat_info.cpp
pdf_mat_info_LDADD = $(SERIAL_LIBS)

pdf_mat_inverse_SOURCES = main_mat_inverse.cpp
pdf_mat_inverse_LDADD = $(SERIAL_LIBS)

pdf_mat_mul_SOURCES = main_mat_mul.cpp
pdf_mat_mul_LDADD = $(SERIAL_LIBS)

pdf_mat_symmetrize_SOURCES = main_mat_symmetrize.cpp
pdf_mat_symmetrize_LDADD = $(SERIAL_LIBS)

pdf_mat_show_SOURCES = main_mat_show.cpp
pdf_mat_show_LDADD = $(SERIAL_LIBS)

pdf_mat_sum_SOURCES = main_mat_sum.cpp
pdf_mat_sum_LDADD = $(SERIAL_LIBS)

pdf_mat_resize_SOURCES = main_mat_resize.cpp
pdf_mat_resize_LDADD = $(SERIAL_LIBS)

pdf_mat_select_SOURCES = main_mat_select.cpp
pdf_mat_select_LDADD = $(SERIAL_LIBS)

pdf_mat_top_SOURCES = main_mat_top.cpp
pdf_mat_top_LDADD = $(SERIAL_LIBS)

pdf_mat_transpose_SOURCES = main_mat_transpose.cpp
pdf_mat_transpose_LDADD = $(SERIAL_LIBS)

pdf_i2pq_show_SOURCES = main_i2pq_show.cpp
pdf_i2pq_show_LDADD = $(SERIAL_LIBS)

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

pdf_rotate_lcao_SOURCES = \
	TlRotateLCAO.h		\
	TlRotateLCAO.cpp	\
	main_rotate_lcao.cpp
pdf_rotate_lcao_LDADD = $(SERIAL_LIBS)

pdf_reconf_vectormatrix_SOURCES = main_reconf_vectormatrix.cpp
pdf_reconf_vectormatrix_LDADD = $(SERIAL_LIBS)

pdf_trans_densmat_SOURCES = main_trans_densmat.cpp
pdf_trans_densmat_LDADD = $(SERIAL_LIBS)

pdf_vtr_diff_SOURCES = main_vtr_diff.cpp
pdf_vtr_diff_LDADD = $(SERIAL_LIBS)

pdf_vtr_show_SOURCES = main_vtr_show.cpp
pdf_vtr_show_LDADD = $(SERIAL_LIBS)

if PARALLEL
bin_PROGRAMS = $(SERIAL_PROGS) $(PARALLEL_PROGS)
else
bin_PROGRAMS = $(SERIAL_PROGS)
EXTRA_PROGRAMS = $(PARALLEL_PROGS)
endif
