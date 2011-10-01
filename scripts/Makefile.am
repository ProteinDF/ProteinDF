bin_script_src = \
	pdf_common.sh		\
	pdf_setup.sh		\
	pdf_serial.sh		\
	pdf_parallel.sh		\
	PDF.sh			\
	PPDF.sh			\
	\
	pdf_mkdir.sh		\
	pdf_clean.sh		\
	print_matrix.py		\
	print_vector.py		\
	\
	mpac2yml.py		\
	yml2mpac.py		


bin_SCRIPTS = $(bin_script_src)
dist_bin_SCRIPTS = $(bin_script_src)


install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for PROG in $(bin_script_src); do \
		echo $${PROG%.*}; \
		$(LN_S) $${PROG} $${PROG%.*}; \
	done

