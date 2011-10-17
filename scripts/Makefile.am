bin_script_src = \
	pdfcommon.sh		\
	pdfsetup.sh		\
	PDF.sh			\
	PPDF.sh			\
	pdfclean.sh


bin_SCRIPTS = $(bin_script_src)
dist_bin_SCRIPTS = $(bin_script_src)


install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for PROG in $(bin_script_src); do \
		echo $${PROG%.*}; \
		$(LN_S) -f $${PROG} $${PROG%.*}; \
	done

