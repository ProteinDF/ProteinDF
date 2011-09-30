bin_script_src = \
	pdf_common.sh		\
	pdf_setup.sh		\
	pdf_serial.sh		\
	pdf_parallel.sh		\
	PDF.sh			\
	PPDF.sh			\
	\
	diff_matrix.py		\
	pdf_stat.py		\
	print_matrix.py		\
	print_vector.py		\
	\
	pdf_mkdir.sh		\
	pdf_clean.sh		\
	\
	pdf_xml.py		\
	pdf_xml_diff.py		\
	pdfdata.py		\
	PdfOutput.py		\
	PdfXml.py		\
	\
	pdf_xml14.py		\
	\
	brposition.py		\
	bratom.py		\
	bratomgroup.py		\
	brperiodictable.py	\
	brpdb.py		\
	\
	dfdata.py		\
	dfoutput.py		\
	dfvector.py		\
	dfmatrix.py		\
	viewmpac.py		\
	mpac2yml.py		\
	yml2mpac.py		\
	tlgraph.py		\
	makepack.py		\
	makegraph.py		\
	reg_harris.py		\
	pdf_extract.py		\
	pdf_report.py		\
	pdf_diff_results.py


bin_SCRIPTS = $(bin_script_src)
dist_bin_SCRIPTS = $(bin_script_src)


install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for PROG in $(bin_script_src); do \
		echo $${PROG%.*}; \
		$(LN_S) $${PROG} $${PROG%.*}; \
	done

