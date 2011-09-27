
TARGET = \
	01-DiPep 			\
	02-Met				\
	03-Asp				\
	04-1HRC.20_22			\
	05-Gly.rks.rho			\
	06-Gly.rks.lcao			\
	07-Gly.uks.rho			\
	08-Gly.uks.lcao			\
	\
	10-Gly.rks.rho.ls_closed	\
	11-Gly+.uks.rho.ls_virtual	\
	12-Gly+.uks.rho.ls_open		\
	13-Gly+.uks.lcao		\
	14-Gly+.uks.lcao.ls_open	\
	15-CoF6				\
	16-H2+.uks.rho			\
	17-H2+.uks.lcao

.PHONY: all clean serial parallel 

all: serial

$(TARGET): 
	@make $($@:%=%.log)

#%.log:
#	(cd $*; $(PDF_HOME)/bin/pdf_mkdir.sh; make -f ../Makefile.common parallel)
#	touch $@

serial: $(TARGET:%=%.serial)

%.serial:
	(cd $*; $(PDF_HOME)/bin/pdf_mkdir.sh; make -f ../Makefile.common serial) 
	@touch $@
	@echo

parallel: $(TARGET:%=%.parallel)

%.parallel:
	(cd $*; $(PDF_HOME)/bin/pdf_mkdir.sh; make -f ../Makefile.common parallel) 
	@touch $@
	@echo

clean:
	#@$(RM) $(TARGET:%=%.log)
	@$(RM) $(TARGET:%=%.serial)
	@$(RM) $(TARGET:%=%.parallel)

	@for i in $(TARGET); do\
	(cd $$i; make -f ../Makefile.common clean;); \
	done


