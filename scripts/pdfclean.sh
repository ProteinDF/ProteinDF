#!/bin/sh

check_pdfhome()
{
    if [ x${PDF_HOME} = x ]; then
	echo "PDF_HOME environment variable is not set."
	echo "please set the parameter."
	return 1
    else
	return 0
    fi
}

# main
check_pdfhome
. ${PDF_HOME}/bin/pdf_common.sh

FULL_CLEANUP=false
KEEP_INTEGRALS=true

if [ x${1} = xall ]; then
    FULL_CLEANUP=true
    KEEP_INTEGRALS=false
fi

# cleanup directories
for DIR in ${PDF_WORK_DIRS}; do
    if [ -d ${DIR} ]; then
	if [ x${KEEP_INTEGRALS} = xtrue ]; then
	    find ${DIR} -regex ".*[0-9]+" -and -type f -print0 | xargs -0 --no-run-if-empty ${RM}
	else
	    find ${DIR} -maxdepth 1 -name "*" -and -type f -print0 | xargs -0 --no-run-if-empty ${RM}
	fi
    fi
done

# cleanup files
PDF_OUTPUT_FILES="fl_Out_Std pdfparam.mpac"
if [ x${FULL_CLEANUP} = xtrue ]; then
    for FILE in ${PDF_OUTPUT_FILES}; do
	if [ -f $FILE ]; then
	    ${RM} $FILE
	fi
    done
fi

