#!/bin/sh

RM="rm"

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

PDF_WORK_DIRS="fl_Work fl_Input fl_Table fl_Temp"

# main
check_pdfhome

FULL_CLEANUP=false
KEEP_INTEGRALS=true

if [ x${1} = xall ]; then
    FULL_CLEANUP=true
    KEEP_INTEGRALS=false
fi

# cleanup directories
XARGS_NO_RUN_IF_EMPTY=''
if [ `uname` = Linux ]; then
    XARGS_NO_RUN_IF_EMPTY='--no-run-if-empty'
fi
for DIR in ${PDF_WORK_DIRS}; do
    if [ -d ${DIR} ]; then
	    if [ x${KEEP_INTEGRALS} = xtrue ]; then
	        find ${DIR} -regex ".*[0-9]+" -and -type f -print0 | xargs -0 ${XARGS_NO_RUN_IF_EMPTY} ${RM}
	    else
	        find ${DIR} -maxdepth 1 -name "*" -and -type f -print0 | xargs -0 ${XARGS_NO_RUN_IF_EMPTY} ${RM}
	    fi
    fi
done

# cleanup files
PDF_OUTPUT_FILES="fl_Out_Std pdfparam.mpac pdfresults.db"
if [ x${FULL_CLEANUP} = xtrue ]; then
    for FILE in ${PDF_OUTPUT_FILES}; do
	if [ -f $FILE ]; then
	    ${RM} $FILE
	fi
    done
fi

