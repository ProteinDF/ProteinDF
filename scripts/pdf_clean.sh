#!/bin/sh

RM="rm"

while getopts fk opt; do
case ${opt} in
    f)
	    # fully clean up
	    FULLY_CLEAN_UP=true;;
    k)
	    # keep integrals
	    KEEP_INTEGRALS=true;;
    *)
	    echo "unkown option: ${opt}"
	    exit 1;;
esac
done

for list in fl_Input fl_Table fl_Temp; do
    if [ -d ${list} ]; then
	find ${list} -maxdepth 1 -name "*" -and -type f -print0 | xargs -0 --no-run-if-empty ${RM}
    fi
done

if [ -d fl_Work ]; then
    if [ x${KEEP_INTEGRALS} = xtrue ]; then
	find fl_Work -regex ".*[0-9]+" -print | xargs -0 --no-run-if-empty rm
    else
	find fl_Work -maxdepth 1 -name "*" -and -type f -print0 | xargs -0 --no-run-if-empty ${RM}
    fi
fi

for target in fl_Out_Std fl_Out_Arc fl_Out_Sys PDF.pid fl_Plot; do
    if [ -f $target ]; then
	${RM} $target
    fi
done

if [ x${FULLY_CLEAN_UP} = xtrue ]; then
    for target in pdfparam.mpac; do
        if [ -f $target ]; then
            ${RM} $target
	fi
    done
fi
