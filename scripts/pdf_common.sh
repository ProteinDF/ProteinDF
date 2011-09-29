#!/bin/sh

# parameter
RM="rm"
PDF_WORK_DIRS="fl_Input fl_Table fl_Temp fl_Work"

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


# ProteinDFを実行する前に必要なディレクトリを作成する
prepare_pdf()
{
    TARGET_DIR=.
    if [ x${1} != x ]; then
	TARGET_DIR=${1}
    fi
	
    for d in ${PDF_WORK_DIRS}; do
	TARGET=${TARGET_DIR}/${d}
	if [ ! -d ${TARGET} ]; then
	    mkdir -p ${TARGET}
	fi
    done
}


# ProteinDFで用いた作業ファイルを削除する
pdf_clean()
{
    
    for DIR in ${PDF_WORK_DIRS}; do
	if [ x${KEEP_INTEGRALS} = xtrue ]; then
	    find ${DIR} -regex ".*[0-9]+" -print | xargs -0 --no-run-if-empty ${RM}
	else
	    find ${DIR} -maxdepth 1 -name "*" -and -type f -print0 | xargs -0 --no-run-if-empty ${RM}
	fi
    done

    for FILE in fl_Out_Std pdfparam.mpac; do
	if [ -f ${FILE} ]; then
	    ${RM} ${FILE}
	fi
    done
}

