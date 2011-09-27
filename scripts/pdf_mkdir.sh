#/bin/sh +x

PDF_WORK_DIR="fl_Input fl_Table fl_Temp fl_Work"

TARGET_DIR=.
if [ x${1} != x ]; then
	TARGET_DIR=${1}
fi
	
for i in ${PDF_WORK_DIR}; do
	TARGET=${TARGET_DIR}/${i}
	if [ ! -d ${TARGET} ]; then
		mkdir -p ${TARGET}
	fi
done



