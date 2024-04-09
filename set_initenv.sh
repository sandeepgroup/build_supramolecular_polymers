#!/bin/bash
# set the environment and creates the commands
# stackgen_run, stackgen_clean
# 

STACKGEN_SRC='/home/arunima/stackgen_updated'


# check if STACKGEN source directory is set properly or not; 

if [ ! -f "${STACKGEN_SRC}/main.py" ]; then
	echo 'STACKGEN source directory is not set to the correct location'
	echo 'modify this script to set the path'
	echo 'exiting ...'
	return
fi 

# set DFTB_PREFIX variable to the directory where the Slater-Koster files are kept
export DFTB_PREFIX="/home/sandeep/software/dftb-input/3ob-3-1/"
# setting up the environment 

if [ -n "${STACKGEN_SRC}" ]; then 
	export PATH=${STACKGEN_SRC}:$PATH 
	export PATH=${STACKGEN_SRC}/utils:$PATH 
fi

type stackgen_run > /dev/null  2>&1
if [ $? -ne 0 ]; then
	alias stackgen_run="python -u ${STACKGEN_SRC}/main.py"
fi

type stackgen_clean > /dev/null  2>&1
if [ $? -ne 0 ]; then 
	alias stackgen_clean="bash ${STACKGEN_SRC}/utils/clean.sh"
fi



