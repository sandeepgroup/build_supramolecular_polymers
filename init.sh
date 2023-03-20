#!/bin/bash

dir_export='/home/sandeep/build_supramolecular_polymers'

export PATH=$dir_export:$PATH 

bash -c ". ${dir_export}" 

cp ${dir_export}/dftb_in.hsd .
cp ${dir_export}/input.xyz .
cp ${dir_export}/user_input.user .
