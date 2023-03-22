#!/bin/bash

dir_source='/home/sandeep/build_supramolecular_polymers'

export PATH=$dir_source:$PATH 

cp ${dir_source}/dftb_in.hsd .
cp ${dir_source}/input.xyz .
cp ${dir_source}/input.user .
