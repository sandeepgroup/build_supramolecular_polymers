#!/bin/bash
# this script will check the presence of all input files needed and then 
# check dftb+ calculation on a test system

export OMP_NUM_THREADS=1

exe='dftb+' 
output='dftb.out'
dftb_input='dftb_in.hsd'
conv_factor=2625.5 

# check whether input env is set up properly or not 
if [ ! -f "${dftb_input}" ]; then
  	exit 1
fi

which $exe > /dev/null
check=$? 
if [ $check -ne 0 ]; then
  	exit 2
fi

#check whether dftb running fine by running water molecule scf calc 
#mkdir water_check && touch water_check/input.hsd

export OMP_NUM_THREADS=1

scp ${dftb_input} ${dftb_input}.bak
sed -e '/ N /d' -e '/ C /d' ${dftb_input} > tmpdftb
mv tmpdftb ${dftb_input} 
tmp='input1.xyz' 
echo "3" > $tmp
echo  >> $tmp
echo "O  -0.23980816   -1.09112708    0.00000000" >> $tmp
echo "H   0.72019184   -1.09112708    0.00000000" >> $tmp 
echo "H  -0.56026274   -0.18619125    0.00000000" >> $tmp 
$exe > $output 
pe=`grep 'Total Energy: ' $output | awk '{print $3}'`

#clean
mv -f ${dftb_input}.bak ${dftb_input}
rm -f input1.xyz 
rm -f charges.bin
rm -f dftb_pin.hsd
rm -f tmp.energy
rm -f $output
rm -f $tmp

if [ -z $pe ]; then
        exit 3
	echo " ERROR: dftb+ calculation not successful. Check dftb+ input file"
	echo " ERROR: stop ... "
else
	echo " LOG: dftb calculation successful on test (h2o) system. continue..."
fi 


