#!/bin/bash
# uses jq and mapfile commands. Install them if not available

# this script will check the presence of all input files needed and then 
# check dftb+ calculation on a test system
# go through the trajectory file and for each configuration, it calculates energy
# and updates the energy.json file 
# this script also works for restarting 

#energies are in kJ/mol 
function dftb_energy_cal {

	source user_input.user

	export OMP_NUM_THREADS=${nproc_dftb}

	exe='dftb+' 
	output='detailed.out'
	#output='dftb.out'
	dftb_input='dftb_in.hsd'
	conv_factor=2625.5 

	ln -s $1  tmpinput1.xyz  
	#ln -s generated_*_*.xyz  input1.xyz  

	$exe 2&>/dev/null 
	pe=`grep 'Total energy: ' $output | awk '{print $3}'`
	if [ ! -z $pe ]; then #the energy is not equal to 0
	        pe=`echo $pe*${conv_factor} | bc -l` 
		local energy_val=$(echo "$pe" | bc -l)
		echo $energy_val
	else
		echo " ERROR: Some issue with dftb+ calculation " 
	#       echo "failed"
	#       en=`awk '{print $2}' dftb.out | grep '[[:digit:]]' | tail -n 1`
	#       en=`echo $en | awk -v x=$conv_factor '{printf "%0.5f\n",$1*x}'`
	#       local energy_val=$(echo "$en" | bc -l)
	#       echo $energy_val
   	fi
	rm -f tmpinput1.xyz 
	rm -f charges.bin
	rm -f dftb_pin.hsd
	rm -f tmp.energy
	rm -f $output
	#local energy_val=$(echo "$pe" | bc -l)
	#echo $energy_val
	#return "$energy_val"
}

input_file=$1 
energy_value=$(dftb_energy_cal "$input_file") 
echo $energy_value


