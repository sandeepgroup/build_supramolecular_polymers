#!/bin/bash
# uses jq and mapfile commands. Install them if not available

source ../input-3.user

export OMP_NUM_THREADS=$np

exe='dftb+' 
output='dftb.out'
dft_input='dftb_in.hsd'
natoms=`head -n 1 $traj`
nlines=`wc -l $traj | awk '{print $1}'`
nrep=$((nlines/(natoms+2)))


which $exe > /dev/null
check=$? 

if [ ! -f "../${dft_input}" ]; then
  echo dftb+ input file \"${dft_input}\" is not found
  echo stop ... 
  exit 0 
elif [ $check -ne 0 ]; then
  echo dftb+ executable not on the path 
  echo stop ... 
  exit 0 
fi 

ln -fs ../${dft_input} . 

# collect all keys from the json file 

mapfile -t keys < <(jq -r 'keys_unsorted[]' $ip_energy)
mapfile -t values < <(jq -r '.[]' $ip_energy)

# will check energy.json file. If energy is zero, it will start the dftb+ calculation. 
# helpful for restarting 

len=`jq 'length' ${ip_energy}`
len=$((len-1)) 

count=0 
for j in $(seq 0 $len); do      

 count=$((count+1))  
 check=`echo ${values[$j]}==0.0| bc`
 if [ ${check} == 1  ]; then          # when energy is 0.0 
  
   key=${keys[$j]}

   grep -C $((natoms)) "$key" $traj | tail -n $((natoms+2)) > input.xyz 
   $exe > $output 
   pe=`grep 'Total Energy: ' $output | awk '{print $3}'`

# update pe in json file 

   mv ${ip_energy} tmp.energy 
   jq -r '."'${key}'" |= "'${pe}'"' tmp.energy > ${ip_energy}

   printf "%10s%7d%1s%-7d%-12s\n" " LOG: conf" $count "/" $nrep " done"
   #clean
   rm -f input.xyz 
   rm -f charges.bin
   rm -f dftb_pin.hsd
   rm -f tmp.energy

#   sed -i '2s/^.*$//' input.xyz

  else 
    printf "%10s%7d%1s%-7d%12s\n" " LOG: conf" $count "/" $nrep " already done. SKIPPING"
 fi 
done 

rm -f ${dft_input} 


