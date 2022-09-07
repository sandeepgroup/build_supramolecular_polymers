#!/bin/bash
# uses jq and mapfile commands. Install them if not available

# this script will check the presence of all input files needed and then 
# check dftb+ calculation on a test system
# go through the trajectory file and for each configuration, it calculates energy
# and updates the energy.json file 
# this script also works for restarting 

#energies are in kJ/mol 


if [ ! -f "../input-3.user" ]; then
echo " file ../input-3.user not found"
echo " stop ... "
exit 0 
fi

source ../input-3.user

export OMP_NUM_THREADS=$np

exe='dftb+' 
output='dftb.out'
dftb_input='dftb_in.hsd'
conv_factor=2625.5 

# check whether input env is set up properly or not 
which $exe > /dev/null
check=$? 
if [ ! -f "../${dftb_input}" ]; then
  echo " ERROR: dftb+ input file \"${dftb_input}\" is not found"
  echo " ERROR: stop ... "
  exit 0 
elif [ $check -ne 0 ]; then
  echo " ERROR: dftb+ executable not on the path "
  echo " ERROR: stop ... "
  exit 0 
elif [ ! -f "${ip_energy}" ]; then
  echo  " ERROR: file \"${ip_energy}\" not found "
  echo " ERROR: stop ... "
  exit 0 
elif [ ! -f "${traj}" ]; then
  echo " ERROR:file \"${traj}\" not found "
  echo " ERROR: stop ... "
  exit 0 
else 
 # check whether dftb running fine by running water molecule scf calc 

 sed -e '/ N /d' -e '/ C /d' ../${dftb_input} > ${dftb_input}
 tmp='input.xyz' 
 echo "3" > $tmp
 echo  >> $tmp
 echo "O  -0.23980816   -1.09112708    0.00000000" >> $tmp
 echo "H   0.72019184   -1.09112708    0.00000000" >> $tmp 
 echo "H  -0.56026274   -0.18619125    0.00000000" >> $tmp 

 $exe > $output 
 pe=`grep 'Total Energy: ' $output | awk '{print $3}'`

 if [ -z $pe ]; then
  echo " ERROR: dftb+ calculation not successful. Check dftb+ input file"
  echo " ERROR: stop ... "
  exit 0 
 else
  echo " LOG: dftb calculation successful on test (h2o) system. continue..."
 fi 

 #clean
 rm -f input.xyz 
 rm -f charges.bin
 rm -f dftb_pin.hsd
 rm -f tmp.energy
 rm -f ${dftb_input} 
 rm -f $output
fi 

natoms=`head -n 1 $traj`
nlines=`wc -l $traj | awk '{print $1}'`
nrep=$((nlines/(natoms+2)))

# collect all keys from the json file 

mapfile -t keys < <(jq -r 'keys_unsorted[]' $ip_energy)
mapfile -t values < <(jq -r '.[]' $ip_energy)

# will check energy.json file. If energy is zero, it will start the dftb+ calculation. 
# helpful for restarting 
# also, it creates a back of energy.json file for every 30 configurations, by default 

len=`jq 'length' ${ip_energy}`
len=$((len-1)) 

ln -fs ../${dftb_input} . 
count=0 ; count_res=0 
nameres=${ip_energy}
for j in $(seq 0 $len); do      

 count=$((count+1))  
 check=`echo ${values[$j]}==0.0| bc`
 if [ ${check} == 1  ]; then          # when energy is 0.0 
  
   key=${keys[$j]}

   grep -C $((natoms)) "$key" $traj | tail -n $((natoms+2)) > input.xyz 
   $exe > $output 
   pe=`grep 'Total Energy: ' $output | awk '{print $3}'`

# update pe in json file 
   if [ ! -z $pe ]; then 
     mv ${ip_energy} tmp.energy 
     pe=`echo $pe*${conv_factor} | bc -l`         # kJ/mol 
     jq -r '."'${key}'" |= "'${pe}'"' tmp.energy > ${ip_energy}
     printf "%10s%7d%1s%-7d%-12s\n" " LOG: conf" $count "/" $nrep " done"

     count_res=$((count_res+1))

     # creates backup file of energy.json file for every 30 conf
     if [ $count_res == 30 ]; then 
       i=1
       while [ -e \#$nameres.$i\# ] ; do 
         i=$((i+1)) 
       done 
       cp $nameres \#$nameres.$i\#
       count_res=0 
     fi 

   else
     printf "%10s%7d%1s%-7d%-12s\n" " LOG: conf" $count "/" $nrep " failed"
   fi
   #clean
   rm -f input.xyz 
   rm -f charges.bin
   rm -f dftb_pin.hsd
   rm -f tmp.energy
   rm -f $output

#   sed -i '2s/^.*$//' input.xyz

  else 
    printf "%10s%7d%1s%-7d%12s\n" " LOG: conf" $count "/" $nrep " already done. SKIPPING"
 fi 
done 

rm -f ${dftb_input} 


