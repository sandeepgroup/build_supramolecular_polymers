#!/bin/bash 
# this script combines all energy json files and creates the final energy file 

ip_ener='energy.json'
output='combined_energy.json'

file=''
for i in `ls -d ener-* | sort -V`; do 
  file="$i/${ip_ener} $file"
done 

jq -s 'reduce .[] as $item ({}; . * $item)' $file > $output

echo " LOG: all energy.json files are combined. Output file $output" 

# check for entries having values 0 

mapfile -t values < <(jq -r '.[]' $output)
len=`jq 'length' ${output}`
len=$((len-1)) 

count=0
for j in $(seq 0 $len); do      
  value=${values[$j]}
  check=`echo $value == 0 | bc`
  if [ $check -eq 1 ]; then
    count=$((count+1)) 
  fi 
done 

if [ $count -gt 0 ]; then 
  echo " WARNING: $count entries out of $((len+1)) are zero. Check!" 
fi 

