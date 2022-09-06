#!/usr/bin/env python3
# it requires two environment variables

import os
import sys
import re
import glob
import json


print(" Log: Using the ", glob.glob("traj*.xyz"), "as the input file for this script")

filename=glob.glob("traj*.xyz")[0]

output="split_"+filename
output_energy='energy.json'

nline = len(open(filename).readlines())
natoms = int(open(filename).readlines()[0])
nrep=int(nline/(natoms+2))

# read input.user file 
ip={}
ip["nconf_calc"]=1
ip["prefix"]="energy_dftb"

with open('input-2.user','r') as fp:
  for line in fp:
    if not len(line.strip())==0 and not line.startswith("#"):
      name,var = line.partition('=')[::2]
      var=var.strip()
      if re.match("^[0-9-+]*$", var):
        var=int(var)
      elif re.match("^[0-9-+.]*$", var):   
        var=float(var)

      ip[name.strip()]=var

fp.close()

# create a directory prefix_$num,... and start dividing the trajectory
prefix=ip["prefix"]
count=0; cdir=100
fp=open(filename,"r")

for rep in range(nrep): 
  if count == 0:

    cdir+=1 
    dict={}
    dir=prefix+"_"+str(cdir) 
    if not os.path.exists(dir):
      os.mkdir(dir)
    else:
      print(" Warning: directory "+ dir +" already exists  ") 

    fo=open("%s_%s/%s" %(prefix,cdir,output),'w')
    count=ip["nconf_calc"]
  
  fp.readline()
  line=re.split("\s+|_",fp.readline())[2:-1]
  key=''
  for element in line[:-1]: 
    key=key+element+"_"

  key=key+line[-1]
  dict[key]=0.0

  count-=1 

  fo.write(str(natoms)+"\n")
  fo.write(str(key)+"\n")

  for atoms in range(int(natoms)): 
    line = fp.readline().split()
    fo.write("%s " %(line[0])+" ")
    for i in range(1,len(line)):
      fo.write("%-7.3f" %(float(line[i]))+" ")

    fo.write("\n")

# write key and energy file in json format 

  if count ==0:
    fener=open("%s_%s/%s" %(prefix,cdir,output_energy),'w')
    x = json.dumps(dict,indent=4) 
    fener.write(x)
    fener.close()  

  
print(" LOG: coordinates are written to " + output)
print(" LOG: energy json init is written to " + output_energy)
print(" LOG: No. of directories created = " + str(ip["nconf_calc"])) 
print(" LOG: Successful ")


 


