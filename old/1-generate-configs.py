#!/usr/bin/env python3
# code to generate trajectory of configurations for a given order paramters
# translation along three axis and the twist angle are the order parameters 
# it also creates a energy.json file for later use 

# the code creates tmpconfig* files. Hence, do not keep any files with prefix tmpconfig.
# if they are present, they are deleted. 

import os
import sys
import re 
import glob 
import re
import json
from natsort import natsorted

# read input.user file and store the data as dictionary 
#default values

ip={}

ip["label"]="supra"
ip["tx"]=0.0
ip["ty"]=0.0
ip["tz"]=3.5
ip["twist"]=30.0
ip["size"]=2
ip["dtx"]=0.0
ip["dty"]=0.0
ip["dtz"]=0.0
ip["dtwist"]=0.0
ip["ntx"]=1
ip["nty"]=1
ip["ntz"]=1
ip["ntwist"]=1


with open('input-1.user','r') as fp:
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

# get natom, tot_atoms, total no. configurations 

with open(ip["filename"],'r') as fp:
  natom = fp.readline().rstrip()
fp.close()

totatom=int(natom)*ip["size"]
totconf=ip["ntx"]*ip["nty"]*ip["ntz"]*ip["ntwist"]

output_energy='energy.json'
output='traj_'+ip["label"]+'_'+str(ip["size"])+'.xyz'

# clean up any tmpconfig* and output files present in the current working directory

os.system("rm -f tmpconfig*")
os.system("rm -f "+ output)

# re-orient the given configuration so that atom1-atom2 is along x-axis and
# the plane normal is along z-axis 

tmpconfig='tmpconfig'
input_in="python3 ./orient.py "+ ip["filename"] + " -p " + str(ip["atom1"]) +" "+ str(ip["atom2"]) +" " + str(ip["atom3"]) + " -tc > " + tmpconfig +"0.xyz"

os.system(input_in)

# main part - generate oligomer and append the new configuration each time 

tx0=ip["tx"]
ty0=ip["ty"]
tz0=ip["tz"]
twist0=ip["twist"]

dict = {}
nconf=0 
for ntx in range(ip["ntx"]):
  for nty in range(ip["nty"]):
    for ntz in range(ip["ntz"]):
      for ntwist in range(ip["ntwist"]):

        ip["tx"]=tx0+ntx*ip["dtx"]
        ip["ty"]=ty0+nty*ip["dty"]
        ip["tz"]=tz0+ntz*ip["dtz"]
        ip["twist"]=twist0+ntwist*ip["dtwist"]
        nconf+=1 
        print(" LOG: evaluating%7d/%-7d tx=%-5.2f ty=%-5.2f tz=%-5.2f twist=%-5.2f " %(nconf,totconf,ip["tx"],ip["ty"],ip["tz"],ip["twist"]))

        for i in range(1,ip["size"]):
          input_in="python3 ./orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(ip["tz"]) + " -tx " + str(ip["tx"])  + " -ty " + str(ip["ty"])   + " -rz " + str(ip["twist"]) +"  > tmpconfig"+str(i)+".xyz" 

          os.system(input_in)
 
        files = [f for f in os.listdir('.') if os.path.isfile(f) if 'tmpconfig' in f]
        rotated_files = natsorted(files)
        
        with open(output,'a+') as f:
           f.write(str(totatom)+"\n")
           tx=str('{:<0.2f}'.format(ip["tx"]))
           ty=str('{:<0.2f}'.format(ip["ty"]))
           tz=str('{:<0.2f}'.format(ip["tz"]))
           twist=str('{:<0.1f}'.format(ip["twist"]))
           f.write("conf_%-7d %s" %(nconf,tx+"_"+ty+"_"+tz+"_"+twist))
           f.write("\n")
        
           for file in rotated_files:
             with open(file,"r") as fp: 
               fp.readline()
               fp.readline()
               for atoms in range(int(natom)): 
                  line = fp.readline().split()
                  f.write("%s " %(line[0])+" ")
                  for i in range(1,len(line)):
                    f.write("%-7.3f" %(float(line[i]))+" ")
                  f.write("\n")

 
#clean tmp files before changing the parameter 
os.system("rm -f tmpconfig*")
print(" LOG: coordinates are written to " + output)
print(" LOG: energy json init is written to " + output_energy)
print(" LOG: No. of configurations = " + str(nconf)) 
print(" LOG: Successful ")

