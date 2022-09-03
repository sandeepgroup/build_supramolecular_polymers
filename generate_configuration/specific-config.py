#!/usr/bin/env python3
# generates configuration of given size, translation vector and rotations 
# reads input.user 

import os
import sys
import re 
import glob 
import re
from natsort import natsorted

from ase.calculators.dftb import Dftb
from ase.io import write
from ase.build import molecule
from ase import Atoms, Atom

# read input.user file and store the data as dictionary 
#default values

ip={}

ip["label"]="supra"
ip["tx"]=0.0
ip["ty"]=0.0
ip["tz"]=3.5
ip["twist"]=30.0
ip["size"]=2


with open('input.user','r') as fp:
  for line in fp:
    if not len(line.strip())==0:
      name,var = line.partition('=')[::2]
      var=var.strip()
      if re.match("^[0-9-+]*$", var):
        var=int(var)
      elif re.match("^[0-9-+.]*$", var):   
        var=float(var)

      ip[name.strip()]=var
fp.close()

# re-orient the given configuration so that atom1-atom2 is along x-axis and
# the plane normal is along z-axis 

tmpconfig='tmpconfig'
input_in="python3 ../orient.py "+ ip["filename"] + " -p " + str(ip["atom1"]) +" "+ str(ip["atom2"]) +" " + str(ip["atom3"]) + " -tc > " + tmpconfig +"0.xyz"

os.system(input_in)


for i in range(1,ip["size"]):

  input_in="python3 ../orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(ip["tz"]) + " -tx " + str(ip["tx"])  + " -ty " + str(ip["ty"])   + " -rz " + str(ip["twist"]) +"  > tmpconfig"+str(i)+".xyz" 

  os.system(input_in)
 

#combining the files 

with open(ip["filename"],'r') as fp:
  natom = fp.readline().rstrip()
fp.close()

totatom=int(natom)*ip["size"]

files = [f for f in os.listdir('.') if os.path.isfile(f) if 'tmpconfig' in f]
rotated_files = natsorted(files)
#print(rotated_files)


with open('generated_'+ip["label"]+'_'+str(ip["size"])+'.xyz','w') as f:
    f.write(str(totatom)+"\n")
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
 
print(" LOG: coordinates are written in " + 'generated_'+ip["label"]+'_'+str(ip["size"])+'.xyz')

#clean tmp files 
os.system("rm -f tmpconfig*")
 

