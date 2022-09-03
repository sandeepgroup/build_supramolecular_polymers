#!/usr/bin/env python3

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

# if files with rotated*xyz are present, the code will exit. 

if len(glob.glob("rotated*xyz")) != 0:
  print(" ERROR: rotated*xyz files are present in current directory")
  print(" ERROR: Delete rotated*xyz files and run the script again")
  print(" ERROR: Exit")
  sys.exit() 

# read input.user file and store the data as dictionary 

ip={}
with open('input.user','r') as fp:
  for line in fp:
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
input_in="python3 ./orient.py "+ ip["filename"] + " -p " + str(ip["atom1"]) +" "+ str(ip["atom2"]) +" " + str(ip["atom3"]) + " -tc > " + tmpconfig +"0.xyz"

print(os.system(input_in))


for i in range(1,ip["size"]):

    input_in="python3 ./orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(ip["tz"]) + " -tx " + str(ip["tx"])  + " -ty " + str(ip["ty"])   + " -rz " + str(ip["twist"]) +"  > tmpconfig"+str(i)+".xyz"

    print(os.system(input_in))
 

#combining the files 

with open(ip["filename"],'r') as fp:
  natom = fp.readline().rstrip()
fp.close()

totatom=int(natom)*ip["size"]

files = [f for f in os.listdir('.') if os.path.isfile(f) if 'tmpconfig' in f]
rotated_files = natsorted(files)
#print(rotated_files)

coor=[]
with open('output_config.xyz','w') as f:
    f.write(str(totatom)+"\n")
    f.write("\n")
    for file in rotated_files:
        fp = open(file,'r')
        for line in fp:
            if(len(line) ==3 or len(line)==1):
                continue
            else:
                line = line.strip()
                f.write(line)
                f.write("\n")
                coor.append(line)
 
#clean tmp files 
print(os.system("rm -f tmpconfig*"))
sys.exit() 
 
#extracting the coordinates only
coordinates=[]
for line in coor:
    pos = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
    pos = [float(i) for i in pos]
    coordinates.append(pos)
    
#DFTB+ calculation in ASE environment

atoms = Atoms('C72H36N4O8',positions=coordinates) #monomer is C36H18N2O4
calc = Dftb(atoms=atoms,
            label='C72H36N4O8',
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance=1e-8,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_C='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            Hamiltonian_MaxAngularMomentum_N='p',
            Hamiltonian_MaxAngularMomentum_O='p',
            )


atoms.calc = calc
potential_en = atoms.get_potential_energy()
print(potential_en)
#dyn = BFGS(atoms, trajectory='test.traj')
#dyn.run(fmax=0.01)
#write('final1.xyz', atoms)

# we have to write this tetramer configuration along with its energy 
with open("output_config_enery.xyz",'w') as fp:
    fp.write('240')
    fp.write("\n")
    fp.write(str(potential_en))
    fp.write("\n")
    coor_fp = open("output_config.xyz",'r')
    for line in coor_fp:
        if(re.findall(r"[-+]?(?:\d*\.\d+)", line)):
            fp.write(line)
    

