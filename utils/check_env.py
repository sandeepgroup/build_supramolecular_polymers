#!/bin/bash
# this script will check the presence of all input files needed and then 
# check dftb+ calculation on a test system


import os
import shutil
import sys
import glob
import subprocess
from ase.calculators.dftb import Dftb
from xtb.ase.calculator import XTB
from ase.io import read,write


class ExecutableChecker():

   def __init__(self, exe):
        self.exe = exe
        #self.output = output
        #self.input_file = input_file
        #self.conv_factor = conv_factor
# check whether input env is set up properly or not

   def environment_check(self):
     
     
     if(self.exe.lower()=='dftb+'):
       
       #check if executable exsists in path
       
       if not any(os.access(os.path.join(path, self.exe), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
         return 1
       dftb_prefix = os.environ.get('DFTB_PREFIX')
       if dftb_prefix is not None:
         print(" LOG: dftb+ parameters is set to:", dftb_prefix)
       else:
         return 2

     #if(self.exe.lower()=='xtb'):
       #if not any(os.access(os.path.join(path, self.exe), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
           #return 1
       #dftb_prefix = os.environ.get('DFTB_PREFIX')
       #if dftb_prefix is not None:
         #print(" LOG: dftb+ parameters is set to:", dftb_prefix)
       #else:
         #return 2


#check whether dftb running fine by running water molecule scf calc 
   def run_calculation(self):
     os.environ['OMP_NUM_THREADS'] = '1'
     # Create a new input file with atom coordinates for a water molecule
     with open('water.xyz', 'w') as tmp:
       tmp.write("3\n\n")
       tmp.write("O  -0.23980816   -1.09112708    0.00000000\n")
       tmp.write("H   0.72019184   -1.09112708    0.00000000\n")
       tmp.write("H  -0.56026274   -0.18619125    0.00000000\n")

     mol = read('water.xyz')
     if(self.exe.lower()=='dftb+'):
       calc = Dftb(atoms=mol,
	      label='h2o',
              Hamiltonian_SCCTolerance=1e-8,
              Hamiltonian_MaxAngularMomentum_='',
              Hamiltonian_MaxAngularMomentum_O='p',
              Hamiltonian_MaxAngularMomentum_H='s',
              )
       mol.calc = calc
       total_energy = mol.get_potential_energy()
       # Clean up temporary files
       dftb_files = glob.glob("dftb*")
       for files in dftb_files:
          os.remove(files)
       os.remove("detailed.out")
       os.remove("geo_end.gen")
       os.remove("band.out")
     if(sef.exe.lower()=='xtb'):
       calc = XTB(method="GFN2-xTB")
       mol.calc = calc
       total_energy = mol.get_potential_energy()

     if total_energy is None:
       return 3
     else:
       return 0
