#!/bin/bash

# this script will check the presence of all input files needed and then 
# check dftb+ calculation on a test system
# go through the trajectory file and for each configuration, it calculates energy
import os 
import re
import glob
from utils.clean import clean_up_files
from ase.build import molecule
from ase.calculators.dftb import Dftb
from xtb.ase.calculator import XTB
from ase.io import read,write

class energy_calculation():
   def __init__(self,input_file,exe):
     self.exe = exe
     #self.output = output
     self.input_file = input_file
     #self.conv_factor = conv_factor

   
    
   def user_input_read(self):
      input_param = {}
      atom_num = 1 
      with open("input.user", "r") as fp:
        for line in fp:
          if not len(line.strip()) == 0 and not line.lstrip().startswith("#"):
            name, var = line.partition("=")[::2]
            var = var.strip()

            if(len(re.findall("[0-9-+]+", var)) > 1 and name == "atoms"):  # More than one element in atoms
              atom_list = re.findall("[0-9-+]+", var)
              for element in atom_list:
                input_param["atom_" + str(atom_num)] = element
                atom_num += 1

            if re.match("^[0-9-+]*$", var):
              var = int(var)
            elif re.match("^[0-9-+.]*$", var):
              var = float(var)
            input_param[name.strip()] = var
      fp.close()
      return input_param
   def dftb_energy_calculation(self):

      input_param = self.user_input_read()
      input_param['nproc_dftb']='2'
      os.environ['OMP_NUM_THREADS'] = input_param['nproc_dftb']
      mol = read(self.input_file)
      
      calc = Dftb(atoms=mol,
	    label='dftb+',
	    Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance=1.e-5,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            Hamiltonian_MaxAngularMomentum_N='p',
            Hamiltonian_MaxAngularMomentum_C='p',
            )


      mol.calc = calc
      total_energy = mol.get_potential_energy()
      total_energy = total_energy * 96.48
      #clean_up_files()
      '''    
      dftb_files = glob.glob("dftb*")
      for files in dftb_files:
         os.remove(files)

      os.remove("detailed.out")
      os.remove("geo_end.gen")
      os.remove("band.out")
      os.remove("charges.bin")
      '''
      return total_energy
   

   def xtb_energy_calculation(self):
      input_param = self.user_input_read()
      input_param['nproc_dftb']='2'
      os.environ['OMP_NUM_THREADS'] = input_param['nproc_dftb']
      mol = read(self.input_file)
      calc = XTB(method="GFN2-xTB")
      mol.calc = calc
      total_energy = mol.get_potential_energy()
      total_energy = total_energy * 96.48

      return total_energy
 
 



