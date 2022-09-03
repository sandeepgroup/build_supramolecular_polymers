#!/usr/bin/env python3
# it requires two environment variables
# 1. dftb+ should be on the path, otherwise set the path
#      export ASE_DFTB_COMMAND="/path/to/dftb+ > PREFIX.out"
# 2. SK files: 
#      export DFTB_PREFIX=/path/to/mio-0-1/


import sys
import re
import glob

from ase.calculators.dftb import Dftb
from ase.io import read
from ase import Atoms, Atom

print(" Log: Using the ", glob.glob("config*.xyz"), "as the input file for this script")

filename=glob.glob("config*.xyz")[0]

mol=read(filename)
natoms=len(mol.numbers)

atoms = Atoms(mol.symbols,positions=mol.positions) 
calc = Dftb(atoms=atoms,
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

# write this tetramer configuration along with its energy 
with open("extended_"+filename+".xyz",'w') as fp:
    fp.write(str(natoms)+"\n")
    fp.write(str(potential_en)+"\n")

    for i in range(natoms):
      fp.write(mol.symbols[i]+" ")

      for j in range(len(mol.positions[i])):
        fp.write(str(mol.positions[i][j])+" ")
      fp.write("\n")

