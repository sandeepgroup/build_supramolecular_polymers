#!/usr/bin/env python

# mandatory:  pyswarms 

import os.path
import sys
import pyswarm as ps
from pyswarm import pso
import specific_config
import subprocess
import os 
import re
from natsort import natsorted

from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.single.local_best import LocalBestPSO
from os import environ
import re

# set default values  
# no default values for atom1, atom2, and atom3; they must be entered by the user; otherwise
# the program will stop. 

input_param={}

input_param["tx_lower"]=0.0
input_param["tx_upper"]=3.0
input_param["ty_lower"]=0.0
input_param["ty_upper"]=3.0
input_param["tz_lower"]=3.3
input_param["tz_upper"]=3.7
input_param["twist_lower"]=0.0
input_param["twist_upper"]=50.0
input_param["rot_type"]=0
input_param["stack_size"]=3
input_param["input_struct"]='input.xyz'
input_param["pso_type"]='Globalbest'
input_param["ener_tol"]=1.e-5
input_param["maxiterations"]=1000
input_param["nproc_dftb"]=1
input_param["nparticle_pso"]=20
input_param["label"]='supramolecule'


#Reading input from the user 

with open('user_input.user','r') as fp:
  for line in fp:
    if not len(line.strip())==0 and not line.startswith("#"):
      name,var = line.partition('=')[::2]
      var=var.strip()
     
      if re.match("^[0-9-+]*$", var):
        var=int(var)
      elif re.match("^[0-9-+.]*$", var):   
        var=float(var)
      input_param[name.strip()]=var
fp.close()

# check if atom1, atom2, and atom3 are defined; if not, stop the program
if not "atom1" in input_param or not "atom2" in input_param or not "atom3" in input_param:
    print(" Log: Enter atom numbers for the atom1, atom2, and atom3")
    print(" Log: Exit ... ") 
    sys.exit()

# check if atom1, atom2, and atom3 are defined; if not, stop the program
if not os.path.isfile(input_param["input_struct"]):
    print(" Log: xyz file "+input_param["input_struct"]+ " not found")
    print(" Log: Exit ... ") 
    sys.exit()

# print the input parameters along with default settings, if any 
for key in input_param:
  print( " LOG: "+key+" = "+str(input_param[key]))


global atom1,atom2,atom3,input_struct,pso_type,nproc,rotation_type


#renaming the variables for simplicity 

ub = [input_param['tx_upper'],input_param['ty_upper'],input_param['tz_upper'],input_param['twist_upper']]
lb = [input_param['tx_lower'],input_param['ty_lower'],input_param['tz_lower'],input_param['twist_lower']]
rot_type = input_param['rot_type']
atom1 = input_param['atom1']
atom2 = input_param['atom2']
atom3 = input_param['atom3']
input_struct = input_param['input_struct']
rotation_type=input_param["rot_type"]
size = input_param['stack_size']
pso_type = input_param['pso_type']
ftol = input_param['ener_tol']
iterations = input_param['maxiterations'] 
nproc=input_param['nproc_dftb']
label=input_param["label"]

n_particles=input_param["nparticle_pso"]

bounds = (lb,ub)

itera=0

def env_check():
    completedProc = subprocess.run('./check_env.sh')
    if(completedProc.returncode==1):
        print(" ERROR: dftb+ input file is not found")
        exit()
    elif(completedProc.returncode==2):
        print(" ERROR: dftb+ executable not found on the path")
        exit()
    elif(completedProc.returncode==3):
        print(" ERROR: dftb+ test run not successful. Check dftb+ input file")
        print(" ERROR: stop ... ")
        exit()
    else:
        return 1
        
def change(itera):
    itera+=1
    return itera
    
def energy_min(params):
    energy_values_list=[]
    itera=change(itera)
    if itera==1:
        print(" LOG: iteration,pso_particle,tx,ty,tz,twist,Energy")
    x=0
    for parm in range(len(params)):
        x+=1
        tx = params[parm][0]
        ty = params[parm][1]
        tz = params[parm][2]
        twist = params[parm][3]         
        configuration_generate(tx,ty,tz,twist)
        try:
            output = subprocess.check_output('./run_dftb_updated.sh generated_'+label+'_'+str(size)+'.xyz' ,shell=True)
            #output = subprocess.check_output('./run_dftb_updated.sh',shell=True)
        except subprocess.CalledProcessError as grepexc:     
            print(" ERROR: error code", grepexc.returncode, grepexc.output)
        if(re.findall(r"[-+]?(?:\d*\.\d+|[eE][+-]\d+)", str(output.strip()))):
            energy_val = float(re.findall(r"[-+]?(?:\d*\.\d+|[eE][+-]\d+)", str(output.strip()))[0])
            energy_values_list.append(energy_val)   
        print(" LOG: %4d %3d %0.2f %0.2f %0.2f %0.2f %0.5f" %(itera,x,tx,ty,tz,twist,energy_val))
    return energy_values_list

def opt_struct(params):
    print("Generating the optimized structure")
    tx = params[0]
    ty = params[1]
    tz = params[2]
    twist = params[3]
    configuration_generate(tx,ty,tz,twist)


# re-orient the given configuration so that atom1-atom2 is along x-axis and
# the plane normal is along z-axis 

def configuration_generate(tx,ty,tz,twist):
  tmpconfig='tmpconfig'
  input_in="python3 orient.py "+ input_struct + " -p " + str(atom1) +" "+ str(atom2) +" " + str(atom3) + " -tc	> " + tmpconfig +"0.xyz"
  os.system(input_in)
  #zigzag configuration
  if(rotation_type==1):
    rot = twist
    orient = -1
    for i in range(1,size):
        rot = rot * orient
        input_in="python3 orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(tz)  + " -rz " + str(rot) +"  > tmpconfig"+str(i)+".xyz"
        os.system(input_in)

  #changing point of rotation 
  elif(rotation_type==2):
      for i in range(1,size):
          input_in="python3 orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(tz) + " -tx " + str(tx)  + " -ty " + str(ty)   + " -rz " + str(twist) +"  > tmpconfig"+str(i)+".xyz"
          os.system(input_in)

  #normal rotation rotation_type=0 
  else:
      for i in range(1,size):
          input_in="python3 orient.py "+ tmpconfig+str(i-1)+".xyz" + " -tz " + str(tz) + " -tx " + str(tx)  + 	" -ty " + str(ty)   + " -rz " + str(twist) +"  > tmpconfig"+str(i)+".xyz"
          os.system(input_in)

  #combining the files 
  with open(input_struct,'r') as fp:
  	natom = fp.readline().rstrip()
  fp.close()
  
  totatom=int(natom)*size
  
  files = [f for f in os.listdir('.') if os.path.isfile(f) if 'tmpconfig' in f]
  rotated_files = natsorted(files)
  #print(rotated_files)
  
  with open('generated_'+label+'_'+str(size)+'.xyz','w') as f:
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
  
  #print(" LOG: coordinates are written in " + 'generated_'+ip["label"]+'_'+str(ip["size"])+'.xyz')
  #clean tmp files 
  os.system("rm -f tmpconfig*")
        

if(env_check()):
    print(" LOG: All inputs are set")
    
    if(pso_type=='Localbest'):
        print(" LOG: Running LocalbestPSO algorithm")
        options = {'c1': 0.5, 'c2': 0.3, 'w':0.9, 'k': 6, 'p': 2}  
        optimizer = LocalBestPSO(n_particles=n_particles, dimensions=4, options=options,bounds=bounds,ftol=ftol)
        
    elif(pso_type=='Globalbest'):
        print(" LOG: Running GlobalbestPSO algorithm")
        options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}    
        optimizer = GlobalBestPSO(n_particles=n_particles, dimensions=4, options=options,bounds=bounds,ftol=ftol)
    
    cost, pos = optimizer.optimize(energy_min,iterations)
    print(" LOG: the best cost and The best position returned by PSO")
    print(" LOG: ", cost,pos)
    #print(" LOG: %0.2f %0.2f %0.2f %0.2f %0.5e" %(pos,cost))
    opt_struct(pos)



