#!/usr/bin/env python
# coding: utf-8

import os.path
import sys
import subprocess
import os 
import re
import pyswarms as ps
import numpy as np
import time 
from natsort import natsorted
from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.single.local_best import LocalBestPSO
from pyswarms.backend.operators import compute_pbest, compute_objective_function

from collections import deque


class bcolors:
    HEADER = '\033[92m'
    BLUE = '\033[94m'
    WARNING = '\033[91m'
    BOLD = '\033[1m'
    END = '\033[0m'

# print header in the log file 

print(bcolors.HEADER + bcolors.BOLD + r"""
+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
           ____ _____  _    ____ _  ______ _____ _   _
          / ___|_   _|/ \  / ___| |/ / ___| ____| \ | |
          \___ \ | | / _ \| |   | ' / |  _|  _| |  \| |
           ___) || |/ ___ \ |___| . \ |_| | |___| |\  |
          |____/ |_/_/   \_\____|_|\_\____|_____|_| \_|
                 Stack generator for supramolecules
+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
"""
+bcolors.END 
)

# set default values  
# no default values for atom1, atom2, and atom3; they must be entered by the user; otherwise
# the program will stop. 

input_param={}
start_opts={}
end_opts={}
options={}

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
input_param["ener_tol"]=10
input_param["ener_tol_iter"]=5
input_param["maxiterations"]=100
input_param["nproc_dftb"]=1
input_param["nparticle_pso"]=10
input_param['c1']=0.5
input_param['c2']=2.5
input_param['w']=9
input_param['dimensions']=1
input_param['c1_start']=0.5
input_param['c1_end']=2.5
input_param['c2_start']=0.5
input_param['c2_end']=2.5
input_param['w_start']=0.9
input_param['w_end']=0.4
input_param['k']=6
input_param['p']=2
input_param["label"]='supramolecule'


#Reading input from the user 

with open('user_input.user','r') as fp:
  for line in fp:
    if not len(line.strip())==0 and not line.lstrip().startswith("#"):
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


# check if input xyz file is present 
if not os.path.isfile(input_param["input_struct"]):
    print(" Log: xyz file "+input_param["input_struct"]+ " not found")
    print(" Log: Exit ... ") 
    sys.exit()

# print the input parameters along with default settings, if any 
print( " LOG: print all key and values used by the code ")
print( " LOG: note that some of the keys displayed here may not be used by the code")
for key in input_param:
  print( " LOG: "+key+" = "+str(input_param[key]))


global atom1,atom2,atom3,input_struct,pso_type,nproc,rotation_type,neighbour,distance,n_particles,bounds,dimensions,param_len

#renaming the variables for simplicity 

tx_upper = input_param['tx_upper']
tx_lower = input_param['tx_lower']
ty_upper = input_param['ty_upper']
ty_lower = input_param['ty_lower']
tz_upper = input_param['tz_upper']
tz_lower = input_param['tz_lower']
twist_upper = input_param['twist_upper']
twist_lower = input_param['twist_lower']
rot_type = input_param['rot_type']
atom1 = input_param['atom1']
atom2 = input_param['atom2']
atom3 = input_param['atom3']
input_struct = input_param['input_struct']
rotation_type=input_param["rot_type"]
size = input_param['stack_size']
pso_type = input_param['pso_type']
ftol = input_param['ener_tol']
ftol_iter = input_param['ener_tol_iter']
iterations = input_param['maxiterations'] 
nproc=input_param['nproc_dftb']
dimensions=input_param['dimension']
options['c1']=input_param['c1']
options['c2']=input_param['c2']
options['w']=input_param['w']

start_opts['c1'] = input_param['c1_start']
start_opts['c2'] = input_param['c2_start']
start_opts['w'] = input_param['w_start']
end_opts['c1'] = input_param['c1_end']
end_opts['c2'] = input_param['c2_end']
end_opts['w'] = input_param['w_end']

neighbour=input_param['k']
distance=input_param['p']
oh_strategy = eval(input_param['oh_strategy'])
label=input_param["label"]
n_particles=input_param["nparticle_pso"]


if(dimensions==1):
    #tx=0,ty=0,tz=3.5
    tx_upper=tx_lower=0
    ty_upper=ty_lower=0
    tz_upper=tz_lower=3.5
    ub=[twist_upper]
    lb=[twist_lower]
    param_len=1
elif(dimensions==2):
        #tx and ty=0
    tx_upper=tx_lower=0
    ty_upper=ty_lower=0
    ub=[tz_upper,twist_upper]
    lb=[tz_lower,twist_lower]
    param_len=2
elif(dimensions==3):
        #ty=0
    ty_upper=ty_lower=0
    ub=[tx_upper,tz_upper,twist_upper]
    lb=[tx_lower,tz_lower,twist_lower]
    param_len=3
elif(dimensions==4):
    ub = [input_param['tx_upper'],input_param['ty_upper'],input_param['tz_upper'],input_param['twist_upper']]
    lb = [input_param['tx_lower'],input_param['ty_lower'],input_param['tz_lower'],input_param['twist_lower']]
    param_len=4
    
else:
    print(bcolors.WARNING + "ERROR: Wrong value is given.Check user.input file" + bcolors.END)
    
bounds =(lb,ub)    


def env_check():
    completedProc = subprocess.run('./check_env.sh')
    if(completedProc.returncode==1):
        print(bcolors.WARNING + " ERROR: dftb+ input file is not found" + bcolors.END)
        exit()
    elif(completedProc.returncode==2):
        print(bcolors.WARNING + " ERROR: dftb+ executable not found on the path"+ bcolors.END )
        exit()
    elif(completedProc.returncode==3):
        print(bcolors.WARNING +" ERROR: dftb+ test run not successful. Check dftb+ input file"+ bcolors.END)
        print(bcolors.WARNING +" ERROR: check dftb+ parameters directory path"+ bcolors.END)
        print(bcolors.WARNING +" ERROR: note that the directory path should terminate with /"+ bcolors.END)
        print(bcolors.WARNING +" ERROR: stop ... "+ bcolors.END)
        exit()
    else:
        return 1
 

# SKR: ftol_iter stopping criteria implemented 
def optimize(objective_func, maxiters, oh_strategy,start_opts, end_opts):
    global pso_type,neighbour,distance,n_particles,bounds,dimensions
   
    if(pso_type=='Globalbest'):
        opt = ps.single.GlobalBestPSO(n_particles, dimensions=dimensions, options=start_opts,
                                      bounds=bounds, oh_strategy=oh_strategy, ftol=ftol, ftol_iter=ftol_iter)
    elif(pso_type=='Localbest'):
        opt = ps.single.LocalBestPSO(n_particles, dimensions=dimensions, options=start_opts, bounds=bounds, 
                                     oh_strategy=oh_strategy, ftol=ftol, ftol_iter=ftol_iter)
   
    swarm = opt.swarm
    opt.bh.memory = swarm.position
    opt.vh.memory = swarm.position
    swarm.pbest_cost = np.full(opt.swarm_size[0], np.inf)

    ftol_history = deque(maxlen=ftol_iter)
    for i in range(maxiters):
   
        # Compute cost for current position and personal best
        swarm.current_cost =  compute_objective_function(swarm, objective_func)
        swarm.pbest_pos, swarm.pbest_cost = compute_pbest(swarm)

  
        best_cost_yet_found = swarm.best_cost
        if(pso_type=='Globalbest'):
            swarm.best_pos, swarm.best_cost = opt.top.compute_gbest(swarm)
        if(pso_type=='Localbest'):
            swarm.best_pos, swarm.best_cost = opt.top.compute_gbest(swarm,p=distance,k=neighbour)

        print(" LOG: " + pso_type + "\t" + 'best_cost ='+ str(swarm.best_cost)+"\n")

        delta = (
                np.abs(swarm.best_cost - best_cost_yet_found)
                < ftol
            )
        if i < ftol_iter:
            ftol_history.append(delta)
        else:
            ftol_history.append(delta)
            if all(ftol_history):
                break 

        swarm.options = opt.oh( opt.options, iternow=i, itermax=maxiters, end_opts=end_opts )
        print("Iteration:", i," Options: ", swarm.options)    

        swarm.velocity = opt.top.compute_velocity(
            swarm, opt.velocity_clamp, opt.vh, opt.bounds
        )
        swarm.position = opt.top.compute_position(
            swarm, opt.bounds, opt.bh
        )
    final_best_cost = swarm.best_cost.copy()
    final_best_pos = swarm.pbest_pos[
        swarm.pbest_cost.argmin()
    ].copy()
    return final_best_cost, final_best_pos

counter = 0
def count():
    global counter
    counter += 1
    
def energy_cal(tx,ty,tz,twist):
        configuration_generate(tx,ty,tz,twist)
        try:
            output = subprocess.check_output('./run_dftb_updated.sh generated_'+label+'_'+str(size)+'.xyz' ,shell=True)
            #output = subprocess.check_output('./run_dftb_updated.sh',shell=True)
        except subprocess.CalledProcessError as grepexc:     
            print(bcolors.WARNING +" ERROR: error code", grepexc.returncode, grepexc.output+ bcolors.END)
        if(re.findall(r"[-+]?(?:\d*\.\d+|[eE][+-]\d+)", str(output.strip()))):
            energy_val = float(re.findall(r"[-+]?(?:\d*\.\d+|[eE][+-]\d+)", str(output.strip()))[0])
        return energy_val
    
    
def energy_min(params):
    global param_len

    energy_value_list=[]
    count()
    if counter==1:
        print(" LOG: iteration,pso_particle,tx,ty,tz,twist,Energy")
    x=0
    if(param_len==1): #only twist parameter
        for parm in range(len(params)):
            x+=1
            tx=ty=0
            tz=3.5
            twist = params[parm][0]  
            energy_val = energy_cal(tx,ty,tz,twist)
            energy_value_list.append(energy_val)
            print(" LOG: %4d %3d %0.2f %0.2f %0.2f %0.2f %0.6f" %(counter,x,tx,ty,tz,twist,energy_val))
        return energy_value_list
    elif(param_len==2):#tz,twist
        for parm in range(len(params)):
            x+=1
            tx=ty=0
            tz = params[parm][0]
            twist = params[parm][1]
            energy_val = energy_cal(tx,ty,tz,twist)
            energy_value_list.append(energy_val)
            print(" LOG: %4d %3d %0.2f %0.2f %0.2f %0.2f %0.6f" %(counter,x,tx,ty,tz,twist,energy_val))
        return energy_value_list
    elif(param_len==3):#tx,tz,twist 
        for parm in range(len(params)):
            x+=1
            ty=0
            tx = params[parm][0]
            tz = params[parm][1]
            twist = params[parm][2] 
            energy_val = energy_cal(tx,ty,tz,twist)
            energy_value_list.append(energy_val)
            print(" LOG: %4d %3d %0.2f %0.2f %0.2f %0.2f %0.6f" %(counter,x,tx,ty,tz,twist,energy_val))
        return energy_value_list
    elif(param_len==4):#tx,tz,twist
        for parm in range(len(params)):
            x+=1
            tx = params[parm][0]
            ty = params[parm][1]
            tz = params[parm][2]
            twist = params[parm][3] 
            energy_val = energy_cal(tx,ty,tz,twist)
            energy_value_list.append(energy_val)
            print(" LOG: %4d %3d %0.2f %0.2f %0.2f %0.2f %0.6f" %(counter,x,tx,ty,tz,twist,energy_val))
        return energy_value_list

def opt_struct(params):
    print(" Log: Generating the final structure")
    global param_len
    if(param_len==1): #only twist parameter
            tx=ty=0
            tz=3.5
            twist = params[0]  
            configuration_generate(tx,ty,tz,twist)
    elif(param_len==2):
            tx=ty=0
            tz = params[0]
            twist = params[1]
            configuration_generate(tx,ty,tz,twist)         
    elif(param_len==3):#tx,tz,twist 
            ty=0
            tx = params[0]
            tz = params[1]
            twist = params[2] 
            configuration_generate(tx,ty,tz,twist)
    elif(param_len==4):#tx,tz,twist
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
  os.system("rm -f tmpconfig*")
        
if(env_check()):
    print(" LOG: All inputs are set \n")
    
    if(pso_type=='Localbest'):
        print(" LOG: Running LocalbestPSO algorithm")
        options['k']=neighbour
        options['p'] = distance
        #two conditons : oh_strategy=False and oh_Staragy=True
        if(oh_strategy == False): 
            pso_type='Localbest'
            print("LocalbestPSO and oh_strategy=0")
            
            optimizer = LocalBestPSO(n_particles=n_particles, dimensions=dimensions,
                                     options=options,bounds=bounds,ftol=ftol,ftol_iter=ftol_iter)
            cost, pos = optimizer.optimize(energy_min,iterations)
            
        elif(oh_strategy == True):
            print("calculation with oh_strategy and localbestpso")
            pso_type='Localbest'
            start_opts['k']=end_opts['k']=neighbour
            start_opts['p']=end_opts['p']=distance
            oh_strategy={ "w":'exp_decay', "c1":'nonlin_mod',"c2":'lin_variation'}
            cost, pos=optimize(energy_min, iterations, oh_strategy, start_opts, end_opts)
        
        
        
    elif(pso_type=='Globalbest'):
        print(" LOG: Running GlobalbestPSO algorithm")
        
        if(oh_strategy==False): #default
            print("Globalbest with oh_strategy=0")
            pso_type='Globalbest'
            optimizer = GlobalBestPSO(n_particles=n_particles, dimensions=dimensions,options=options,
                                      bounds=bounds,ftol=ftol,ftol_iter=ftol_iter)
            cost, pos = optimizer.optimize(energy_min,iterations)
            
        elif(oh_strategy==True):
            pso_type='Globalbest'
            print("calculation with oh_strategy and globalbest")
            oh_strategy={ "w":'exp_decay', "c1":'nonlin_mod',"c2":'lin_variation'}
            cost, pos=optimize(energy_min, iterations, oh_strategy, start_opts, end_opts)
                
        
     

    print(" LOG: the best cost and The best position returned by PSO")
    print(" LOG: ", cost,pos)
    opt_struct(pos)

    

