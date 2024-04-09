# StackGen - A Supramolecular Stack Generator

<pre>
+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
           ____ _____  _    ____ _  ______ _____ _   _
          / ___|_   _|/ \  / ___| |/ / ___| ____| \ | |
          \___ \ | | / _ \| |   | ' / |  _|  _| |  \| |
           ___) || |/ ___ \ |___| . \ |_| | |___| |\  |
          |____/ |_/_/   \_\____|_|\_\____|_____|_| \_|
                 Stack generator for supramolecules
+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
</pre>


StackGen is a user-friendly open-source framework capable of generating energy-optimized one-dimensional supramolecular stack structures with minimal computational cost. StackGen employs PSO (Particle Swarm Optimization) and a semiempirical quantum mechanical approach to identify stable one-dimensional stack structures from a wide array of potential stack configurations resulting from the translation and twisting of neighboring molecules along all three directions. These low-energy configurations can be used as the initial structures for quantum mechanical calculations and molecular simulations.

## Prerequisites

Before using StackGen, ensure you have the following prerequisites installed:

1. Python 3 or above
2. `natsort`, `ase` and `xtb` module for Python
   - You can install it using the following command:
     ```
     pip install natsort
     pip install ase
     pip install xtb
     ```

3. Install the Quantum mechanical simulation software [DFTB+](https://dftbplus.org/download).

## Getting Started

Follow these steps to use StackGen:

1. In the StackGen folder, edit the `STACKGEN_SRC` variable in the `set_initenv.sh` file to include the path of the StackGen folder as follows:
	```
	STACKGEN_SRC='/path/to/the/folder/containing/main.py/file/'
	```

2. Download the DFTB+ software package and either update the `.bashrc` file by adding the path of DFTB+ as an environment variable or add it to the current shell temporarily.

3. Update the path for DFTB parameters in the `set_initenv.sh` file. 

4. Enter the parameter values in the `input.user` file present in the `examples` folder using the following format:
	```
       tx_upper=2.0 
       ty_lower=0.0 
       ty_upper=2.0 
       tz_lower=3.5 
       tz_upper=3.8 
       twist_lower=120 
       twist_upper=180 
       input_struct=crest_conformer_10.xyz
       atoms=7,45,16 
       stack_size=3
       energy_calculator=dftb+
       label=pdi 
	```
Here, crest_conformer_10.xyz is the input file which contains coordinates of conformers. `energy_calculator` parameter has two option: dftb+ and xtb.

5. Source the `set_initenv.sh` file using the following command:
  	 ```
	source set_initenv.sh 
	```


6. Execute the code using the following command:
	```
	stackgen_run > out
	```

  

The optimized stacked configuration for each of the conformer will be saved in a file named `generated_{label}_{stack_size}_conformer_{conformer_number}.xyz`. Individual traj file and swarm_traj file will be generated for each of the conformer. `stacked_conformers.xyz` will contain the coordinates of all of these optimized stacked configurations. 


