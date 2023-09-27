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
2. `natsort` module for Python
   - You can install it using the following command:
     ```
     pip install natsort
     ```

3. Install the Quantum mechanical simulation software [DFTB+](https://dftbplus.org/download).

## Getting Started

Follow these steps to use StackGen:

1. In the StackGen folder, edit the `STACKGEN_SRC` variable in the `set_initenv.sh` file to include the path of the StackGen folder as follows:
	```
	STACKGEN_SRC='/path/to/the/folder/containing/main.py/file/'
	```

2. Download the DFTB+ software package and either update the `.bashrc` file by adding the path of DFTB+ as an environment variable or add it to the current shell temporarily.

3. Update the path for DFTB parameters in the `dftb_in.hsd` file. You can use the file in the `examples` folder as a starting point for creating the DFTB input file.

4. Enter the parameter values in the `input.user` file present in the `examples` folder using the following format:
	```
       tx_upper=2.0 
       ty_lower=0.0 
       ty_upper=2.0 
       tz_lower=3.5 
       tz_upper=3.8 
       twist_lower=120 
       twist_upper=180 
       input_struct=input.xyz 
       atoms=7,45,16 
       stack_size=3 
       label=pdi 
	```

  These are basic parameters. Advanced parameters can also be inserted. Details of advanced parameters are given in the user manual.

5. Execute the code using the following command:
	```
	stackgen_run > out
	```

6. To clean the temporary files, use the following command:
	```
       stackgen_clean 
	```


The optimized stacked configuration will be saved in a file named `generated_{label}_{stack_size}.xyz`, which can be visualized using software such as VMD (Visual Molecular Dynamics).


