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

STACKGEN is intended to find the optimal structure of different molecules using PSO algorithm for stack formation alongside twisting the molecule and translating it along x, y and z directions. We are using the PySwarms python package for implementing the PSO algorithm for getting the most optimized structure for different configurations.

Prerequisites:

1. Python 3 or above.
2. Install 'pyswarms' and 'natsort' libraries using the below command:<br/>
           
           pip install pyswarms natsort
           
3. Install dftb+ package from the link : https://dftbplus.org/download.

Steps to use the STACKGEN code:

1. In the STACKGEN folder open the file 'set_initenv.sh' and edit the STACKGEN_SRC variable to insert the path of the STACKGEN folder.
2. While staying in the STACKGEN directory in the command prompt run the patch file 'run.patch' using the command<br/>
           
           patch -p0 < run.patch
3. Input the parameters in 'input.user' file present in generate_configuration folder in the below format:<br/>

           filename = input.xyz
           label = pdi
           atom1 = 7
           atom2 = 45
           atom3 = 16
           tx = 1.0
           ty = 0.0
           tz = 3.5
           twist = 30
           size = 3

    These are basic parameters. Advanced parameters can also be inserted. Details of advanced parameters are given in the user manual.
    
4. Finally the code can be run using the following command:<br/>

           python3 main.py > log.txt
   
   The above line will run the main.py file, generate the optimal configuration (in .xyz format) and stores all the results in the log.txt file (the file      will be in the stackgen folder). The generated configuration will be created in a file named as “generated_{label}_{stack_size}.xyz”.
   
Use this code to generate a configuration for the given input parameters such <br />
tx = translation along x-axis <br />
ty = translation along y-axis <br />
tz = translation along z-axis <br />
twist = angle <br />
size = oligomer size 

### use the following command to initialize the environment 

source set_initenv.sh  

edit the 'set_initenv.sh' file to set the dir_source to correct path  


### Apply patch file to pyswarms directory 

#### for our group: 
To create a patch file, use the following command:  
diff -ruN pyswarms/ pyswarms.new/ | sed '/Binary\ files\ /d'   > run.patch

#### for users: 
to apply the patch, use the following command:  
patch -p0 < run.patch 

The above command applies the patch to a folder. Therefore, before applying the patch, move to the directory, where
pyswarms is installed (the command should be executed outside the pyswarms directory). 





