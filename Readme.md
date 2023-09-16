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
StackGen is a user-friendly open-source framework capable of generating energy-optimized one-dimensional supramolecular stack structures with minimal computational cost. StackGen employs PSO and a semiempirical quantum mechanical approach method to identify stable one-dimensional stack structures from a wide array of potential stack configurations resulting from the translation and twisting of neighboring molecules along all three directions. These low-energy configurations can be used as the initial structures for quantum mechanical calculations and molecular simulations. 

Prerequisites:

1. Python 3 or above
2. natsort modules of Python
3. Install 'natsort' module using the below command:<br/>
           
      pip install natsort
           
4. Install Quantum mechanical simulation software dftb+ from https://dftbplus.org/download.

Steps to use the StackGen:

1. In the StackGen folder, edit the <tt>STACKGEN_SRC</tt> variable in the <tt>set_initenv.sh</tt> file to include the path of the StackGen folder as follows: <br/>
      STACKGEN_SRC = ’.../path to the folder containing main.py file/’
2. Download the DFTB+ software package and either update the <tt>.bashrc</tt> file by adding the path of
  DFTB+ as an environment variable or add it to the current shell temporarily.
3. Update the path for DFTB parameters in the <tt>dftb_in.hsd</tt> file. Use the file in <tt>examples</tt> as a 
   starting point for creating the dftb input file.
4. Enter the parameter values in the <tt>input.user</tt> file present in <tt>generate_configuration</tt> folder using      the below format: <br/>
       tx_lower=0.0 <br/>
       tx_upper=2.0 <br/>
       ty_lower=0.0 <br/>
       ty_upper=2.0 <br/>
       tz_lower=3.5 <br/>
       tz_upper=3.8 <br/>
       twist_lower=120 <br/>
       twist_upper=180 <br/>
       input_struct=input.xyz <br/>
       atoms=7,45,16 <br/>
       stack_size=3 <br/>
       label=pdi <br/>
   
 These are basic parameters. Advanced parameters can also be inserted. Details of advanced parameters are given in 
 the user manual. 
 
4. Execute the code using the following command:<br/>
       stackgen_run > out 

5. To clean the temporary files, use the following command: <br/>
       stackgen_clean 

The optimized stacked configuration will be saved in a file named <tt>generated_{label}_{stack_size}.xyz</tt>, which can be visualized using software such as VMD. 
   
