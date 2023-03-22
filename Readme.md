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

Use this code to generate a configuration for the given input parameters such <br />
tx = translation along x-axis <br />
ty = translation along y-axis <br />
tz = translation along z-axis <br />
twist = angle <br />
size = oligomer size 

## use the following command to initialize the environment 

source init.sh 

edit the 'init.sh' file to set the dir_source to correct path 


## Apply patch file to pyswarms directory 

#### for our group: 
To create a patch file, use the following command:
diff -ruN pyswarms/ pyswarms.new/ | sed '/Binary\ files\ /d'   > run.patch

#### for users: 
to apply the patch, use the following command:
patch -p0 < run.patch 

The above command applies the patch to a folder. Therefore, before applying the patch, move to the directory, where
pyswarms is installed. 





