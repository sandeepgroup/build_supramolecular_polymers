#!/bin/bash

import os
import glob 

def clean_up_files():
    files_to_remove = ["band.out","charges.bin","geo_end.gen","dftb_pin.hsd", "tmp.energy", "detailed.out","report.log"]
    for file_name in files_to_remove:
       if os.path.exists(file_name):
          os.remove(file_name)
       else:
          pass
    files_to_remove = glob.glob("tmpconfig*.xyz") + glob.glob("conformer_*.xyz") + glob.glob("dftb*")
    for file_name in files_to_remove:
       if os.path.exists(file_name):
          os.remove(file_name)
    #os.remove("tmpinput1.xyz")
    #os.remove("charges.bin")
    #os.remove("dftb_pin.hsd")
    #os.remove("tmp.energy")
    #os.remove("detailed.out")
    #os.remove("tmpconfig*")

