#!/usr/bin/python

import os
import shutil

octIgnore = [7,9,59,60,61,62,63,64,65,66,67,74,81,88,96,103,110,117] # 67,74,81,88,96,103,110,117 these were already done... 


year = 2012
octet_file_base = None
octet_range = []

if year==2011:
    octet_file_base = "%s/2011-2012/"%os.getenv("OCTET_LIST")
    octet_range = [0,59]
elif year==2012:    
    octet_file_base = "%s/2012-2013/"%os.getenv("OCTET_LIST")
    octet_range = [113,121]
else:
    exit

for oct in range(octet_range[0],octet_range[1]+1,1):
    if oct not in octIgnore:
        #os.system("./XuanStyle_DeltaExpProcessor.exe %i"%oct)
        os.system("./DeltaExpProcessor.exe %i"%oct)

