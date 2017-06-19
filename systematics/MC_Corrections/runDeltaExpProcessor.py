#!/usr/bin/python

import os
import shutil

octIgnore = [7,59,60,61,62,63,64,65,66,67,91,93,101,107,121]#[7,60,61,62,63,64,65,66] # 67,74,81,88,96,103,110,117 these were already done... 


year = 2011
octet_file_base = None
octet_range = []

if year==2011:
    octet_file_base = "%s/2011-2012/"%os.getenv("OCTET_LIST")
    octet_range = [51,59]
elif year==2012:    
    octet_file_base = "%s/2012-2013/"%os.getenv("OCTET_LIST")
    octet_range = [101,102]# (67,75) (76,84) (85,93) (94,102) (103,111) (112,121) 
else:
    exit

for octet in range(octet_range[0],octet_range[1]+1,1):
    if octet not in octIgnore:
        #os.system("./XuanStyle_DeltaExpProcessor.exe %i"%oct)
        os.system("./DeltaExpProcessor.exe %i"%octet)

