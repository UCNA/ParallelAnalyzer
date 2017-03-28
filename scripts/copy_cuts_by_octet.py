#!/usr/bin/python

import os
import shutil

year = 2012
octet_file_base = None
octet_range = []
cuts_base = os.getenv("CUTS")
if year==2011:
    octet_file_base = "%s/2011-2012/"%os.getenv("OCTET_LIST")
    octet_range = [0,59]
elif year==2012:
    octet_file_base = "%s/2012-2013/"%os.getenv("OCTET_LIST")
    octet_range = [60,121]
else:
    exit

for i in range(octet_range[0],octet_range[1]+1,1):
    octFileName = octet_file_base+"octet_list_%i.dat"%i
    octFile = open(octFileName,'r')
    octRuns = []
    if os.path.isfile(octFileName):
        for line in octFile:
            entries = line.split()
            octRuns.append(int(entries[1]))
            print int(entries[1])
        refRun = None
        for run in octRuns:
            if os.path.isfile(cuts_base+"cuts_%i.dat"%run):
                refRun = run
                continue

        if refRun:
            for run in octRuns:
                if run!=refRun:
                    shutil.copy(cuts_base+"cuts_%i.dat"%refRun, cuts_base+"cuts_%i.dat"%run)


            
        
        
