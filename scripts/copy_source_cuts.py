#!/usr/bin/python

import os
import shutil

year = 2012
source_file_base = "../run_lists/"
source_range = []

ref_runs = ( { 1:17238 , 2:17370 , 3:17521 , 4:17892 , 5:18361 , 6:18621 , 7:18749,
               8:19232 , 9:19359 , 10:19511 , 11:19857 , 12:19899 , 13:20519,
               14:20820 , 15:20905 , 16:21091 , 17:21315 , 18:21683 , 19:21918,
               20:22219 , 21:22298 , 22:22441 , 23:22771 , 24:22925} )

cuts_base = os.getenv("CUTS")
if year==2011:
    source_range = [1,2,3,4,5,6,7,8,9,10,11,12]
elif year==2012:
    source_range = [13,14,15,16,17,18,19,20,21,22,23,24]
else:
    exit

for i in source_range:
    srcFileName = source_file_base+"Source_Calibration_Run_Period_%i.dat"%i
    srcFile = open(srcFileName,'r')
    srcRuns = []
    if os.path.isfile(srcFileName):
        for line in srcFile:
            #entries = line.split()
            srcRuns.append(int(line))
            print int(line)
        refRun = ref_runs[i]
        if not os.path.isfile(cuts_base+"cuts_%i.dat"%refRun):
            print "NO CUTS FOR REFERENCE RUN %i IN PERIOD %i"%(refRun,i)
            exit

        for run in srcRuns:
            if run!=refRun:
                shutil.copy(cuts_base+"cuts_%i.dat"%refRun, cuts_base+"cuts_%i.dat"%run)
                #print cuts_base+"cuts_%i.dat"%run

            
        
        
