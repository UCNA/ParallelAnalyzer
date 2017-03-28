#!/usr/bin/python

import os


with open("blind_times_2011-2012.dat",'w') as output:
    for run in range(15000, 20000, 1):
        path = "%s/runInfo_%i.dat"%(os.getenv("RUN_INFO_FILES"),run)
        if os.path.isfile(path):
            infile = open(path,'r')
            lines = infile.readlines()
            output.write("%i\t%f\t%f\n"%(run,float(lines[0].split()[1]),float(lines[1].split()[1])))
            infile.close()


