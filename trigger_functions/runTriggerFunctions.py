#!/usr/bin/python

import os
import sys

xenonList = [2,3,4,5,7]#[8,9,10]
sourceList = [1,2,3,4,5,6,7,8,9,10,11,12]#[13,14,15,16,17,18,19,20,21,22,23,24]
octetRange = [30, 59]


#for index in range(octetRange[0], octetRange[1], 1):
#    os.system("./findADCthreshold.exe octet %i"%index)

#for index in xenonList:
#    os.system("./findADCthreshold.exe xenon %i"%index)

for index in sourceList:
    os.system("./findADCthreshold.exe source %i"%index)
