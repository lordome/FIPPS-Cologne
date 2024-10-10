#! /usr/bin/python3

import os
import sys
import subprocess

if( len(sys.argv) != 6 ): 
	print("Automatization of the AnaRaw calibration. \n Please refer to AnaRaw for the parameter description.")
	print("WARNING : The output files get saved inside the folder ./outputAnaRaw. They might get overwritten.")					
	print("Usage: autoAnaRaw.py runBegin runEnd source polyDegree verboseLevel")
	print("Example > autoAnaRaw.py 207421 207583 152Eu 1 2")
	exit()
        
vecRun = range(int(sys.argv[1]), int(sys.argv[2]))

source = " -source " + str(sys.argv[3])
polyDegree = " -poly " +str(sys.argv[4])
verboseLevel = " -v " + str(sys.argv[5])

calibconf = source + polyDegree + verboseLevel

for i, run in enumerate(range(int(sys.argv[1]), int(sys.argv[2]))):

    cmd2 = "AnaRaw -calib -file ./LSTSort/" + str(int(vecRun[i])) + ".root -calibconf" + calibconf + " > ./outputAnaRaw/" + str(int(vecRun[i])) + ".txt ^M" 

    print(cmd2)
    os.system(cmd2)

