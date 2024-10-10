#! /opt/anaconda3/bin/python

import os
import sys
import subprocess

if( len(sys.argv) != 5 ): 
	print("Automatization of the AnaRaw calibration. \n Please refer to AnaRaw for the parameter description.")
	print("WARNING : The output files get saved inside the folder ./outputAnaRaw. They might get overwritten.")					
	print("Usage: autoAnaRaw.py runBegin runEnd source polyDegree ")
	print("Example > autoAnaRaw.py 207421 207583 152Eu 1 ")
	exit()
        
vecRun = range(int(sys.argv[1]), int(sys.argv[2]))

source = " -source " + str(sys.argv[3])
polyDegree = " -poly " +str(sys.argv[4])


calibconf = source + polyDegree 

for i, run in enumerate(range(int(sys.argv[1]), int(sys.argv[2]))):

    cmd2 = "AnaRaw -calib -file ./LSTSort/" + str(int(vecRun[i])) + ".root -id 1 -calibconf" + calibconf + " ^M" 

    print(cmd2)
    os.system(cmd2)

