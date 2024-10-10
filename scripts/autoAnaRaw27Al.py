#! /usr/bin/python3

import os
import sys
import subprocess

if( len(sys.argv) != 6 ): 
	print("Automatization of the AnaRaw calibration FOR 27Al runs. \n Please refer to AnaRaw for the parameter description.")
	print("WARNING : The output files get saved inside the folder ./outputAnaRaw. They might get overwritten.")					
	print("Usage: autoAnaRaw.py runBegin runEnd source polyDegree verboseLevel")
	print("Example > autoAnaRaw27Al.py 207421 207583 152Eu 1 2")
	exit()
        
vecRun = range(int(sys.argv[1]), int(sys.argv[2]))

source = " -ener 30.638 -ener 104.6234 -ener 212.017 -ener 271.175 -ener 314.395 -ener 831.450 -ener 846.764 -ener 983.020 -ener 1622.870 -ener 1778.969 -ener 1810.726 -ener 2282.773 -ener 2590.244 -ener 3033.890 -ener 3465.070 -ener 4133.408 -ener 4259.539 -ener 4733.847 -ener 6702.034 -ener 7213.034 -ener 7724.034 "
polyDegree = " -poly " +str(sys.argv[4])
verboseLevel = " -v " + str(sys.argv[5])

calibconf = source + polyDegree + verboseLevel

for i, run in enumerate(range(int(sys.argv[1]), int(sys.argv[2]))):

    cmd2 = "AnaRaw -calib -file ./LSTSort/" + str(int(vecRun[i])) + ".root -calibconf" + calibconf + " > ./outputAnaRaw/" + str(int(vecRun[i])) + ".txt ^M" 

    print(cmd2)
    os.system(cmd2)

