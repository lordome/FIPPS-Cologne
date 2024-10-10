#! /opt/anaconda3/bin/python

import os
import sys
import subprocess


if( len(sys.argv) != 3 ): 
    print("Usage: scriptLinearity27Al.py runBegin runEnd ")
    exit()

runBegin = str(int(sys.argv[1]))
runEnd = str(int(sys.argv[2]))



for i in range(64):
    
    
	cmd1 = "screen -dmS run" + str(i) + "  "
	print(cmd1)
	os.system(cmd1)


	cmd2 = "screen -S run" + str(i) + " -X stuff 'scriptLinearity27Al.py " + runBegin + " " + runEnd + "  " + str(i) +  "^M'" 
	print(cmd2)
	os.system(cmd2)

