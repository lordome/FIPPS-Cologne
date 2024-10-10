#! /opt/anaconda3/bin/python

import os
import sys
import subprocess

if( len(sys.argv) != 6 ): 
	print("Distribution of EventBuilder over different terminals. \n Please refer to EventBuilder for the parameter description.\n")

	print("WARNING: the processing takes around 6GB of Memory for each terminal - mind the resources available.\n")

	print("Usage: autoEventBuilder.py runBegin runEnd numTerminals confFile overWrite(bool)")
	print("Example > autoEventBuilder.py 207421 207583 5 EventBuilder.conf 0")
	exit()
        
    
def linspace(start, stop, n):
    if n == 1:
        yield int(stop)
        return
    h = (stop - start) / (n - 1)
    for i in range(n):
        yield int(start + h * i)
        
vecRun = list(linspace(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])+1))


confFile = " -conf " + str(sys.argv[4])
overWrite = " "
if(int(sys.argv[5])):
	print("\n \n \t WARNING : OVERWRITE OPTION USED \n \n")
	overWrite = " -ow "

sub = 1
for i in range(int(sys.argv[3])):
    
	if( i == int(sys.argv[3])-1 ): 
		sub = 0
    
	cmd1 = "screen -dmS run" + str(i+100) + "  "
	print(cmd1)
	os.system(cmd1)

	if(int(vecRun[i]) == int(vecRun[i+1])-sub):
		cmd2 = "screen -S run" + str(i+100) + " -X stuff 'EventBuilder " + confFile + " -batch " + overWrite + " -run " + str(int(vecRun[i])) + "^M'" 
		print(cmd2)
		os.system(cmd2)
		continue

	cmd2 = "screen -S run" + str(i+100) + " -X stuff 'EventBuilder " + confFile + " -batch " + overWrite + " -run " + str(int(vecRun[i])) + "-" + str(int(vecRun[i+1])-sub) + "^M'" 
	print(cmd2)
	os.system(cmd2)

