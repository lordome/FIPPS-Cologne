#! /usr/bin/python3

import os
import sys
import subprocess

if( len(sys.argv) != 6 ): 
	print("Distribution of LSTCloverBuilder over different terminals. \n Please refer to LSTCloverBuilder for the parameter description.\n")

	print("WARNING: the processing takes around 6GB of Memory for each terminal - mind the resources available.\n")

	print("Usage: autoLSTCloverBuilder.py runBegin runEnd numTerminals confFile overWrite(bool)")
	print("Example > autoLSTCloverBuilder.py 207421 207583 5 LSTClovBuilder.conf 0")
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
    
	cmd1 = "screen -dmS run4" + str(i) + "  "
	print(cmd1)
	os.system(cmd1)

	if(int(vecRun[i]) == int(vecRun[i+1])-sub):

		cmd2 = "screen -S run4" + str(i) + " -X stuff ' module load exp_3-17-58 ^M'" 
		print(cmd2)
		os.system(cmd2)
		cmd2 = "screen -S run4" + str(i) + " -X stuff 'LSTClovBuilder " + confFile + " -batch " + overWrite + " -run " + str(int(vecRun[i])) + "^M'" 
		print(cmd2)
		os.system(cmd2)
		continue


	cmd2 = "screen -S run4" + str(i) + " -X stuff ' module load exp_3-17-58 ^M'" 
	print(cmd2)
	os.system(cmd2)

	cmd2 = "screen -S run4" + str(i) + " -X stuff 'LSTClovBuilder " + confFile + " -batch " + overWrite + " -run " + str(int(sys.argv[1])) + "-" + str(int(sys.argv[2])) + "^M'" 
	print(cmd2)
	os.system(cmd2)

