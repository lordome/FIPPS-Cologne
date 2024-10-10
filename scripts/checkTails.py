#!/opt/anaconda3/bin/python3


import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if( len(sys.argv) != 2 ): 
    print("Usage: scriptLinearity.py runBegin ")
    exit()

runBegin = int(sys.argv[1])


append = False
runningCount = 0
lines = []

countClover = 0
countDetector = 0
countDetectorTotale = 0

df = np.zeros((6,10))

list_files = [str(i) + ".txt" for i in range(runBegin, runBegin+1)]

list_files = ["outputAnaRaw/"+i for i in list_files]

df_1778 = []
df_7724 = []

df_121 =[]
df_778 = []        
df_1408   =[]      


for countFile, iterfiles in enumerate(list_files):
        
    with open(iterfiles, 'r') as fp:


        for i, line in enumerate(fp):

            if(line.strip().endswith('calibrated energy= 0')):   
                    runningCount = 200
                    
            if(line.strip().endswith("rChi2%")):
                countDetector += 1 
                countDetectorTotale += 1
                append = True
                continue
                #print("appendTrue")
                
            if(append == True ):

                if(1): 
                    lines.append(line.strip())
                    #print(line.strip())
                    
                runningCount +=1 
                
                if(line.strip().endswith('calibrated energy= 0')): 
                    lines = lines[:-1]
                    runningCount = 100

                        

fig, ax = plt.subplots( figsize=(8, 4 ) )

        
ax.set_yticks(np.arange(0.,7, 0.2))
    
df_plot = pd.DataFrame(df_plot)
ax.bar(df_plot.iloc[:32,0], df_plot.iloc[:32,1], alpha = 0.5, label = 'FIPPS')
ax.bar(df_plot.iloc[32:,0], df_plot.iloc[32:,1], alpha = 0.5, label = 'IFIN')
ax.axhline(y=2, xmin = 0, xmax = 32.5/66, color='b', linestyle='-', label = '2.0 FIPPS reference')
ax.axhline(y=2.2, xmin = 32.5/66, xmax = 1, color='tab:orange', linestyle='-', label = '2.2 IFIN reference' )

ax.legend(prop={'size': 5})
ax.set_ylabel("FWHM")
ax.set_xlabel("Clover ID", fontsize = 10)
ax.set_xlim( -1, 65 )

ax.set_title("RunNumber: ")

ax.grid( 'off', axis='x', which='minor', alpha = 0.5 )
ax.grid( 'off', axis='y', alpha = 0.2)

figName = "outputPictures/outputFWHM"  + ".png"
plt.savefig(figName)

plt.show()

