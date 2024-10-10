#!/opt/anaconda3/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if( len(sys.argv) != 2 ): 
    print("Usage: python3 displayFWHM.py runNum")
    exit()
    
    
append = False
runningCount = 0
lines = []

countClover = 0
countDetector = 0

df_plot = []

countDetectorTotale = 0

filenumber = sys.argv[1]
filename = r"outputAnaRaw/" + str(int(filenumber)) + ".txt"

with open(filename, 'r') as fp:

    for i, line in enumerate(fp):
        
        if(line.strip().endswith('calibrated energy= 0')):   
                runningCount = 100
       
        if(append == True): 
            lines.append(line.strip())
            runningCount +=1 
            
            if(line.strip().endswith('calibrated energy= 0')): 
                lines = lines[:-1]
                runningCount = 100

        if(line.strip().endswith("rChi2%")):
            countDetector += 1 
            countDetectorTotale += 1
        
        if(line.strip().endswith("Residue (keV)")):
            append = True
        
        if(runningCount > 9): 
            runningCount = 0
            append = False
            
            a = np.array((([np.fromstring(i[3:], dtype=float, sep=' ') for i in lines[0:]], )), dtype=object)
            if ( a.shape != (1,10,6) ):
                a = np.zeros((10,6))
                a[:,4] = df[4]
            df = pd.DataFrame( a.reshape(10,6) )
            df_plot.append([countDetectorTotale-1, df.iloc[9,3]])            
            
            lines = []
            
            if(countDetector > 3):
                countDetector = 0 
                countClover +=1

fig, ax = plt.subplots( figsize=(8, 4 ) )
ax.set_yticks( np.linspace(0, 200, 11 ) )

vheight = []
xticks= []
xlbls = []
for i in range(64):
    if((i-1)%4 == 0):
        xlbls.append(str(int(i/4 + 1)))
        xticks.append(i+0.5)
        vheight.append(0)

xticks_minor = np.linspace(0,64,17, dtype = int)
ax.set_xticks( xticks )
ax.set_xticklabels( xlbls, fontsize = 12 )
ax.set_xticks( np.array(xticks_minor) - 0.52, minor=True )
ax.tick_params( axis='x', which='minor', length=10) 
ax.tick_params( axis='x', which='major', bottom='off', top='off', labelsize=8, length = 0)

ax.set_yticks(np.arange(0.,2.8, 0.2))
    
df_plot = pd.DataFrame(df_plot)
ax.bar(df_plot.iloc[:32,0], df_plot.iloc[:32,1], alpha = 0.5, label = 'FIPPS')
ax.bar(df_plot.iloc[32:,0], df_plot.iloc[32:,1], alpha = 0.5, label = 'IFIN')
ax.axhline(y=2, xmin = 0, xmax = 32.5/66, color='b', linestyle='-', label = '2.0 FIPPS reference')
ax.axhline(y=2.2, xmin = 32.5/66, xmax = 1, color='tab:orange', linestyle='-', label = '2.2 IFIN reference' )

ax.legend(prop={'size': 5})
ax.set_ylabel("FWHM")
ax.set_xlabel("Clover ID", fontsize = 10)
ax.set_xlim( -1, 65 )

ax.set_title("RunNumber: " + str(int(filenumber)))

ax.grid( 'off', axis='x', which='minor', alpha = 0.5 )
ax.grid( 'off', axis='y', alpha = 0.2)

figName = "outputPictures/outputFWHM" + str(int(filenumber)) + ".png"
plt.savefig(figName)

plt.show()
