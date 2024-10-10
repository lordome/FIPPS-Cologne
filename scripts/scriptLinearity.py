#!/opt/anaconda3/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if( len(sys.argv) != 4 ): 
    print("Usage: scriptLinearity.py runBegin runEnd detToAnalyse")
    exit()

runBegin = int(sys.argv[1])
runEnd = int(sys.argv[2])
detToAnalyse = int(sys.argv[3])

append = False
runningCount = 0
lines = []

countClover = 0
countDetector = 0
countDetectorTotale = 0

df = np.zeros((6,10))

list_files = [str(i) + ".txt" for i in range(runBegin, runEnd+1)]

list_files = ["outputAnaRaw/"+i for i in list_files]

df_121 = []
df_778 = []
df_1408 = []

for countFile, iterfiles in enumerate(list_files):
        
    with open(iterfiles, 'r') as fp:


        for i, line in enumerate(fp):

            if(line.strip().endswith('calibrated energy= 0')):   
                    runningCount = 200
                    
            if(line.strip().endswith("rChi2%")):
                countDetector += 1 
                countDetectorTotale += 1
                #print(countDetectorTotale)

                
            if(line.strip().endswith("Residue (keV)") and countDetectorTotale == detToAnalyse):
                lines = []
                append = True
                #print("appendTrue")
                
            if(append == True and countDetectorTotale == detToAnalyse): 
                lines.append(line.strip())
                runningCount +=1 
                
                if(line.strip().endswith('calibrated energy= 0')): 
                    lines = lines[:-1]
                    runningCount = 100

            if(runningCount == 100):
                print(f"Error in analysed detector - check inside input file {list_files[countFile]}. Or change detector")
                runningCount = 0
                lines = []
                append = False
                continue
            
            if(runningCount == 200):
                #print(f"Skipping a detector not calibrated")
                runningCount = 0
                lines = []
                append = False
                    
                
             
            if(runningCount > 10 and countDetectorTotale == detToAnalyse): 
                runningCount = 0
                append = False
                
                a = np.array(([np.fromstring(i[3:], dtype=float, sep=' ') for i in lines[0:]]) , dtype=object)[1:]            
                a = np.array([i.flatten() for i in a])
                
                if ( a.shape != (10,) and a[0].shape!=(6,) ):
                    a = np.zeros((10,6))
                    a[:,4] = df[4]
                df = pd.DataFrame( a ) 
                #print(df)
	
                
                #plt.scatter(df[4], df[5], linewidth = None, label = f"Clover {countClover}, detector {countDetector}")
                df_121.append([countFile, df.iloc[0,1]]) 
                df_778.append([countFile, df.iloc[5,1]])         
                df_1408.append([countFile, df.iloc[9,1]])         
                lines = []
                break
                
        countClover = 0
        countDetector = 0
        countDetectorTotale = 0
        continue

df_121 = pd.DataFrame(df_121)
df_778 = pd.DataFrame(df_778)
df_1408 = pd.DataFrame(df_1408)


#print(df_121)

vecFIPPSFill = [39, 132]
vecIFINFill  = [53, 148]
                        
            
fig, axs = plt.subplots( 1,2, figsize = (15,5))


axs[0].scatter(x=df_121[0], y = df_121[1]/df_121.iloc[0,1], label = f"121")
axs[0].scatter(x=df_778[0], y = df_778[1]/df_778.iloc[0,1], label = f"778")
axs[0].scatter(x=df_1408[0], y = df_1408[1]/df_1408.iloc[0,1], label = f"1408")

axs[0].set_title("Comparison of ratios for three peaks")
axs[0].set_xlabel("Run number (starting from zero)")
axs[0].set_ylabel("Ratio of channels energy for each peak to first run")



for i in vecFIPPSFill:
	axs[0].axvline(i, c = 'red', label = "FIPPS Fill" )
for i in vecIFINFill:
	axs[0].axvline(i, c = 'black', label = "IFIN Fill")

axs[0].legend()

axs[1].scatter(x=df_121[0], y = df_121[1]-df_121.iloc[0,1], label = f"121")
axs[1].scatter(x=df_778[0], y = df_778[1]-df_778.iloc[0,1], label = f"778")
axs[1].scatter(x=df_1408[0], y = df_1408[1]-df_1408.iloc[0,1], label = f"1408")

axs[1].set_title("Comparison of differences for three peaks")
axs[1].set_xlabel("Run number (starting from zero)")
axs[1].set_ylabel("Difference of channels energy for each peak to first run")



for i in vecFIPPSFill:
	axs[1].axvline(i, c = 'red', label = "FIPPS Fill" )
for i in vecIFINFill:
	axs[1].axvline(i, c = 'black', label = "IFIN Fill" )

axs[1].legend()

fig.suptitle(f"{runBegin}-{runEnd}_det{detToAnalyse}")
plt.savefig(f"outputPictures/outputLinearity/{runBegin}-{runEnd}_det{detToAnalyse}.png")
#plt.show()


