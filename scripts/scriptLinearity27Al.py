#!/opt/anaconda3/bin/python3


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
                #print(countDetectorTotale)

                
            if(line.strip().endswith("Residue (keV)") and countDetectorTotale == detToAnalyse):
                lines = []
                append = True
                continue
                #print("appendTrue")
                
            if(append == True and countDetectorTotale == detToAnalyse):

                if(line.strip().startswith('#2')): 
                    lines.append(line.strip())
                    #print(line.strip())
                    
                runningCount +=1 
                
                if(line.strip().endswith('calibrated energy= 0')): 
                    lines = lines[:-1]
                    runningCount = 100

            if(runningCount == 100):
                print("Error in analysed detector - check inside input file", list_files[countFile],"  Or change detector")
                runningCount = 0
                lines = []
                append = False
                continue
            
            if(runningCount == 200):
                #print(f"Skipping a detector not calibrated")
                runningCount = 0
                lines = []
                append = False
                    
                
             
            if(runningCount > 20 and countDetectorTotale == detToAnalyse): 
                runningCount = 0
                append = False
                
                a = np.array(([np.fromstring(i[3:], dtype=float, sep=' ') for i in lines[0:]]) , dtype=object)[1:]            
                a = np.array([i.flatten() for i in a])
                
                df = pd.DataFrame( a ) 

                for ind in range(len(df.index)) :
                    if(df.iloc[ind, 5]< 1790 and df.iloc[ind, 5] > 1770):
                        df_1778.append([countFile, df.iloc[ind,2]])

                    if(df.iloc[ind, 5]< 7730 and df.iloc[ind, 5] > 7720):
                        df_7724.append([countFile, df.iloc[ind,2]])
                
                #plt.scatter(df[4], df[5], linewidth = None, label = f"Clover {countClover}, detector {countDetector}")     
                lines = []
                break
                
        countClover = 0
        countDetector = 0
        countDetectorTotale = 0
        continue

df_7724 = pd.DataFrame(df_7724)
df_1778 = pd.DataFrame(df_1778)


vecFIPPSFill = []
vecIFINFill  = []
                        
            
fig, axs = plt.subplots( 1,2, figsize = (15,5))


axs[0].scatter(x=df_1778[0], y = df_1778[1]/df_1778.iloc[0,1], label = "1778")
axs[0].scatter(x=df_7724[0], y = df_7724[1]/df_7724.iloc[0,1], label = "7724")
axs[0].set_title("Comparison of ratios for three peaks")
axs[0].set_xlabel("Run number (starting from zero)")
axs[0].set_ylabel("Ratio of channels energy for each peak to first run")


for i in vecFIPPSFill:
	axs[0].axvline(i, c = 'red', label = "FIPPS Fill" )
for i in vecIFINFill:
	axs[0].axvline(i, c = 'black', label = "IFIN Fill")

axs[0].legend()

axs[1].scatter(x=df_1778[0], y = df_1778[1]-df_1778.iloc[0,1], label = "778")
axs[1].scatter(x=df_7724[0], y = df_7724[1]-df_7724.iloc[0,1], label = "7724")

axs[1].set_title("Comparison of differences for three peaks")
axs[1].set_xlabel("Run number (starting from zero)")
axs[1].set_ylabel("Difference of channels energy for each peak to first run")



for i in vecFIPPSFill:
	axs[1].axvline(i, c = 'red', label = "FIPPS Fill" )
for i in vecIFINFill:
	axs[1].axvline(i, c = 'black', label = "IFIN Fill" )

axs[1].legend()


pString = str(runBegin)+"-"+str(runEnd)+"_det"+str(detToAnalyse)
fig.suptitle(pString)

pString = "outputPictures/outputLinearity/" + str(runBegin) + "-"+ str(runEnd) + "_det"+ str(detToAnalyse) + ".png"
plt.savefig(pString)
plt.show()


