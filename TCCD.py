#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 08:32:47 2020

@author: Mathew
"""

import numpy as np
import matplotlib.pyplot as plt
import csv



# The thresholds and filenames etc. need to be placed below

path = r"/Users/Mathew/Documents/Edinburgh Code/TCCD Analysis/Files/"           # This is the folder that needs to be analysed

file_stem="DNA"         # This is the filename before the underscore. 
number_of_files=2       # Number of files in the folder (could make this automatic in the future).
channelA_thresh=5      # Threshold for Channel A (Green).
channelB_thresh=12      # Threshold for Channel B (Red).
channelA_AF=1.16        # Autofluorescence
channelB_AF=1.13
xtalk=0.237             # Cross-talk from A to B

def load_files(number_of_files):
    
    channelA=[]             # Where channel A data will be stored
    channelB=[]             # Where channel B data will be stored
    for i in range(0,number_of_files):
        
        if(i<1):
            filename=path+file_stem+".txt"
        elif(i<10):
            filename=path+file_stem+"_0"+str(i+1)+".txt"
        else:
            filename=path+file_stem+"_"+str(i+1)+".txt"
        a=0                                                                             # Row counter
        with open(filename) as csvDataFile:                                                # Opens the file as a CSV
            csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
            for row in csvReader:
                channelA.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
                channelB.append(row[1])
                a+=1
        
        print("Loaded %s, which contains %s rows."%(filename,a))   
        
    channelA_arr=np.asarray(channelA,dtype=np.float32)                              # Converts these to numpy arrays for vector calcs.
    channelB_arr=np.asarray(channelB,dtype=np.float32)
    return channelA_arr,channelB_arr

channelA_arr,channelB_arr=load_files(number_of_files)



# Now need to account for autofluorescence and crosstalk etc. 

channelB_arr=(channelB_arr-xtalk*channelA_arr)-channelB_AF
channelA_arr=channelA_arr-channelA_AF

#This part is for the thresholding:

channelA_only_events=channelA_arr[(channelA_arr>channelA_thresh)]                       # Total A events
channelB_only_events=channelB_arr[(channelB_arr>channelB_thresh)]                       # Total B events
channelA_events=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # A coincident events             
channelB_events=channelB_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # B coincident events

# Now need to account for chance events:

channelB_shuffle=channelB_arr.copy()
np.random.shuffle(channelB_shuffle)

channelA_chance=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_shuffle>channelB_thresh)]    # These are the chance events    
channelB_chance=channelB_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_shuffle>channelB_thresh)]

# Now need to calculate Q value:

var_real_events=float(len(channelA_events))
var_A_events=float(len(channelA_only_events))
var_B_events=float(len(channelB_only_events))
var_chance_events=float(len(channelA_chance))
Q=float((var_real_events-var_chance_events)/(var_A_events+var_B_events-(var_real_events-var_chance_events)))

print(('There were %s A events, %s B events, %s coincidence events, and %s chance events. Q = %f.')%(var_A_events,var_B_events,var_real_events,var_chance_events,Q))


# Now want histograms etc. 

ln_events=np.log(channelB_events/channelA_events)
ln_chance=np.log(channelB_chance/channelA_chance)

textstr='Q = %.3f'%Q

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = "12"
plt.figure(figsize=(8, 6))
plt.hist(ln_events, bins = 60,range=[-3,3], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,label="Real Events")
plt.hist(ln_chance, bins = 60,range=[-3,3], rwidth=0.9,ec='black',color='#cccccc',alpha=0.5,label="Chance Events")
plt.text(0.05,0.90, textstr,transform=plt.gca().transAxes)
plt.legend(loc='upper right')         
plt.xlabel('Z=$ln(I_{B}/I_{A}$)')
plt.ylabel('Number of events')
plt.savefig(path+'/'+'lnz.pdf')
plt.show()



def maxQ():
    q_vals = np.zeros(shape=(20,20))

    for A in range(20):
        for B in range(20):
            channelA_only_events=channelA_arr[(channelA_arr>A)]                       # Total A events
            channelB_only_events=channelB_arr[(channelB_arr>B)]                       # Total B events
            channelA_events=channelA_arr[np.logical_and(channelA_arr>A, channelB_arr>B)]  # A coincident events             
            
            
            # Now need to account for chance events:
            
            channelB_shuffle=channelB_arr.copy()
            np.random.shuffle(channelB_shuffle)
            
            channelA_chance=channelA_arr[np.logical_and(channelA_arr>A, channelB_shuffle>B)]    # These are the chance events    
            
            # Now need to calculate Q value:
            
            var_real_events=float(len(channelA_events))
            var_A_events=float(len(channelA_only_events))
            var_B_events=float(len(channelB_only_events))
            var_chance_events=float(len(channelA_chance))
            Q=float((var_real_events-var_chance_events)/(var_A_events+var_B_events-(var_real_events-var_chance_events)))
    
            q_vals[A][B]=Q
            

    maximum_Q=np.amax(q_vals)
    result=np.where(q_vals == np.amax(q_vals))
    ThresholdA,ThresholdB=result
    
    print('The maximum value of Q is %.3f, with a threshold of %s in channel A, and %s in channel B.'%(maximum_Q,str(ThresholdA),ThresholdB))
    
    
    contourplot = plt.contourf(q_vals,20,origin='lower')
    cbar = plt.colorbar(contourplot)
    plt.xlabel("Channel B Threshold")
    plt.ylabel("Channel A Threshold")
    cbar.ax.set_ylabel('Q')

