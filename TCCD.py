#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 08:32:47 2020

@author: Mathew
"""

import numpy as np
import pandas as pd
import csv


# The thresholds and filenames etc. need to be placed below

path = r"/Users/Mathew/Dropbox (Cambridge University)/Ed Code/TCCD Analysis/Files/"           # This is the folder that needs to be analysed

file_stem="DNA"         # This is the filename before the underscore. 
number_of_files=2       # Number of files in the folder (could make this automatic in the future).
channelA_thresh=10      # Threshold for Channel A (Green).
channelB_thresh=10      # Threshold for Channel B (Red).
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
        
        print "Loaded %s, which contains %s rows."%(filename,a)   
        
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

print ('There were %s A events, %s B events, %s coincidence events, and %s chance events. Q = %f.')%(var_A_events,var_B_events,var_real_events,var_chance_events,Q)

