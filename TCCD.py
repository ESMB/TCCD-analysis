#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 08:32:47 2020

@author: Mathew
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import os
pathlist=[]

# Where to store the overall file containing means etc. for each experiment.
path_root=r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/"

# Foldert to analyse here:
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/20pM488_1pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/PBS Only/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/30pM488_4pM647_200nMasyn/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/35pM488_5pM647_200nMasyn_1/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/30pM488_4pM647_20nMasyn_1/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/25pM488_1pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/70pM488_10pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/30pM488_4pM647_200nMasyn_7mins/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/20pM488_3pM647_200nMasyn/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/5pM488_1pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/35pM488_5pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/30pM488_4pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/30pM488_4pM647_20nMasyn/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/7pM488_1pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/10pM488_1pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/20pM488_3pM647/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20210712_UCBhighaffinityantibodies/Setup experiments/Concentration 50_50ratio/35pM488_5pM647_200nMasyn/")

file_stem="Ab"          # This is the part of the filename that will be searched for in each folder.
# number_of_files=5       # Number of files in the folder (could make this automatic in the future).

# Thresholds and other parameters:
    
channelA_thresh=10      # Threshold for Channel A (Green).
channelB_thresh=10      # Threshold for Channel B (Red).
channelA_AF=1.16        # Autofluorescence
channelB_AF=1.13
xtalk=0.237             # Cross-talk from A to B

def load_files(filename_contains,path):
    print(path)
    num=0
    channelA_sample=[]             # Where channel A data will be stored
    channelB_sample=[]             # Where channel B data will be stored
    for root, dirs, files in os.walk(path):
      for name in files:
              # print(name)
              if filename_contains in name:
                  resultsname = name
                  print(name)
                  num+=1
                  a=0
                  with open(path+name) as csvDataFile:                                                # Opens the file as a CSV
                        csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
                        for row in csvReader:
                            channelA_sample.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
                            channelB_sample.append(row[1])
                            a+=1
        
                        print ("Loaded %s, which contains %s rows."%(resultsname,a))
    rows=len(channelA_sample)
    print("Loaded %s files in total, with a total of %s rows"%(num,rows))
        
    channelA_arr_sample=np.asarray(channelA_sample,dtype=np.float32)                              # Converts these to numpy arrays for vector calcs.
    channelB_arr_sample=np.asarray(channelB_sample,dtype=np.float32)
    return channelA_arr_sample,channelB_arr_sample,num

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
    
Output_all = pd.DataFrame(columns=['Path','Number_of_files','Threshold_A','Threshold_B','Events_A','Events_B','Events_coincindent',
                                       'Events_chance','Q','Total_Intensity_mean','Total_Intensity_SD','Total_Intensity_med','Intensity_A_mean','Intensity_A_SD','Intensity_A_med','Intensity_B_mean','Intensity_B_SD','Intensity_B_med'])

for path in pathlist:
    channelA_arr,channelB_arr,num=load_files(file_stem,path)
    
    
    
    # Now need to account for autofluorescence and crosstalk etc. 
    
    channelB_arr=(channelB_arr-xtalk*channelA_arr)-channelB_AF
    channelA_arr=channelA_arr-channelA_AF
    
    #This part is for the thresholding:
    
    channelA_only_events=channelA_arr[(channelA_arr>channelA_thresh)]                       # Total A events
    channelB_only_events=channelB_arr[(channelB_arr>channelB_thresh)]                       # Total B events
    channelA_events=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # A coincident events             
    channelB_events=channelB_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # B coincident events
    
    channelA_mean=channelA_only_events.mean()
    channelA_SD=channelA_only_events.std()
    channelA_med=np.median(channelA_only_events)
    
    channelB_mean=channelB_only_events.mean()
    channelB_SD=channelB_only_events.std()
    channelB_med=np.median(channelB_only_events)
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
    
    
    channelB_arr_inv=channelB_arr*(-1)
    
    plt.plot(channelA_arr,color='green')
    plt.plot(channelB_arr_inv,color='red')
    plt.xlabel('Bin number')
    plt.ylabel('Intensity (photons)')
    plt.xlim(0,8000)
    plt.ylim(-200,200)
    plt.savefig(path+'/'+'example_trace.pdf')
    plt.show()
    
    total_intensity=channelB_events+channelA_events
    plt.hist(total_intensity, bins = 60,range=[0,1000], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,)
    plt.yscale('log')
    plt.xlabel('Total intensity (photons)')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'Intensity.pdf')
    plt.show()
    
    total_mean=total_intensity.mean()
    total_SD=total_intensity.std()
    total_med=np.median(total_intensity)
    Output_all= Output_all.append({'Path':path,'Number_of_files':num,'Threshold_A':channelA_thresh,'Threshold_B':channelB_thresh,'Events_A':var_A_events,'Events_B':var_B_events,'Events_coincindent':var_real_events,'Q':Q,
                                       'Events_chance':var_chance_events,'Total_Intensity_mean':total_mean,'Total_Intensity_SD':total_SD,'Totla_Intensity_med':total_med,'Intensity_A_mean':channelA_mean,'Intensity_A_SD':channelA_SD,'Intensity_A_med':channelA_med,'Intensity_B_mean':channelB_mean,'Intensity_B_SD':channelB_SD,'Intensity_B_med':channelB_med},ignore_index=True)

Output_all.to_csv(path_root + '/' + 'all_metrics.csv', sep = '\t')


