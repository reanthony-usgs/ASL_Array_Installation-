#!/usr/bin/env python


# This script loads in a day-long PSD and calculates the peridogram
# Estimate for several randomly selected windows. 

# Confidence intervals are then assigned based on the distribution of 
# these peridogram estimates in various bins.



from obspy.core import Stream, read, UTCDateTime, inventory 
from obspy.core.inventory.inventory import read_inventory 
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
from scipy.signal import welch
from obspy.signal.invsim import evalresp
from math import pi
import sys
from itertools import combinations
import math, copy
import random
import csv


debug = True

net = "GS"
sta = "VEA1"
loc = "00"
chan = "HH2"



# Parameters for Code ##################################################

# input Sample Rate
fs = 100

# subwindow size and overlap
nfft = 2**16
windlap = 0.75 
# size of window and overlap for PSD
window = 3600
overlap = 0.5

# Data completeness for each window
Comp_Thres = 0.9

########################################################################

seed_ID = net + "." + sta + "." + loc + "." + chan
sday = UTCDateTime("2019-238T00:00:00")


# Get the response

resppath = "/APPS/metadata/RESPS/RESP." + net + "." + sta + "." + loc + "." + chan
inv = read_inventory(resppath)
inv_Resp = inv.get_response(seed_ID,sday)
resp, freqR = inv_Resp.get_evalresp_response(t_samp=1./fs, nfft=nfft, output='ACC')
resp = resp[1:]

# read in several months of data
st = Stream()

#for day in range(99,235):
#    try:
#        st += read("/msd/" + net + "_" + sta + "/" + str(sday.year) + "/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
#    except:
#        print("No Data Here")
    
# Add in GS Data   
    
for day in range(238,263):
    st += read("/msd/" + net + "_" + sta + "/2019/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
    

if debug:
    print(st)
    
st.detrend('constant')

#st.merge(fill_value=0.)


# Start the loop of PSD calculations 

if debug:
    print('Heading into the loop')

total_count = window*fs

for stT in st.slide(window,window*(1.-overlap)):
    
    stT.merge(fill_value=0.)
    data_count = np.count_nonzero(stT[0].data)
        
    if data_count/total_count >= Comp_Thres:
        try:
        
            tr = stT[0]
        
            # Get all timing information 
            
            year = tr.stats.starttime.year
            julday = tr.stats.starttime.julday
            hour = tr.stats.starttime.hour
            minute = tr.stats.starttime.minute
            second = tr.stats.starttime.second
            
            if hour==10:
                print(julday)
            
            # Default uses a Hann Window and demeans the data
            freq, power = welch(tr.data, fs=fs, nperseg = nfft, 
                            noverlap = nfft*windlap, detrend = 'linear')
            
            freq = freq [1:]
            power = power[1:] 
            
            power = np.abs(power)
            power = 10.*np.log10(power/(np.abs(resp)**2))
            
            # get rid of the extra digits
            power = np.round(power,3)
            
            
             # Assemble timing information Vector
            date_vec = [year, julday, hour, minute, second]
        
            # Assemble the Matrix of PSDs and dates
        
        
            if "Power_M" not in vars():
                Power_M = power
                date_vec_M = date_vec
                
            else:
                Power_M = np.vstack((Power_M,power))
                date_vec_M = np.vstack((date_vec_M,date_vec))
        except:
            print('Not Enough Data')
    
# Flat File 
#freq = np.matrix.transpose(freq)

#f=open('HHZ_freqs_1Hr.txt','w')
#for fw in freq:
#    f.write(str(fw) + '\n')
#f.close()

f2=open('PSD_Times_VEA1_HH2_Month.txt','w')
writer = csv.writer(f2, delimiter='\t')
writer.writerows(date_vec_M)
f2.close()

f3=open('PSDs_VEA1_HH2_Month.txt','w')
writer = csv.writer(f3, delimiter='\t')
writer.writerows(Power_M)
f3.close()
    
# get the mean and STD

Power_Avg = np.mean(Power_M,axis=0)


print(Power_Avg.shape)        


# Plot up the mean PSD and 95% confidence intervals (2*STD)

NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()


fig = plt.figure(1)

plt.semilogx(1./freq, Power_Avg, linewidth=5., color='k', label='Mean PSD')
plt.semilogx(NLNMper, NLNMpower, linewidth=2., color='k')
plt.semilogx(NHNMper, NHNMpower, linewidth=2., color='k',label='NLNM/NHNM')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB)')
plt.xlim((0.05,500.))
plt.ylim((-200., -80.))
plt.legend()
plt.show()                                                                    
