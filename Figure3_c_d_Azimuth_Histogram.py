#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
from scipy.signal import coherence
from statistics import median as median
import sys
from scipy.optimize import fmin
import numpy as np
from matplotlib.mlab import csd
from math import pi
import sys
from obspy.signal.invsim import evalresp
import matplotlib.pyplot as plt
from obspy.signal.rotate import rotate2zne
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
from multiprocessing import Pool 
import csv


# Histogram of best fit orientations of a test sensor with a reference using sliding data windows. 

debug = True


def RobrotData(angle, trR, st, nol, nwin, debug = False):
    
    # Band which to maximize coherence over 
    minper = 4.
    maxper = 8. 
    
    stTemp = st.copy()
    stTemp.rotate('NE->RT',angle)
    trR2 = stTemp.select(component='R')[0]
    (fre, cxy) =  coherence(trR.data,trR2.data, \
            fs = 1./trR.stats.delta, noverlap = int(nol*nwin), nfft = nwin, nperseg = nwin)
    fre = fre[1:]
    cxy = cxy[1:] 
    mask = (fre <= 1./minper) & (fre >= 1./maxper)
    newcoh = np.abs(np.mean(cxy[mask])-1.)
    # print(newcoh)
    return newcoh
    


## Start of Actual Code #################################################


# Enter info for reference station (aligned to North)
net = "GS"
sta = "ALQ1"
loc = "00"
chan = "LH*"


# Test station
net2 = "GS"
sta2 = "ALQ2"
loc2 = "00"
# initial guess of station's orientaion (use Sensor Test Suite)
Or_guess = 180.


# Pressure Data
netp = "IU"
stap = "ANMO"
loc3 = "50"
chan3 = "LDO"


## Windowing parameters 

Window = 3600.
Overlap = 0.5

## Parameters to calculate coherance within windows
Coh_Length = 2**9
Coh_Overlap = 0.75

# Data completeness for each window
Comp_Thres = 0.9


# Set the start day - we will set various start dates in the window below
sday = UTCDateTime("2019-001T00:00:00")


# read in several months of data

st = Stream()
st2 = Stream()
stP = Stream()
Daynum = 0.

# Make a stream for all of the reference station data

for day in range(239,340):
    try:
        st += read("/msd/" + net + "_" + sta + "/" + str(sday.year) + "/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
    except: 
        print('Day',day, 'has no data for reference sensor')    
    Daynum=Daynum+1.
    
    
# Make a stream for all of the test station data
for day in range(239,340):
    try:
        st2 += read("/msd/" + net2 + "_" + sta2 + "/" + str(sday.year) + "/" + str(day).zfill(3) + "/" + loc2 + "_" + chan + "*.seed")
    except: 
        print('Day',day, 'has no data for test sensor') 
print(Daynum)


for day in range(239,340):
    try:
        stP += read("/msd/" + netp + "_" + stap + "/" + str(sday.year) + "/" + str(day).zfill(3) + "/" + loc3 + "_" + chan3 + "*.seed")
    except: 
        print('Day',day, 'has no data for test sensor') 
print(Daynum)



# Rename data for rotations and organize in reverse 
for tr in st:
    if tr.stats.channel == 'LH1':
        tr.stats.channel = 'LHN'
    if tr.stats.channel == 'LH2':
        tr.stats.channel = 'LHE'
st.detrend()        
st.merge()

st.sort(reverse=True)


for tr in st2:
    if tr.stats.channel == 'LH1':
        tr.stats.channel = 'LHN'
    if tr.stats.channel == 'LH2':
        tr.stats.channel = 'LHE'
st2.detrend()        
st2.merge(fill_value=0.)
st2.sort(reverse=True)

# Data is in ZNE format
if debug:
    print(st)
    print(st2)

# Wrap-around issues near 0 degree. Hack and rotate test station 180

st2 = st2.rotate('NE->RT',180.)
st2 = st2.rotate('RT->NE',0.)



# Setup the variable arrays  
angF = [] 
coh = []
Pres = []

syear = []
smonth = []
sday = []
shour = []
sminute = []
ssecond = []

fs = st2[0].stats.sampling_rate
if debug:
    print("Sample rate is ", fs)
total_count = Window*fs

# Start of Sliding window loop to calculate orientation angle
for stW in st.slide(Window, (1.-Overlap)*Window):
    print('On day: ' + str(stW[0].stats.starttime))
    
    # Get start and end times and trim accordingly 
    starttime = stW[0].stats.starttime
    
    endtime = starttime+Window
    
    st2W = st2.copy()
    st2W.trim(starttime,endtime)
    
    
    stPW = stP.copy()
    stPW.trim(starttime,endtime)
    
    # Fill data gaps and detrend 
    stW.merge(fill_value=0.)
    st2W.merge(fill_value=0.)
    stPW.merge(fill_value=0.)
    
    data_count = np.count_nonzero(st2W[0].data)
    if data_count/total_count >= Comp_Thres:
    
        try:
        
            stW.detrend
            st2W.detrend
        
            # Convert from counts to Pascals 
            Pressure = stPW[0].data*10.
            
            # get angle and coherence, append pressure 
            ang = fmin(RobrotData, Or_guess, args=(stW.select(component='N')[0], st2W, Coh_Overlap, Coh_Length))
            angF.append(ang)
            coh.append(RobrotData(ang, stW.select(component='N')[0], st2W, Coh_Overlap,Coh_Length))
            Pres.append(median(Pressure))
            
            # save start and end times
            syear.append(starttime.year)
            smonth.append(starttime.month)
            sday.append(starttime.day)
            shour.append(starttime.hour)
            sminute.append(starttime.minute)
            ssecond.append(starttime.second)
        
        except Exception as e:
            print(e)
 
 
# Make all variables we wish to store into arrays            
angF = 360-np.asarray(angF) 
coh = 1.-np.asarray(coh)
Pres = np.asarray(Pres)

np.asarray(syear)
np.asarray(smonth)
np.asarray(sday)
np.asarray(shour)
np.asarray(sminute)
np.asarray(ssecond)



# Concatinate the variables and save
Output = np.c_[syear,smonth,sday,shour,sminute,ssecond,angF,coh,Pres]

print(Output)

f3=open('Orientation_ALQ2_SM_Pres_GS_2019.txt','w')

writer = csv.writer(f3, delimiter='\t')
writer.writerows(Output)
f3.close()
      
  
fig = plt.figure(1)
plt.hist(angF,100, range=(175,185),density=True,facecolor='0.6')
plt.ylabel('Normalized Counts')
plt.xlabel('Orientation Angle')
plt.legend()

fig = plt.figure(2)
plt.scatter(angF,coh)
plt.ylabel('Coherance')
plt.xlabel('Orientation Angle')
plt.xlim((177.,183.))

plt.show()  
