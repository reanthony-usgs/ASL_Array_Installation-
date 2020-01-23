#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
from obspy.io.xseed import Parser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Code to plot the seismogram for the local event 

debug = True 

# User inputs 

net = "GS"
sta = "ALQ1"
loc = "00"
chan = "HHZ"

starttime = (2019,10,1)
stime=UTCDateTime("2019-253T10:00:00")
etime = stime + 12.*3600

if debug:
    print(stime.day)

# Load in the data and get just the station you want 
st = read("/msd/" + net + "_" + sta + "/2019/" + str(stime.julday).zfill(3) + "/" + loc + "_" + chan + "*.seed")


if debug:
    print(st)

tr = st[0]
tr.detrend('constant')


path = '/APPS/metadata/SEED/' + net + '.' + sta + '.dataless'    
parser = Parser(path)

tr.simulate(seedresp={'filename':parser, 'units':"VEL"})
tr.filter("bandpass",freqmin = 3.1 ,freqmax= 3.2, corners=4.)
    


# Plotting format
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

fig = plt.figure(2)
tr.trim(stime, etime)
t=np.arange(0, len(tr.data))/(tr.stats.sampling_rate*3600.)
plt.plot(t,tr.data*1000000,linewidth=1,c='k')
plt.xlabel('Time (Hrs)')
plt.ylabel('Velocity (mm/s)')
plt.xlim((0, 12))
plt.show()
