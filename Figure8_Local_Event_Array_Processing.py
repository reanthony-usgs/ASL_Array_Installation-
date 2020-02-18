#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize


import obspy
from obspy.io.xseed import Parser
from obspy import read_inventory, UTCDateTime
from obspy.core.util import AttribDict
from obspy.imaging.cm import obspy_sequential
from obspy.signal.invsim import corn_freq_2_paz
from obspy.signal.array_analysis import array_processing

# Import/apply font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)


# Load data
st = obspy.read("Array_Sep_Event.mseed")

debug = True 

# Global parameters 

starttime = (2019,10,1)
stime=obspy.UTCDateTime("2019-253T18:11:15")
etime = stime + 20.



# Instrument correction to 1Hz corner frequency
# paz1hz = corn_freq_2_paz(1.0, damp=0.707)

# Now we need to loop through all of the channels in the stream,
# grab the dataless, and setup dictionaties for each channel

for tr in st: 
    net = tr.stats.network
    sta = tr.stats.station 
    loc = tr.stats.location 
    chan = tr.stats.channel 
    
    if debug:
        print('sample rate of', sta, 'is:', tr.stats.sampling_rate)
    
    # read in the dataless for each station as an inventory object    
    path = '/APPS/metadata/SEED/' + net + '.' + sta + '.dataless'
    
    inv = read_inventory(path)
    
    # scrub all the excess informaion 
    inv_chan = inv.select(station=sta,location=loc,channel=chan, starttime=starttime)
    
    #  ALQ1 for plotting
    if sta == 'ALQ1':
        
        # Prepare ALQ1 trace for plotting 
        tr_plot = tr.copy()
        tr_plot.detrend('constant')
        parser = Parser(path)
    
        #remove response from trace
        tr_plot.simulate(seedresp={'filename':parser, 'units':"VEL"})
        tr_plot.filter('bandpass',freqmin=1.0,freqmax=15.,corners=4)
        
        
    else:
    
        # Get the desired parameters and throw into dictionaries 
        tr.stats.coordinates = AttribDict({
        'latitude': inv_chan[0][0][0].latitude,
        'elevation': inv_chan[0][0][0].elevation/1000.,
        'longitude': inv_chan[0][0][0].longitude})
    
    parser = Parser(path)
    
    #remove response from trace
    tr.simulate(seedresp={'filename':parser, 'units':"VEL"})
    tr.filter('bandpass',freqmin=0.5,freqmax=15.)
    
    
# Execute array_processing
kwargs = dict(
    # slowness grid: X min, X max, Y min, Y max, Slow Step
    sll_x=-3.0, slm_x=3.0, sll_y=-3.0, slm_y=3.0, sl_s=0.03,
    # sliding window properties
    win_len=0.5, win_frac=0.1,
    # frequency properties
    frqlow=0.5, frqhigh=15.0, prewhiten=0,
    # restrict output
    semb_thres=-1e9, vel_thres=-1e9,
    stime=stime,etime=etime)
    
out = array_processing(st, **kwargs)

# Plot

# prepare timeseries 



cmap = obspy_sequential

# make output human readable, adjust backazimuth to values between 0 and 360
t, rel_power, abs_power, baz, slow = out.T
baz[baz < 0.0] += 360

t = t-t[0]
t *= 3600.*24.

if debug:
    print(t)


# Creating figure and moving subplots to have no space between them
fig = plt.figure(1, figsize=(12,16))
plt.subplots_adjust(hspace=0.001)

# Define labels for the legend
labels = ['Ground Velocity (mm/s)', 'Relative Power', 'Backazimuth ($^\circ$)', 'Apparent\nVelocity (km/s)']

# Making time vector 
tr_plot.trim(stime, etime)
times = range(tr_plot.stats.npts)
times = np.asarray(times)/(tr_plot.stats.sampling_rate) 



# First subplot (bandpass filtered Seismogram) 
plt.subplot(411)
plt.plot(times,tr_plot.data*1000000.,linewidth=1,c='k')
plt.ylabel('Ground\nVelocity ($\mu$m/s)', fontsize=18)
plt.tick_params(labelsize = 15) 
plt.text(-3,2.5,'(a)', fontsize=20)
plt.ylim(-2.5, 2.5)
plt.xlim(0, 20)
plt.xticks([])

# 2nd Subplot gives relative power
plt.subplot(412)
plt.scatter(t,rel_power, c=rel_power, alpha=0.8,
               edgecolors='none', cmap=obspy_sequential)
plt.ylabel('Relative Power', fontsize=18)  
plt.tick_params(labelsize = 15) 
plt.text(-3,1.0,'(b)', fontsize=20)
plt.xlim(0, 20)
plt.xticks([])

# 3rd plot gives backazimuth 
plt.subplot(413)
plt.scatter(t,baz, c=rel_power, alpha=0.8,
               edgecolors='none', cmap=obspy_sequential)
plt.ylabel('Back Azimuth ($^\circ$)', fontsize=18) 
plt.ylim(0, 360)  
plt.xlim(0, 20)
plt.xticks([])
plt.yticks([90,180,270])
plt.tick_params(labelsize = 15) 
plt.text(-3,360.0,'(c)', fontsize=20)
plt.text(-3,0,'(d)', fontsize=20)

# 4th plot gives apparent velocity  
plt.subplot(414)

plt.scatter(t,1./slow, c=rel_power, alpha=0.8,
               edgecolors='none', cmap=obspy_sequential)
plt.ylabel('Apparent\nVelocity (km/s)', fontsize=18)
plt.ylim(2, 10.0)  

             


# Horizontal axis labeling
plt.xlabel('Time (s)', fontsize = 18)  
plt.tick_params(labelsize = 15)  
plt.xlim(0, 20)

# Save or show figure
plt.show()













   
