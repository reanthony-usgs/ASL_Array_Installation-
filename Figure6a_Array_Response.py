#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from obspy import read_inventory, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.imaging.cm import obspy_sequential
from obspy.signal.array_analysis import array_transff_wavenumber

def pol2cart(theta, rho):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y

debug = False 

# Make a list of all the stations in the Array 
stas = ['ASA1', 'ASA2', 'ASA3', 'ASA4', 'ASA5', 'ASA6', 'ALQ1', 'ASL9','ANMO']


if debug:
    stas = ['ASA1', 'ASA2', 'ASA3']

chan = 'HHZ'
starttime = UTCDateTime(2019,10,1)
# get the dataless to obtain the lats, lons, and elevations 
lats, lons, elevs = [],  [], []

for sta in stas:
    if sta == 'ANMO':
        net, loc = 'IU', '10'
    else:
        net, loc = 'GS', '00'
    

    else:

        # read in the dataless for each station as an inventory object    
        path = '/APPS/metadata/SEED/' + net + '.' + sta + '.dataless'
        inv = read_inventory(path)
    
        # scrub all the excess informaion 
        inv_chan = inv.select(station=sta,location=loc,channel = chan, starttime=starttime)
    
        lats.append(inv_chan[0][0][0].latitude)
        lons.append(inv_chan[0][0][0].longitude)
        elevs.append(inv_chan[0][0][0].elevation)

# Now Convert to NP array   
lats = np.asarray(lats)
lons = np.asarray(lons)
elevs = np.asarray(elevs)


# generate array coordinates
coords = np.array([lons,lats,elevs])
coords = np.transpose(coords)

print(coords)

# coordinates in km
#coords /= 1000.

# set limits for wavenumber differences to analyze
klim = 30.
kxmin = -klim
kxmax = klim
kymin = -klim
kymax = klim
kstep = klim / 200.

# compute transfer function as a function of wavenumber difference

# transfer function is in units of power - so convert to dB using 10*log(10)
transff = 10*np.log10(array_transff_wavenumber(coords, klim, kstep, coordsys='lonlat'))

# plot

dKX =  np.arange(kxmin, kxmax + kstep , kstep) 
dKY =  np.arange(kymin, kymax + kstep , kstep) 

X_Grid, Y_Grid = np.meshgrid(dKX, dKY)

if debug:
    print(np.shape(dKX), np.shape(transff))

plt.pcolor(dKX, dKY, transff.T, cmap=obspy_sequential)   
plt.tick_params(length=6, width=3, labelsize=12)     

cbar = plt.colorbar()
cbar.set_label('Relative Beam Power (dB)', rotation=270, fontsize=18, labelpad=9)
cbar.ax.tick_params(length=6, width=3, labelsize=14) 
plt.clim(vmin=-10., vmax=0.)

plt.contour(X_Grid,Y_Grid,transff.T, [-9.,-6.,-3.], colors='w')   
plt.xlim(kxmin, kxmax)
plt.ylim(kymin, kymax)
plt.xlabel('E-W Wavenumber Difference',fontsize=18)
plt.ylabel('N-S Wavenumber Difference',fontsize=18)


plt.show()
