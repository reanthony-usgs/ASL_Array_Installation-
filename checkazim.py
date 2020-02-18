#!/usr/bin/env python
from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.signal.polarization import flinn, particle_motion_odr
import numpy as np
from sklearn.decomposition import PCA




model = TauPyModel(model="iasp91")
client = Client("IRIS")


def decomppca(st, phase, debug = True):
    X = np.stack((st.select(component="Z")[0].data, st.select(component="R")[0].data), axis=0)
    if debug:
        print('Here is the covariance matrix:' + str(np.cov(X)))
    eig_vals, eig_vecs = np.linalg.eig(np.cov(X))
    if debug:
        print(eig_vals)
        print(eig_vecs)
    eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]
    # Sort from high to low
    eig_pairs.sort(key=lambda x: x[0], reverse=True)
    if debug:
        print(eig_pairs)
    if phase == 'P':
        angle = np.rad2deg(np.arccos(eig_pairs[0][1][0]/np.sqrt(np.sum(eig_pairs[0][1]**2))))
    else:
        angle = np.rad2deg(np.arccos(eig_pairs[1][1][0]/np.sqrt(np.sum(eig_pairs[1][1]**2))))
    scale = eig_pairs[0][0]/(eig_pairs[0][0]+eig_pairs[1][0])   
    if angle > 90.:
        angle = 180. - angle 
    return angle, scale


def getangles(net, sta, loc, stime, etime, debug = True):
    # Input net, station, location, starttime, endtime
    # grab the station info, then the events
    # We loop through the events and do the calculations on each of the events for the station
    # we append the results to the mess file and finally print it out
    
    mess = []
    # Get the latitude and longitude
    inv = client.get_stations(network=net, station=sta,
                                    channel = 'BH*', level="response",
                                    location=loc, starttime=stime, endtime=etime)
                                    
    stalat = inv[0][0][0].latitude
    stalon = inv[0][0][0].longitude
    staele = inv[0][0][0].elevation
    if debug:
        print('Station lat:' + str(stalat))
        print('Station lon:' + str(stalon))
        print('Station ele:' + str(staele))
        
    # Get our list of events
    cat = client.get_events(starttime=stime, minmagnitude=6., latitude=stalat, 
                        longitude=stalon, maxradius=90., minradius = 30., endtime=etime, mindepth=30)
                        
    for idx, eve in enumerate(cat):
        print('On event: ' + str(idx+1) + ' of ' + str(len(cat)))
        (dis, bazi, azi) = gps2dist_azimuth(stalat, stalon, eve.origins[0].latitude,eve.origins[0].longitude)
        
        # bazi is the station to event azimuth (back azimuth)
        # azi is the azimuth from the event to the station 
        if debug:
            print('Here is the back-azimuth: ' + str(bazi))
            print('Here is the azimuth: ' + str(azi))
        dis = kilometer2degrees(dis/1000.)
        arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis, phase_list=['P', 'S'])
        # If we have multiple phase skip to the next event
        if len(arrivals) != 2:
            continue
        if debug:
            print(arrivals)
            print(eve)
        
        # We need to also get a Noise window
        phasestime = (eve.origins[0].time + arrivals[0].time) - 10. 
        phaseetime = (eve.origins[0].time + arrivals[0].time) -5. 
        try:
            st = client.get_waveforms(net, sta, loc, 'BH*', phasestime,
                                    phaseetime)
        except:
            continue
        st.detrend('constant')

        st.remove_sensitivity(inventory=inv)


        if debug:
            print(phasestime)
            print(phaseetime)

        # from here we can estimate the SNR
        noise = sum(st.std())
        if debug:
            print('Here is the noise:' + str(noise))
        # So now we have two events one P and one S
        for arrival in arrivals:
            phasestime = (eve.origins[0].time + arrival.time) +0. 
            phaseetime = (eve.origins[0].time + arrival.time) +5. 
            if debug:
                print(phasestime)
                print(phaseetime)
            try:
                st = client.get_waveforms(net, sta, loc, 'BH*', phasestime,
                                      phaseetime)
            except:
                continue
            st.detrend('constant')
            #st.detrend('linear')
            #st.merge(fill_value=0)
            st.remove_sensitivity(inventory=inv)

            #st.taper(0.05)
            st.rotate(method="->ZNE", inventory=inv)
            
            signal = sum(st.std())
            if debug:
                print('Here is the signal:' + str(signal))
                print(st)
                print('Here is the SNR:' + str(signal/noise))
            azies, inces, aziesE, incesE = particle_motion_odr(st)
            if debug:
                print('Here is the azies: ' + str(azies) + ' +/-' + str(aziesE))
                print('Here is the inces: ' + str(inces) + ' +/-' + str(inces))
            
            st.rotate(method="NE->RT", back_azimuth=bazi)
            angle, scale  = decomppca(st, arrival.name)

            if debug:
                print('Here is PCA: ' + str(angle) + ' Here is odr:' + str(azies) )
            if debug:
                print('Here is the azies:' + str(azies) + ' here is azimuth:' + str(azi))
                print('Here is the inces:' + str(inces) + ' here is incident:' + str(arrival.incident_angle))
            print(arrival.name)
            mess.append({'SNR': signal/noise, 'azi': azi, 'Phase': arrival.name,  'Incident': arrival.incident_angle, 
            'IncidentE': inces, 'AziE': azies  , 'Year': eve.origins[0].time.year, 'Jday': eve.origins[0].time.julday,
            'PCA': angle, 'Lambda': + scale, 'Ray_Param': arrival.ray_param_sec_degree})
        
    return mess
    
    
    
    
    

sta = 'ANMO'
net ='IU'
loc ='10'
stime = UTCDateTime('2010-001T00:00:00.0')
etime = UTCDateTime('2019-001T00:00:00.0')

mess = getangles(net, sta, loc, stime, etime) 

f = open(net +'_' + sta + '_' + loc + 'Results','w')
f.write('Year, doy, phase, SNR, Aiz, AziE, Incident, IncidentE, PCA, Lambda, Ray Param\n')

for ele in mess:
    f.write(str(ele['Year']) + ', ' + str(ele['Jday']) + ', ' + str(ele['Phase']) + ', ' + str(ele['SNR']) + 
            ', ' + str(ele['azi']) + ', ' + str(ele['AziE']) +  ', ' + str(ele['Incident']) + ', ' + str(ele['IncidentE']) +
            ', ' + str(ele['PCA']) + ', ' + str(ele['Lambda'])  + ', ' + str(ele['Ray_Param']) + '\n')

f.close()


          

