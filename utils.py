#!/usr/bin/env python

def get_parameters(phase):
    paramdic = {}
    if phase == 'P':
        
        paramdic['station_radius'] = 2.0
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 1./20.
        paramdic['fmax'] =1./5.
    elif phase == 'Rayleigh':
        paramdic['station_radius'] = 1.5
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 0.1
        paramdic['fmax'] =1.
    elif  phase == 'S':
        paramdic['station_radius'] = 1.5
        paramdic['min_radius'] = 30.
        paramdic['max_radius'] = 80.
        paramdic['min_mag'] = 6.
        paramdic['max_mag'] = 8.5
        paramdic['length'] = 50
        paramdic['fmin'] = 0.1
        paramdic['fmax'] =1.
    else:
        pass
    
    return paramdic
