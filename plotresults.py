##!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
# Import/apply font
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)
from scipy.optimize import fmin
import sys
alpha= 4.1
beta = 2.7

net = 'IU'
sta = 'ANMO'
loc = '10'
fname = net + '_' + sta + '_' + loc + 'Results2'

pp = []
psIE = []
pg =[]
ps = []
ssIE = []
sg =[]
ins =[]
inp =[]
pvec = []
snrp =[]
snrs =[]
with open(fname) as f:
    next(f)
    for line in f:
        line = line.replace(' ','')
        line = line.split(',')

        if line[2] == 'P':
            # This is the ray parameter in degrees/s
            pp.append(float(line[10]))
            # Here is theta bar
            psIE.append(float(line[8]))
            # Here is lambda
            pg.append(float(line[9]))
            
            # Here is the calculated incidenct angle
            inp.append(float(line[6]))
            snrp.append(float(line[3]))
        else:
            
            # This is the array parameter in degrees/s
            ps.append(float(line[10]))
            # Here is phi bar
            ssIE.append(float(line[8]))
            # Here is lambda
            sg.append(float(line[9]))
            # Here is the calculated incident angle
            ins.append(float(line[6]))
            snrs.append(float(line[3]))

print(len(pp))
print(len(ps))
snrs = np.asarray(snrs)
snrp = np.asarray(snrp)
pp = np.asarray(pp)
pg = np.asarray(pg)
ps = np.asarray(ps)
sg = np.asarray(sg)
# Now we want to invert for a best beta and alpha    
# Here we make functions 6 and 12

# Here is 6
def theta(b, pv):

    return 2.*np.arcsin(b*pv*np.pi/180.)

# Here is 12
def phi(b, a, pv):
    inside = (2.*(pv*np.pi/180.)*b**2)*np.sqrt(1.-((pv*np.pi/180.)*a)**2)/(a*1.-2.*(b*(pv*np.pi/180.))**2)
    return np.arctan(inside)


# We want to find b that minimizes theta and theta bar
# This is 14 with only the P-wave component
#def cost_funP(b):
    #res = np.sum(pg*(psIE -theta(b,pp))**2)/np.sum(pg)
    #print(res)
    #return res


# The function cost_funP is unstable so we need to use phi as well
def cost_fun(banda, pg, sg):
    # Here we do a joint inversion for b and a
    b = banda[0]
    a = banda[1]
    
    pterm = pg*(psIE - theta(b, pp))**2  
    sterm = sg*(ssIE - phi(b,a,ps))**2
    # Need to remove bad points coming from p
    ptermF = pterm[~np.isnan(pterm)] 
    stermF = sterm[~np.isnan(pterm)]
    pgF = pg[~np.isnan(pterm)]
    sgF = sg[~np.isnan(pterm)]
    # Need to remove bad points coming from s
    ptermF = ptermF[~np.isnan(stermF)]
    pgF = pgF[~np.isnan(stermF)]
    sgF = sgF[~np.isnan(stermF)]
    stermF = stermF[~np.isnan(stermF)]
    res = np.sum(ptermF + stermF)/np.sum(pgF+sgF)
    return res

def cost_fun_good(banda):
    return cost_fun(banda,pg, sg)

print(cost_fun_good([beta,alpha]))


for a in np.arange(0.001, 5., 0.01):
    for b in np.arange(.001,7., 0.01):
        if b >= np.sqrt(3.)*a/2.:
            continue
        else:
            if 'newsol' in vars():
                if (cost_fun_good([b,a]) < cost_fun_good(newsol)):
                    newsol = [b,a]
            else:
                newsol = [b,a]




print(newsol)




fig = plt.figure(1, figsize=(12,12))
plt.subplot(2,1,1)
plt.title(net + ' ' + sta + ' ' + loc + ' P:' + str(round(newsol[1],2)) + ' km/s S:' + str(round(newsol[0],2)) + 'km/s')
plt.scatter(pp, psIE, c=pg)
plt.plot(pp, inp, c='k')

plt.ylabel('$ \overline{\\theta} (^{\circ})$')
plt.xlabel('$p(s/  ^{\circ})$')
plt.colorbar()
#plt.xlim((min(pp), max(pp)))
plt.subplot(2,1,2)
plt.scatter(ps, ssIE, c=sg)
plt.plot(ps, ins, c='k')
plt.ylabel('$ \overline{\phi} (^{\circ})$')
plt.xlabel('$p (s/ ^{\circ})$')
#plt.xlim((min(ps), max(ps)))
plt.colorbar()
plt.savefig(net + '_' + sta + '_' + loc + '2.png', dpi=400)
plt.show()

