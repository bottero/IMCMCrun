# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 13:01:38 2014

Downsample (or upsample) a curve defined as :
  |       v[0]       |       v[1]       |        ...        |      v[nz-1]      |
z[0]               z[1]               z[2]                z[nz-1]             z[nz]
           |                  |        ...        |                   |
           v                  v                   v                   v
         zp[0]              zp[1]              zp[nz-2]              zp[nz-1]

To obtain :
  |         vp[0]         |         vp[0]        |  ... |     vp[nzFilt-1]      |
zFilt[0]                zFilt[1]               zFilt[2] zFilt[nzFilt-1]      zFilt[nzFilt]
              |                       |     ...      |              |
              v                       v              v              v
           zFiltp[0]                 zFiltp[1]  zFiltp[nzFilt-2] zFiltp[nzFilt-1]

If nzFilt < nz => downsample 
If nzFilt > nz => upsample
 
vp in intrapolated by a linear approximation

@author: alexis dot bottero at gmail dot com
"""

import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import matplotlib.pyplot as plt

def find_interval(array, value):
    """ Returns the index idxInf and idxSup of array verifying : 
        array[idxInf] < value < array[idxSup] """
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    idxInf=-1
    idxSup=-1
    if value < array.min() or value > array.max():
        idxInf=-1
        idxSup=-1
    elif array[idx] >= value and idx != 0:
        idxInf=idx-1
        idxSup=idx
    else:
        idxInf=idx
        idxSup=idx+1

    return idxInf,idxSup

def reSample(zp,v,nzFilt):
    """ Downsample (or upsample) a curve defined as :
      |       v[0]       |       v[1]       |        ...        |      v[nz-1]      |
    z[0]               z[1]               z[2]                z[nz-1]             z[nz]
               |                  |        ...        |                   |
               v                  v                   v                   v
             zp[0]              zp[1]              zp[nz-2]              zp[nz-1]
    
    To obtain :
      |         vp[0]         |         vp[0]        |  ... |     vp[nzFilt-1]      |
    zFilt[0]                zFilt[1]               zFilt[2] zFilt[nzFilt-1]      zFilt[nzFilt]
                  |                       |     ...      |              |
                  v                       v              v              v
               zFiltp[0]                 zFiltp[1]  zFiltp[nzFilt-2] zFiltp[nzFilt-1]
    
    If nzFilt < nz => downsample 
    If nzFilt > nz => upsample
     
    vp in intrapolated by a linear approximation
    
    @author: alexis dot bottero at gmail dot com
    """
    zFilt=np.linspace(zp.min(),zp.max(),nzFilt)
    zFiltp=np.zeros(nzFilt-1)
    for i in np.arange(1,len(zFilt)):
        zFiltp[i-1]=zFilt[i-1]+(zFilt[i]-zFilt[i-1])/2.0
    vFilt=np.zeros(nzFilt-1)
    for i in np.arange(0,len(vFilt)):
        idxInf,idxSup=find_interval(zp,zFiltp[i])
        vFilt[i]=(v[idxSup]-v[idxInf])/(zp[idxSup]-zp[idxInf])*(zFiltp[i]-zp[idxInf])+v[idxInf]
    return zFiltp,vFilt

#zmin=0.0
#zmax=10.0
#nz=150    # Number of points describing the curve that we want to downsample
#nzFilt=50 # Number of border points in the downsampled curve (nzFilt < nz)
#z=np.linspace(zmin,zmax,nz)
#zp=np.zeros(nz-1)
#for i in np.arange(1,len(z)):
#    zp[i-1]=z[i-1]+(z[i]-z[i-1])/2.0
#v=3.0*np.cos(zp/2.5)+4000.0;
#
#zFiltp,vFilt = reSample(zp,v,nzFilt)
#
#plt.hold(True)
#plt.plot(zp,v,'x',color='blue')
#plt.plot(zFiltp,vFilt,'x',color='red')

logP=np.loadtxt("VP-botteroRS4.txt")
logS=np.loadtxt("VS-botteroRS4.txt")
zmin=0.0
zmax=1219.051181
nz=len(logP)
z=np.linspace(zmin,zmax,nz)
zp=np.zeros(nz)
dz=(z[1]-z[0])/2.0
for i in np.arange(0,len(z)):
    zp[i]=z[i]+dz

nzFilt=2049
zFiltp,logPfilt = reSample(zp,logP,nzFilt)
zFiltp,logSfilt = reSample(zp,logS,nzFilt)

plt.figure()
plt.hold(True)
plt.plot(zp,logP,color='blue')
plt.plot(zFiltp,logPfilt,color='red')
plt.figure()
plt.hold(True)
plt.plot(zp,logS,color='blue')
plt.plot(zFiltp,logSfilt,color='red')

np.savetxt("logPreal.txt", np.vstack([zFiltp,logPfilt]).T, delimiter=" ")
np.savetxt("logSreal.txt", np.vstack([zFiltp,logSfilt]).T, delimiter=" ")



