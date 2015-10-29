# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 16:15:45 2015

@author: bottero
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
NPU=10
ztop=-100
zbottom=-500
N=512
N2=10
zp=np.linspace(ztop,zbottom,N)
vpPrior=np.cos(zp/50)
#vpPrior = np.zeros(512)+8000
#for i in np.arange(295):
#    if i < 56:
#        vpPrior[i] = 6000
#    elif i >= 56 and i < 295:
#        vpPrior[i] = 9000

optimumZlayers=np.zeros(NPU-1)
zLayers=np.zeros(NPU-1)
indices=np.zeros(NPU-1)
vp=np.zeros(NPU)

dzNPU=(zbottom-ztop)/NPU
print "dz:",dzNPU
optimumZlayers[0]=ztop+dzNPU
print "optimumZlayers:",optimumZlayers
indices[0]=int(np.searchsorted(-zp,-optimumZlayers[0]))
print "indices:",indices
zLayers[0]=zp[indices[0]]
print "zLayers:",zLayers
vp[0]=np.mean(vpPrior[0:indices[0]])
print "vp:",vp

for i in 1+np.arange(NPU-2):
    optimumZlayers[i]=optimumZlayers[i-1]+dzNPU
    indices[i]=np.searchsorted(-zp,-optimumZlayers[i])
    zLayers[i]=zp[indices[i]]
    vp[i]=np.mean(vpPrior[indices[i-1]:indices[i]+1])

vp[NPU-1]=np.mean(vpPrior[indices[NPU-2]:-1])

print "optimumZlayers:",optimumZlayers
print "indices:",indices
print "zLayers:",zLayers
print "vp:",vp

zStepCenter=np.zeros(NPU)
zStepCenter[0:-1]=zLayers-dzNPU/2
zStepCenter[-1]=zLayers[-1]+dzNPU/2


vpFilt = np.zeros(N)
for i in np.arange(N):
    if i < indices[0]:
        vpFilt[i] = vp[0]
    elif i >= indices[-1]:
        vpFilt[i] = vp[-1]
    else:
        for j in np.arange(NPU-2):
            if i >= indices[j] and i < indices[j+1]:
                vpFilt[i] = vp[j+1]

print "zStepCenter:",zStepCenter

plt.figure()
plt.hold(True)
plt.plot(zp,vpPrior,'+-')
#for i in np.arange(NPU-2)+1:
#    rinf=(zp[0]-zLayers[i-1])/(zp[0]-zp[-1])
#    rsup=(zp[0]-zLayers[i])/(zp[0]-zp[-1])
#    plt.axhline(vp[i],xmin=rinf,xmax=rsup)
for zz in zLayers:
    plt.axvline(zz)
plt.plot(zStepCenter,vp,'o')
plt.plot(zp,vpFilt,'+-')