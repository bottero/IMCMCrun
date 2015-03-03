# -*- coding: utf-8 -*-
"""
Created on Fri Dec 8 10:21:27 2014

Script to plot and represent our data for realistic IMCMC simulation

@author: alexis dot bottero At gmail dot com
"""

import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

plt.close('all')
stats=np.loadtxt("../data/coordStats.txt")
shots=np.loadtxt("../data/coordShots.txt")
logP=np.loadtxt("../data/falseLogP4096.txt")
logS=np.loadtxt("../data/falseLogS4096.txt")
realVelP=np.loadtxt("../data/realVelP.txt")
realVelS=np.loadtxt("../data/realVelS.txt")
if (realVelP.shape != realVelS.shape):
    print "!! realVelP doesn't have the same size than realVelS !!"
z=realVelP[:,0];

falseLogP=np.loadtxt("../../Data/falseLogP4096.txt")
falseLogS=np.loadtxt("../../Data/falseLogS4096.txt")
z=np.linspace(-999,0,4096)
np.savetxt("firstGuessP.txt", np.vstack([z[::-1],falseLogP]).T, delimiter=" ")
np.savetxt("firstGuessS.txt", np.vstack([z[::-1],falseLogS]).T, delimiter=" ")

#z2=np.linspace(0,1000,256)
#homoVelP=np.zeros(256)+1.0
#homoVelS=np.zeros(256)+2.0
#np.savetxt("homoVelP256.txt", np.vstack([z2,homoVelP]).T, delimiter=" ")
#np.savetxt("homoVelS256.txt", np.vstack([z2,homoVelP]).T, delimiter=" ")

x=np.linspace(0,500,4096)
#ax = Axes3D(plt.gcf())
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.hold(True)
#ax.plot(np.zeros(4096)+100,logP[0]/16*np.ones(4096),z,zdir='z',c='g',linestyle='--')
ax.plot(np.zeros(4096)+100,np.zeros(4096)+100,z,zdir='z',c='b',linestyle='--',label='Receivers')
ax.plot(np.zeros(4096)+200,np.zeros(4096)+400,z,zdir='z',c='b',linestyle='--')
ax.scatter(stats[:,0],stats[:,1],-stats[:,2],zdir='z',s=20,c='b')
ax.scatter(shots[:,0],shots[:,1],-shots[:,2],zdir='z',s=20,c='r',marker='^')
ax.plot(x,np.zeros(4096)+275,np.zeros(4096)-900,zdir='z',c='r',linestyle='--',label='Sources')
#ax.plot(np.zeros(4096),logP[::-1]/16,z,c='g',zdir='z')
#ax.plot(logS[::-1]/16,np.zeros(4096),z,c='pink',zdir='z')
ax.set_xlim3d(0, 500)
ax.set_ylim3d(0,500)
ax.set_zlim3d(-1000,0)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.legend()

plt.show()

fig2 = plt.figure()
plt.hold(True)
plt.plot(logP,z,color=(0.5,0.5,0.95))
plt.plot(logS,z,color=(0.5,0.95,0.5))
plt.plot(realVelP[:,1],z,color=(0,0,0.5),linewidth=4,label='Real P wave velocity')
plt.plot(realVelS[:,1],z,color=(0,0.5,0),linewidth=4,label='Real S wave velocity')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel(r'Wave speed ($m.s^{-1}$)',fontsize='14')
plt.ylabel(r'Depth ($m$)',fontsize='14')
plt.legend()

# How I have built the real velocity files :
#realVelP=np.loadtxt("../data/velpobs.txt")
#realVelS=np.loadtxt("../data/velsobs.txt")
#z=np.linspace(-999.5,0,4096)
#np.savetxt("realVelP.txt", np.vstack([z[::-1],realVelP]).T, delimiter=" ")
#np.savetxt("realVelS.txt", np.vstack([z[::-1],realVelS]).T, delimiter=" ")