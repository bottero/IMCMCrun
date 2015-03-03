# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:52:12 2014
Script to create synthetic profiles

@author: bottero
"""
import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import matplotlib.pyplot as plt
npoints=512
zmin=100
zmax=500
profileP=np.zeros(npoints)+1
profileS=np.zeros(npoints)+0.5
guessProfileP=np.zeros(npoints)+1
guessProfileS=np.zeros(npoints)+0.5
z=np.linspace(zmin,zmax,npoints)

#for i in np.arange(npoints):
#    if i< npoints/2:
#        profileP[i]=4000
#        profileS[i]=2500
#    else:
#        profileP[i]=5000
#        profileS[i]=3500
#
#for i in np.arange(npoints):
#    if i< npoints/2.5:
#        guessProfileP[i]=4500
#        guessProfileS[i]=2000
#    else:
#        guessProfileP[i]=6000
#        guessProfileS[i]=3000

np.savetxt("realPhomo512.txt", np.vstack([z,profileP]).T, delimiter=" ")
np.savetxt("realShomo512.txt", np.vstack([z,profileS]).T, delimiter=" ")
np.savetxt("guessPhomo512.txt", np.vstack([z,guessProfileP]).T, delimiter=" ")
np.savetxt("guessShomo512.txt", np.vstack([z,guessProfileS]).T, delimiter=" ")
