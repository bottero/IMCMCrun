#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:08:14 2014

@author: bottero

Script to watch the results of a run.

TODO : _load prior profiles as well
       _

"""
### --- MODULES AND PACKAGES --- ###
import os, sys
import argparse # To deal with arguments :
# https://docs.python.org/2/library/argparse.html
import numpy as np # NumPy (multidimensional arrays, linear algebra, ...)
import glob # Unix style pathname pattern expansion
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt # Matplotlib's pyplot: MATLAB-like syntax

def representsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def representsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

### --- Parser arguments --- ###
parser = argparse.ArgumentParser(description='Watch the results of a run')
parser.add_argument("pathToDir",
                    help="Path to result directory (ex : OUTPUT_FILES/XXX)")
parser.add_argument("-v","--verbose", help="Increase output verbosity",
                    action="store_true")
parser.add_argument("-a","--all", help="Plot everything",
                    action="store_true", default=False)
parser.add_argument("--recalculate_t0", help="Recalculate t0",
                    action="store_true")
parser.add_argument("-d","--data", help="Plot the first guess curve",
                    action="store_true")
parser.add_argument("-t","--times", help="Plot the calculated times",
                    action="store_true")
parser.add_argument("-g","--geometry", help="Plot the geometry used (sources and receivers)",
                    action="store_true")
parser.add_argument("--homogeneousVelocityModel", help="Consider that the velocity model is constant",
                    action="store_true")
args = parser.parse_args()

### --- Test arguments --- ###

if not os.path.isdir(args.pathToDir): # If the path does not exist
    print "Directory ",args.pathToDir," not found."
    parser.print_help()
    sys.exit(0)
if args.pathToDir[-1:] != '/': # add a / at the end of the path if necessary
    args.pathToDir+='/'
code = glob.glob1(args.pathToDir,"calculatedTimes.*")[0].split('.')[1]

if not representsInt(code):
    print "Directory ",args.pathToDir," does not seem to be a correct IMCMC \
directory... (no stats0 found)"
    sys.exit(0)

if not (args.all or args.data or args.geometry or args.times):
    print "Nothing has to be done!"
    sys.exit(0)

### --- Load files --- ###

# Extract informations from config.XXX.dat :
# TODO : read the cfg file

T = np.zeros(0)
with open(args.pathToDir+'config.cfg') as configFile:
    for line in configFile:
        if line.split('#')[0].strip(): # If the line is not a comment
            if 'SWAVES' in line:
                if line.split(" = ")[1].split("#")[0].strip() == "1":
                    swaves=True
                else:
                    swaves=False
            if 'ANALYTICAL_RUN' in line:
                if line.split(" = ")[1].split("#")[0].strip() == "1":
                    analytical=True
                else:
                    analytical=False
            if 'NAME_OF_FIRST_GUESS_P_FILE' in line:
                nameOfFirstGuessP = line.split(" = ")[1].split("#")[0].strip()
            if 'NAME_OF_FIRST_GUESS_S_FILE' in line:
                nameOfFirstGuessS = line.split(" = ")[1].split("#")[0].strip()
            if 'NAME_OF_STATIONS_FILE' in line:
                nameOfStationsFile = line.split(" = ")[1].split("#")[0].strip()
            if 'NAME_OF_SHOTS_FILE' in line:
                nameOfShotsFile = line.split(" = ")[1].split("#")[0].strip()
            if 'NAME_OF_TIMES_FILE =' in line:
                nameOfTimesFile = line.split(" = ")[1].split("#")[0].strip()
            if 'SIGMAP =' in line:
                sigmaP = line.split(" = ")[1].split("#")[0].strip()
                if representsFloat(sigmaP): # Verify that the string extracted is a int
                    sigmaP = float(sigmaP)
            if 'SIGMAS =' in line:
                sigmaS = line.split(" = ")[1].split("#")[0].strip()
                if representsFloat(sigmaS): # Verify that the string extracted is a int
                    sigmaS = float(sigmaS)
            if 'NX_DEFAULT' in line:
                nx = line.split(" = ")[1].split("#")[0].strip()
                if representsInt(nx): # Verify that the string extracted is a int
                    nx = int(nx)
            if 'NY_DEFAULT' in line:
                ny = line.split(" = ")[1].split("#")[0].strip()
                if representsInt(ny): # Verify that the string extracted is a int
                    ny = int(ny)
            if 'COORD_TOL =' in line:
                coord_tol = line.split(" = ")[1].split("#")[0].strip()
                if representsFloat(coord_tol): # Verify that the string extracted is a int
                    coord_tol = float(coord_tol)
if args.verbose:
    print "Watching the results of run : ",code,"..."
    if analytical:
        print "This is an analytical run"
    if swaves:
        print "S waves arrival times were calculated"
    print
    print "Loading files ..."

### --- Load files --- ###
# declare empty list to store files :
timesData=[0]
coordShots=np.loadtxt(args.pathToDir+nameOfShotsFile)
coordStats=np.loadtxt(args.pathToDir+nameOfStationsFile)
firstGuessP=np.loadtxt(args.pathToDir+nameOfFirstGuessP)
if swaves:
    firstGuessS=np.loadtxt(args.pathToDir+nameOfFirstGuessS)
calculatedTimes=np.loadtxt(args.pathToDir+"calculatedTimes."+code+".dat")
if not analytical:
    timesData=np.loadtxt(args.pathToDir+nameOfTimesFile)

if args.verbose:
    print "Loading done !"

xmin=-1e99
ymin=-1e99
zmin=-1e99
xmax=1e99
ymax=1e99
zmax=1e99
### --- Analyses --- ###
if np.size(coordShots) > 3:
    xmin = coordShots[:,0].min()
    ymin = coordShots[:,1].min()
    zmin = coordShots[:,2].min()
    xmax = coordShots[:,0].max()
    ymax = coordShots[:,1].max()
    zmax = coordShots[:,2].max()
else:
    xmin = coordShots[0]
    ymin = coordShots[1]
    zmin = coordShots[2]
    xmax = coordShots[0]
    ymax = coordShots[1]
    zmax = coordShots[2]

if np.size(coordStats) > 3:
    xmin = min(xmin,coordStats[:,0].min())
    ymin = min(ymin,coordStats[:,1].min())
    zmin = min(zmin,coordStats[:,2].min())
    xmax = max(xmax,coordStats[:,0].max())
    ymax = max(ymax,coordStats[:,1].max())
    zmax = max(zmax,coordStats[:,2].max())
else:
    xmin = min(xmin,coordStats[0])
    ymin = min(ymin,coordStats[1])
    zmin = min(zmin,coordStats[2])
    xmax = max(xmax,coordStats[0])
    ymax = max(ymax,coordStats[1])
    zmax = max(zmax,coordStats[2])

xmin2 = xmin - (xmax-xmin)*coord_tol;
ymin2 = ymin - (ymax-ymin)*coord_tol;
xmax2 = xmax + (xmax-xmin)*coord_tol;
ymax2 = ymax + (ymax-ymin)*coord_tol;
xmin = xmin2
ymin = ymin2
xmax = xmax2
ymax = ymax2

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);

z=firstGuessP[:,0]
nz = len(z) + 1
dz=z[1]-z[0]
zminProfile=z[0]-dz/2
zmaxProfile=z[-1]+dz/2
if zmax > zmaxProfile:
    print "PROBLEM zmax > zmaxProfile (zmax: ",zmax," zmaxProfile: ",zmaxProfile,")"
zmax = zmaxProfile
if zmin < zminProfile:
    print "PROBLEM zmin < zminProfile (zmin: ",zmin," zminProfile: ",zminProfile,")"
zmin = zminProfile

nStats=coordStats.size/3
nShots=coordShots.size/3

if args.verbose:
    print "nx : ",nx," ny : ",ny," nz : ",nz
    print "dx : ",dx," dy : ",dy," dz : ",dz

epsilonX = 0
epsilonY = 0
epsilonZ = 0

if nShots == 1:
    epsilonX = coordShots[0] - xmin - np.floor((coordShots[0]-xmin)/dx)*dx;
    epsilonY = coordShots[1] - ymin - np.floor((coordShots[1]-ymin)/dy)*dy;
    epsilonZ = coordShots[2] - zmin - np.floor((coordShots[2]-zmin)/dz)*dz;
    xmin = xmin + epsilonX;
    ymin = ymin + epsilonY;
    zmin = zmin + epsilonZ;
    if args.verbose:
        print "xmin : ",xmin," ymin : ",ymin," zmin : ",zmin
else:
    epsilonX = np.zeros(nShots)
    epsilonY = np.zeros(nShots)
    epsilonZ = np.zeros(nShots)
    xminShots = np.zeros(nShots)
    yminShots = np.zeros(nShots)
    zminShots = np.zeros(nShots)
    for ishot in np.arange(nShots):
        epsilonX[ishot] = coordShots[ishot,0] - xmin - np.floor((coordShots[ishot,0]-xmin)/dx)*dx;
        epsilonY[ishot] = coordShots[ishot,1] - ymin - np.floor((coordShots[ishot,1]-ymin)/dy)*dy;
        epsilonZ[ishot] = coordShots[ishot,2] - zmin - np.floor((coordShots[ishot,2]-zmin)/dz)*dz;
        xminShots[ishot] = xmin + epsilonX[ishot];
        yminShots[ishot] = ymin + epsilonY[ishot];
        zminShots[ishot] = zmin + epsilonZ[ishot];
        if args.verbose:
            print "Shot number ",ishot
            print "xmin : ",xminShots[ishot]," ymin : ",yminShots[ishot]," zmin : ",zminShots[ishot]

plt.close('all')

if args.all:
    args.geometry=True
    args.data=True
    args.times=True

if args.geometry:
    if swaves:
        timesP=calculatedTimes[:,0]
    else:
        timesP=calculatedTimes
    if args.homogeneousVelocityModel:
        speed = firstGuessP[0,1]
        realTimes=np.sqrt((coordStats - coordShots)**2)/speed
        realTimes=realTimes[realTimes!=0]
        residuals=timesP-realTimes
    fig = plt.figure()
    ax = fig.gca(projection='3d') #Axes3D(fig)
    ax.hold(True)
    for stat in np.arange(len(coordStats[:,0])):
        if args.homogeneousVelocityModel:
            if (abs(residuals[stat]) > 0.01/speed):
                ax.scatter(coordStats[stat,0],coordStats[stat,1],coordStats[stat,2],zdir='z',s=20,c='r')
        else:
            ax.scatter(coordStats[stat,0],coordStats[stat,1],coordStats[stat,2],zdir='z',s=20,c='g')
    if (coordShots.size>3):
        ax.scatter(coordShots[:,0],coordShots[:,1],coordShots[:,2],zdir='z',s=20,c='b',marker='^')
    else:
        ax.scatter(coordShots[0],coordShots[1],coordShots[2],zdir='z',s=200,c='b',marker='^')
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    ax.set_zlim3d(zmin,zmax)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Geometry')
    ax.invert_zaxis()

if args.data:
    fig2 = plt.figure()
    plt.hold(True)
    plt.plot(firstGuessP[:,1],z,color=(0.5,0.5,0.95))
    if (swaves):
        plt.plot(firstGuessS[:,1],z,color=(0.5,0.95,0.5))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(r'Wave speed ($m.s^{-1}$)',fontsize='14')
    plt.ylabel(r'Depth ($m$)',fontsize='14')
    if swaves:
        plt.title(r'First guess velocity profiles',fontsize='14')
    else:
        plt.title(r'First guess velocity profile',fontsize='14')
    plt.ylim(ymax=z.max())
    plt.gca().invert_yaxis()

if args.times:
    if not analytical:
        diffDataBestModel=calculatedTimes-timesData
        if args.recalculate_t0:
            for i in np.arange(nShots):
                diffPshoti=diffDataBestModel[i*nStats:(i+1)*nStats,0][timesData[i*nStats:(i+1)*nStats,0]>0]
                t0ShotsPi=diffPshoti.mean()
                if args.verbose:
                    print "t0P[",i,"] = ",t0ShotsPi
                diffDataBestModel[i*nStats:(i+1)*nStats,0]=diffDataBestModel[i*nStats:(i+1)*nStats,0]-t0ShotsPi
                if swaves:
                    diffSshoti=diffDataBestModel[i*nStats:(i+1)*nStats,1][timesData[i*nStats:(i+1)*nStats,1]>0]
                    t0ShotsSi=diffSshoti.mean()
                    if args.verbose:
                        print "t0S[",i,"] = ",t0ShotsSi
                    diffDataBestModel[i*nStats:(i+1)*nStats,1]=diffDataBestModel[i*nStats:(i+1)*nStats,1]-t0ShotsSi
        diffP2=diffDataBestModel[:,0]
        diffP2[timesData[:,0]<0] = -0.0005
        diffS2=diffDataBestModel[:,1]
        diffS2[timesData[:,1]<0] = -0.0005
        np.savetxt(args.pathToDir+"residuals_model.txt",np.dstack((diffP2,diffS2))[0])
        diffP=diffDataBestModel[:,0][timesData[:,0]>0]
        diffS=diffDataBestModel[:,1][timesData[:,1]>0]
        fig = plt.figure()
        plt.hold(True)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+sigmaP,'b--',linewidth=2)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-sigmaP,'b--',linewidth=2)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+2*sigmaP,'--',color=(0.3,0.3,1),linewidth=1.5)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-2*sigmaP,'--',color=(0.3,0.3,1),linewidth=1.5)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+3*sigmaP,'--',color=(0.5,0.5,1))
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-3*sigmaP,'--',color=(0.5,0.5,1))
        plt.plot(np.arange(len(diffP)),diffP,'g+')
        plt.ylim([-5*sigmaP,5*sigmaP])
        plt.xlim([0,len(diffP)-1])
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.text(0.6666*len(diffP), 1.1*sigmaP, r'$1\sigma$',fontsize='30',color='b')
        plt.text(0.7708*len(diffP), 2.1*sigmaP, r'$2\sigma$',fontsize='30',color=(0.3,0.3,1))
        plt.text(0.875*len(diffP), 3.1*sigmaP, r'$3\sigma$',fontsize='30',color=(0.5,0.5,1))
        plt.xlabel(r'Receiver number',fontsize='14')
        plt.ylabel(r'P waves arrival times residuals ($s$)',fontsize='14')
        plt.title(r'Best model residuals',fontsize='14')
        if swaves:
            fig2 = plt.figure()
            plt.hold(True)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+sigmaS,'b--',linewidth=2)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-sigmaS,'b--',linewidth=2)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+2*sigmaS,'--',color=(0.3,0.3,1),linewidth=1.5)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-2*sigmaS,'--',color=(0.3,0.3,1),linewidth=1.5)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+3*sigmaS,'--',color=(0.5,0.5,1))
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-3*sigmaS,'--',color=(0.5,0.5,1))
            #plt.plot(np.arange(len(diffS)),diffDataBestModel[:,1],'g+')
            plt.plot(np.arange(len(diffS)),diffS,'g+')
            plt.ylim([-5*sigmaS,5*sigmaS])
            plt.xlim([0,len(diffS)-1])
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.text(0.6666*len(diffS), 1.1*sigmaS, r'$1\sigma$',fontsize='30',color='b')
            plt.text(0.7708*len(diffS), 2.1*sigmaS, r'$2\sigma$',fontsize='30',color=(0.3,0.3,1))
            plt.text(0.875*len(diffS), 3.1*sigmaS, r'$3\sigma$',fontsize='30',color=(0.5,0.5,1))
            plt.xlabel(r'Receiver number',fontsize='14')
            plt.ylabel(r'S waves arrival times residuals ($s$)',fontsize='14')
            plt.title(r'Best model residuals',fontsize='14')
            E=sum((diffP/sigmaP)**2/2)+sum((diffS/sigmaS)**2/2) #+ep
        else:
            E=sum((diffP/sigmaP)**2/2) #+ep
        if args.verbose:
            print "Energy of best model (without prior) :",E
    else:
        fig = plt.figure()
        if swaves:
            timesP=calculatedTimes[:,0]
        else:
            timesP=calculatedTimes
        plt.hold(True)
        plt.plot(np.arange(len(timesP)),timesP,'g+')
        plt.xlim([0,len(timesP)-1])
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel(r'Receiver number',fontsize='14')
        plt.ylabel(r'P waves arrival times ($s$)',fontsize='14')
        if swaves:
            timesS=calculatedTimes[:,1]

plt.show()