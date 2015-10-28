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

def pdense(x, y, sigma, M=1000):
    """ Plot probability density of y with known stddev sigma
    """
    assert len(x) == len(y) and len(x) == len(sigma)
    N = len(x)
    # TODO: better y ranging
    ymin, ymax = min(y - 2 * sigma), max(y + 2 * sigma)
    yy = np.linspace(ymin, ymax, M)
    a = [np.exp(-((Y - yy) / s) ** 2) / s for Y, s in zip(y, sigma)]
    A = np.array(a)
    A = A.reshape(N, M)
    plt.imshow(-A, cmap='gray', aspect='auto',
               origin='upper', extent=(ymin,ymax,max(x),min(x)))
   # plt.title('Density plot')

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
parser.add_argument("-d","--data", help="Plot the first guess curve (and real curves for analytical runs))",
                    action="store_true")
parser.add_argument("-g","--geometry", help="Plot the geometry used (sources and receivers)",
                    action="store_true")
parser.add_argument("-e","--energies", help="Plot the energies of the chains",
                    action="store_true")
parser.add_argument("-b","--best", help="Plot the some good models found and the residuals of the best",
                    action="store_true")
parser.add_argument("--dont_show_guess", help="Don't show first guess model on plots",
                    action="store_true", default=False)
parser.add_argument("--show_ranges", help="Show exploration ranges on plots",
                    action="store_true", default=False)
parser.add_argument("--show_averages", help="Show average models on plots",
                    action="store_true", default=False)
parser.add_argument("--no_density_plots", help="Represent uncertainties by a range",
                    action="store_true", default=False)
parser.add_argument("-r","--results",type=int, 
                    help="Plot the results from the inversion for given \
chain number",default=-1)
parser.add_argument("-t","--treshold",type=int, 
                    help="Iteration after which we show the model calculated",
                    default=0)
parser.add_argument("--recalculate_t0", help="Force recalculate t0 (even if it had not been chosen during the run)",
                    action="store_true")
parser.add_argument("-s","--swaps", help="Plot the swaps",
                    action="store_true")
parser.add_argument("--vpvs", help="Plot mean Vp/Vs ratio plus uncertainties",
                    action="store_true")
args = parser.parse_args()

### --- Test arguments --- ###

if not os.path.isdir(args.pathToDir): # If the path does not exist
    print "Directory ",args.pathToDir," not found."
    parser.print_help()
    sys.exit(0)
if args.pathToDir[-1:] != '/': # add a / at the end of the path if necessary
    args.pathToDir+='/'
code = glob.glob1(args.pathToDir,"stats0.*")[0].split('.')[1] 
#args.pathToDir.split("/")[-2] # just keep the name of the directory
if not representsInt(code):
    print "Directory ",args.pathToDir," does not seem to be a correct IMCMC \
directory... (no stats0 found)"
    sys.exit(0)

if not (args.all or args.data or args.geometry or args.energies or args.best or (args.results != -1) or args.swaps or args.vpvs):
    print "Nothing has to be done!"
    sys.exit(0)

### --- Load files --- ###

# Extract informations from config.XXX.dat :
# TODO : read the cfg file

T = np.zeros(0)
with open(args.pathToDir+'config.cfg') as configFile:
    for line in configFile:
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
        if 'RECALCULATE_T0' in line:
            if line.split(" = ")[1].split("#")[0].strip() == "1":
                recalculate_t0=True
            else:
                recalculate_t0=False
        if 'QP =' in line:
            qp = line.split(" = ")[1].split("#")[0].strip()
        if 'NAME_OF_FIRST_GUESS_P_FILE' in line:
            nameOfFirstGuessP = line.split(" = ")[1].split("#")[0].strip()
        if 'NAME_OF_FIRST_GUESS_S_FILE' in line:
            nameOfFirstGuessS = line.split(" = ")[1].split("#")[0].strip()
        if 'NAME_OF_REAL_PROFILE_FILE_P' in line:
            nameOfrealP = line.split(" = ")[1].split("#")[0].strip()
        if 'NAME_OF_REAL_PROFILE_FILE_S' in line:
            nameOfrealS = line.split(" = ")[1].split("#")[0].strip()
        if 'N_PRIOR_PROFILES' in line:
            nPriorProfiles = line.split(" = ")[1].split("#")[0].strip()
            if representsInt(nPriorProfiles): # Verify that the string extracted is a int
                nPriorProfiles = int(nPriorProfiles)
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

with open(args.pathToDir+'config.'+code+'.dat') as outConfigFile:
    for line in outConfigFile:
        if 'Energy of the prior : ' in line:
            ep=line.split(":")[1].strip("\n").strip()
            if representsFloat(ep): # Verify that the string extracted is a int
                ep = float(ep)
        if 'temperatures : ' in line:
            # Extract the number of temperatures from the line :
            nbt=(line.split(":")[-1]).split("\n")[0].strip() 
            if representsInt(nbt): # Verify that the string extracted is a int
                nbt = int(nbt)
                T=np.resize(T,nbt)
        if 'Temperature max : ' in line:
            # Extract temperature max from the line :
            tmax=(line.split(":")[-1]).split("\n")[0].strip() 
            if representsInt(tmax): # Verify that the string extracted is a int
                tmax = int(tmax)
        if 'xmin:' in line:
            xmin=line.split(":")[1].strip().split()[0]
            if representsFloat(xmin): # Verify that the string extracted is a int
                xmin = float(xmin)
            xmax=line.split(":")[2].strip().split()[0]
            if representsFloat(xmax): # Verify that the string extracted is a int
                xmax = float(xmax)
        if 'ymin:' in line:
            ymin=line.split(":")[1].strip().split()[0]
            if representsFloat(ymin): # Verify that the string extracted is a int
                ymin = float(ymin)
            ymax=line.split(":")[2].strip().split()[0]
            if representsFloat(ymax): # Verify that the string extracted is a int
                ymax = float(ymax)
        if 'zmin:' in line:
            zmin=line.split(":")[1].strip().split()[0]
            if representsFloat(zmin): # Verify that the string extracted is a int
                zmin = float(zmin)
            zmax=line.split(":")[2].strip().split()[0]
            if representsFloat(zmax): # Verify that the string extracted is a int
                zmax = float(zmax)
        if 'Temperature ladder : ' in line: # Extract temperature ladder from line
            T[nbt-1]=line.split("T["+str(nbt-1)+"] = ")[1].strip("\n").strip()
            for i in np.arange(nbt-1):
                T[i]=line.split("T["+str(i)+"] = ")[1].split("T["+str(i+1)+"] = ")[0].strip()
                
if args.verbose:
    print "Watching the results of run : ",code,"..."
    if analytical:
        print "This is an analytical run"
    if swaves:
        print "S waves arrival times were calculated"
    print "There are ",nbt," temperatures (tmax =",tmax,") :",
    for i in np.arange(nbt):
        print " T[",i,"] = ",T[i],
    print
    print "Loading files ..."
    
### --- Load files --- ###
# declare empty list to store files :
M=1000 # For density plot
averagesP=[0]*nbt # declare empty list to store files :
averagesS=[0]*nbt
chains=[0]*nbt
varPs=[0]*nbt
varSs=[0]*nbt
qInfPs=[0]*nbt
qSupPs=[0]*nbt
qInfSs=[0]*nbt
qSupSs=[0]*nbt
minP=[0]*nbt
minS=[0]*nbt
maxP=[0]*nbt
maxS=[0]*nbt
priorP=[0]*nPriorProfiles
priorS=[0]*nPriorProfiles
timesData=[0]

# Loop on temperature chains
for i in np.arange(nbt):
    averagesP[i]=np.loadtxt(args.pathToDir+"averageP"+str(i)+"."+code+".dat")
    chains[i]=np.loadtxt(args.pathToDir+"chain"+str(i)+"."+code+".dat")
    varPs[i]=np.loadtxt(args.pathToDir+"varP"+str(i)+"."+code+".dat")
    qSupPs[i]=np.loadtxt(args.pathToDir+"qSupP"+str(i)+"."+code+".dat")
    qInfPs[i]=np.loadtxt(args.pathToDir+"qInfP"+str(i)+"."+code+".dat")
    minP[i]=np.loadtxt(args.pathToDir+"minP."+str(i)+"."+code+".dat")
    maxP[i]=np.loadtxt(args.pathToDir+"maxP."+str(i)+"."+code+".dat")
    if swaves:
        averagesS[i]=np.loadtxt(args.pathToDir+"averageS"+str(i)+"."+code+".dat")
        varSs[i]=np.loadtxt(args.pathToDir+"varS"+str(i)+"."+code+".dat")
        qSupSs[i]=np.loadtxt(args.pathToDir+"qSupS"+str(i)+"."+code+".dat")
        qInfSs[i]=np.loadtxt(args.pathToDir+"qInfS"+str(i)+"."+code+".dat")
        minS[i]=np.loadtxt(args.pathToDir+"minS."+str(i)+"."+code+".dat")
        maxS[i]=np.loadtxt(args.pathToDir+"maxS."+str(i)+"."+code+".dat")
globalMaxP=np.loadtxt(args.pathToDir+"maxP."+code+".dat") #TODO don't forget them
globalMinP=np.loadtxt(args.pathToDir+"minP."+code+".dat")
if swaves:
    globalMaxS=np.loadtxt(args.pathToDir+"maxS."+code+".dat")
    globalMinS=np.loadtxt(args.pathToDir+"minS."+code+".dat")
coordShots=np.loadtxt(args.pathToDir+nameOfShotsFile)
coordStats=np.loadtxt(args.pathToDir+nameOfStationsFile)
firstGuessP=np.loadtxt(args.pathToDir+nameOfFirstGuessP)
firstGuessS=np.loadtxt(args.pathToDir+nameOfFirstGuessS)
if analytical:
    realP=np.loadtxt(args.pathToDir+nameOfrealP)
    realS=np.loadtxt(args.pathToDir+nameOfrealS)
    timesData=np.loadtxt(args.pathToDir+"calculatedTimes."+code+".dat")
else:
    timesData=np.loadtxt(args.pathToDir+nameOfTimesFile)
if os.path.isfile(args.pathToDir+"bestModelTimes."+code+".dat"):
    bestModelCalculated=True
    bestModelTimes=np.loadtxt(args.pathToDir+"bestModelTimes."+code+".dat")
else:
    bestModelCalculated=False
ll=np.loadtxt(args.pathToDir+"ll."+code+".dat")

if nPriorProfiles > 0:
    for i in np.arange(nPriorProfiles):
        priorP[i]=np.loadtxt(args.pathToDir+"priorProfiles"+code+"/priorProfileP."+code+"."+str(i)+".dat")
        priorS[i]=np.loadtxt(args.pathToDir+"priorProfiles"+code+"/priorProfileS."+code+"."+str(i)+".dat")

nit=len(chains[0])

if args.verbose:
    print "Loading done !"
    print
    print "During this simulation the",nbt,"chains have run during",nit,"steps"

if nit < 50:
    print "Take care!! Below 50 iterations min and max profiles don't make sense!"


### --- Analyses --- ###

z=firstGuessP[:,0]
zFilt=varPs[0][:,0]
nStats=coordStats.size/3
nShots=coordShots.size/3
plt.close('all')

if args.all:
    args.geometry=True
    args.energies=True
    args.best=True
    args.data=True
    if swaves:
        args.vpvs=True
   # args.swaps=True
    
if args.energies:
    plt.hold(True)
    ii=0
    nBest=len(glob.glob1(args.pathToDir,"bestPprofile*"))
    E=[0]*nBest    
    iterationBest=[0]*nBest
    chain=[0]*nBest
    for bestModel in glob.glob1(args.pathToDir,"bestSprofile*"):
        ii=ii+1
        iterationBest[ii-1]=int(bestModel.split("idx")[1].split(".")[0])
        chain[ii-1]=int(bestModel.split("chain")[1].split(".")[0])
        E[ii-1]=float(bestModel.split("E")[1].split(code)[0].strip(".")) 
        # TODO : this does not work if the energy value contains the code of the run : ex code=745 energy=745.23
    if args.verbose:
        print "Models kept after iteration : "+str(args.treshold)+" will be shown"
    from operator import itemgetter
    idxBest=min(enumerate(E), key=itemgetter(1))[0] # index of best model
    itBestE=iterationBest[idxBest]
    chainBest=chain[idxBest]
    chain=[chain[i] for i in np.arange(nBest) if iterationBest[i] > args.treshold]
    iterationBest=[i for i in iterationBest if i>args.treshold]
    iteration=np.arange(nit)
    for i in np.arange(nbt):
        plt.semilogy(iteration,chains[i][:,-1]*T[i])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    flagForLabel=True
    for j in np.arange(len(iterationBest)):
        if flagForLabel:
            plt.semilogy(iteration[iterationBest[j]], T[chain[j]]*chains[chain[j]][:,-1][iterationBest[j]], 'bD', label="Best models saved")
            flagForLabel=False
        else:
            plt.semilogy(iteration[iterationBest[j]], T[chain[j]]*chains[chain[j]][:,-1][iterationBest[j]], 'bD')
    plt.semilogy(itBestE, T[chainBest]*chains[chainBest][:,-1][itBestE], 'rD', label="Best model")
    if recalculate_t0 is True:
        if swaves:
            plt.semilogy(iteration,np.zeros(nit)+nStats*nShots+ep,'b--',linewidth=2,label=r'Behind that line every model can be acceptable ($1\sigma$ misfit for each measurement)')
        else:
            plt.semilogy(iteration,np.zeros(nit)+nStats*nShots/2.0+ep,'b--',linewidth=2)
    plt.legend()
    plt.semilogy(iteration,np.zeros(nit)+ep)
    plt.xlim(xmax=iteration.max())
    plt.rc('font', family='serif')
    plt.xlabel('Iteration number',fontsize='14')
    plt.ylabel('Energy',fontsize='14')

if args.geometry:
    fig = plt.figure()
    ax = fig.gca(projection='3d') #Axes3D(fig)
    ax.hold(True)
    ax.scatter(coordStats[:,0],coordStats[:,1],coordStats[:,2],zdir='z',s=20,c='b')
    if (coordShots.size>3):
        ax.scatter(coordShots[:,0],coordShots[:,1],coordShots[:,2],zdir='z',s=20,c='r',marker='^')
    else:
        ax.scatter(coordShots[0],coordShots[1],coordShots[2],zdir='z',s=200,c='r',marker='^')
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    ax.set_zlim3d(zmin,zmax)
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
    if analytical:
        plt.plot(realP[:,1],z,color=(0,0,0.5),linewidth=4)
        if swaves:
            plt.plot(realS[:,1],z,color=(0,0.5,0),linewidth=4)
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
    if nPriorProfiles > 0:
        plt.figure()
        plt.hold(True)
        for i in np.arange(nPriorProfiles):
            plt.plot(priorP[i][:,1],z)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel(r'P wave velocity profiles from prior ($m.s^{-1}$)',fontsize='14')
        plt.ylabel(r'Depth ($m$)',fontsize='14')
        plt.ylim(ymin=z.max())
        plt.ylim(ymax=z.min())
        if swaves:
            plt.figure()
            plt.hold(True)
            for i in np.arange(nPriorProfiles):
                plt.plot(priorS[i][:,1],z)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.xlabel(r'S wave velocity profiles from prior ($m.s^{-1}$)',fontsize='14')
            plt.ylabel(r'Depth ($m$)',fontsize='14')
            plt.ylim(ymin=z.max())
            plt.ylim(ymax=z.min())
                
        
if args.all:
  chainsToPlot=np.arange(nbt)
elif args.results >= 0:
    if args.results > nbt-1:
        print "There were just ",nbt," chains running!"
        print "-> maximum index : ",nbt-1
        sys.exit(0)
    chainsToPlot=np.array([args.results])

if (args.results >= 0) or args.all:
    lb=qp+"\% confidence interval"
    maxiP=globalMaxP[:,1].max()
    miniP=globalMinP[:,1].min()
    maxiS=globalMaxS[:,1].max()
    miniS=globalMinS[:,1].min()
    dp=(maxiP-miniP)/10
    ds=(maxiS-miniS)/10
    for i in chainsToPlot:
        plt.figure()
        plt.hold(True)
        if args.show_ranges:
            plt.plot(maxP[i][:,1],zFilt,color=(0.4,0.8,0.8),label="Range investigated by chain "+str(i))
            plt.plot(minP[i][:,1],zFilt,color=(0.4,0.8,0.8))
        if not args.dont_show_guess:
            plt.plot(firstGuessP[:,1],z,color=(0.5,0.5,0.95),linewidth=4,label="First guess velocity profile")
        if analytical:
            plt.plot(realP[:,1],z,color=(0,0,0.5),linewidth=4,label="Real velocity profile")
        if not args.no_density_plots:
            pdense(zFilt,(qSupPs[i][:,1]+qInfPs[i][:,1])/2,(qSupPs[i][:,1]-qInfPs[i][:,1])/2,M)
        else:
            plt.plot(qSupPs[i][:,1],zFilt,linestyle='--',color=(0.3,0.3,0.7),label=lb)
            plt.plot(qInfPs[i][:,1],zFilt,linestyle='--',color=(0.3,0.3,0.7))
        if args.show_averages:
            plt.plot(averagesP[i][:,1],zFilt,label="Average model")
     #   plt.plot(averagesP[i][:,1]+np.sqrt(varPs[i][:,1]),zFilt,color=(0.5,0.5,0),label="Standard deviation")
     #   plt.plot(averagesP[i][:,1]-np.sqrt(varPs[i][:,1]),zFilt,color=(0.5,0.5,0))
        if args.show_ranges:
            plt.plot(globalMaxP[:,1],zFilt,color=(1,0,0),label="Range investigated by all chains")
            plt.plot(globalMinP[:,1],zFilt,color=(1,0,0))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.title(r'Chain '+str(i),fontsize='14')
        plt.xlabel(r'P waves velocity ($m.s^{-1}$)',fontsize='14')
        plt.ylabel(r'Depth ($m$)',fontsize='14')
        plt.xlim(miniP-dp, maxiP+dp)  
        plt.legend()
        plt.ylim(ymin=zFilt.max())
        plt.ylim(ymax=zFilt.min())
        if swaves:
            plt.figure()
            if args.show_ranges:
                plt.plot(maxS[i][:,1],zFilt,color=(0.4,0.8,0.8),label="Range investigated by chain "+str(i))
                plt.plot(minS[i][:,1],zFilt,color=(0.4,0.8,0.8))
            if not args.dont_show_guess:
                plt.plot(firstGuessS[:,1],z,color=(0.5,0.95,0.5),linewidth=4,label="First guess velocity profile")
            if analytical:
                plt.plot(realS[:,1],z,color=(0,0.5,0),linewidth=4,label="Real velocity profile")
            if not args.no_density_plots:
                pdense(zFilt,(qSupSs[i][:,1]+qInfSs[i][:,1])/2,(qSupSs[i][:,1]-qInfSs[i][:,1])/2,M)
            else:
                plt.plot(qSupSs[i][:,1],zFilt,linestyle='--',color=(0.3,0.3,0.7),label=lb)
                plt.plot(qInfSs[i][:,1],zFilt,linestyle='--',color=(0.3,0.3,0.7))                
            if args.show_averages:
                plt.plot(averagesS[i][:,1],zFilt,label="Average model")
     #       plt.plot(averagesS[i][:,1]+np.sqrt(varSs[i][:,1]),zFilt,color=(0.5,0.5,0),label="Standard deviation")
      #      plt.plot(averagesS[i][:,1]-np.sqrt(varSs[i][:,1]),zFilt,color=(0.5,0.5,0))

            if args.show_ranges:
                plt.plot(globalMaxS[:,1],zFilt,color=(1,0,0),label="Range investigated by all chains")
                plt.plot(globalMinS[:,1],zFilt,color=(1,0,0))
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.title(r'Chain '+str(i),fontsize='14')
            plt.xlabel(r'S waves velocity ($m.s^{-1}$)',fontsize='14')
            plt.ylabel(r'Depth ($m$)',fontsize='14')
            plt.xlim(miniS-ds, maxiS+ds)  
            plt.legend()
            plt.ylim(ymin=zFilt.max())
            plt.ylim(ymax=zFilt.min())

if args.best:
    if bestModelCalculated:
        diffDataBestModel=bestModelTimes-timesData
        if recalculate_t0 or args.recalculate_t0:
            for i in np.arange(nShots):
                diffPshoti=diffDataBestModel[i*nStats:(i+1)*nStats,0][timesData[i*nStats:(i+1)*nStats,0]>0]
                t0ShotsPi=diffPshoti.mean()
                diffDataBestModel[i*nStats:(i+1)*nStats,0]=diffDataBestModel[i*nStats:(i+1)*nStats,0]-t0ShotsPi
                if swaves:
                    diffSshoti=diffDataBestModel[i*nStats:(i+1)*nStats,1][timesData[i*nStats:(i+1)*nStats,1]>0]
                    t0ShotsSi=diffSshoti.mean()
                    diffDataBestModel[i*nStats:(i+1)*nStats,1]=diffDataBestModel[i*nStats:(i+1)*nStats,1]-t0ShotsSi
        diffP=diffDataBestModel[:,0][timesData[:,0]>0]
        diffS=diffDataBestModel[:,1][timesData[:,1]>0]
        fig = plt.figure()
        plt.hold(True)
        plt.grid(True)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+sigmaP,'b--',linewidth=2)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-sigmaP,'b--',linewidth=2)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+2*sigmaP,'--',color=(0.3,0.3,1),linewidth=1.5)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-2*sigmaP,'--',color=(0.3,0.3,1),linewidth=1.5)
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))+3*sigmaP,'--',color=(0.5,0.5,1))
        plt.plot(np.arange(len(diffP)),np.zeros(len(diffP))-3*sigmaP,'--',color=(0.5,0.5,1))
        plt.plot(np.arange(len(diffP)),diffP,'g+')
        plt.ylim([-20*sigmaP,20*sigmaP])
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
            plt.grid(True)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+sigmaS,'b--',linewidth=2)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-sigmaS,'b--',linewidth=2)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+2*sigmaS,'--',color=(0.3,0.3,1),linewidth=1.5)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-2*sigmaS,'--',color=(0.3,0.3,1),linewidth=1.5)
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))+3*sigmaS,'--',color=(0.5,0.5,1))
            plt.plot(np.arange(len(diffS)),np.zeros(len(diffS))-3*sigmaS,'--',color=(0.5,0.5,1))
            #plt.plot(np.arange(len(diffS)),diffDataBestModel[:,1],'g+')
            plt.plot(np.arange(len(diffS)),diffS,'g+')
            plt.ylim([-20*sigmaP,20*sigmaP])
            plt.xlim([0,len(diffS)-1])
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.text(0.6666*len(diffS), 1.1*sigmaS, r'$1\sigma$',fontsize='30',color='b')
            plt.text(0.7708*len(diffS), 2.1*sigmaS, r'$2\sigma$',fontsize='30',color=(0.3,0.3,1))
            plt.text(0.875*len(diffS), 3.1*sigmaS, r'$3\sigma$',fontsize='30',color=(0.5,0.5,1))
            plt.xlabel(r'Receiver number',fontsize='14')
            plt.ylabel(r'S waves arrival times residuals ($s$)',fontsize='14')
            plt.title(r'Best model residuals',fontsize='14')
            E=sum((diffP/sigmaP)**2/2)+sum((diffS/sigmaS)**2/2)+ep
        else:
            E=sum((diffP/sigmaP)**2/2)+ep
        if args.verbose:
            print "Energy of best model :",E
    else:
        print "The best model has not been calculated"
    nBest=len(glob.glob1(args.pathToDir,"bestPprofile*"))
    if args.verbose:
        print "Number of good model kept : ",nBest
    bestP=[0]*nBest
    bestS=[0]*nBest
    ii=0
    iterationBest=[0]*nBest
    EP=[0]*nBest
    ES=[0]*nBest
    for bestModel in glob.glob1(args.pathToDir,"bestPprofile*"):
        ii=ii+1
        chain=bestModel.split("chain")[1].split(".")[0]
        idx=bestModel.split("idx")[1].split(".")[0]
        iterationBest[ii-1]=int(idx)
        EP[ii-1]=float(bestModel.split("E")[1].split(code)[0].strip("."))
        bestP[ii-1]=np.loadtxt(args.pathToDir+bestModel)
        if args.verbose:
            print "Model number : ",ii, " -> ",bestModel," generated by chain ",chain," at iteration ",idx," (energy "+str(EP[ii-1])+")"
    ii=0
    for bestModel in glob.glob1(args.pathToDir,"bestSprofile*"):
        ii=ii+1
        chain=bestModel.split("chain")[1].split(".")[0]
        idx=bestModel.split("idx")[1].split(".")[0]
        ES[ii-1]=float(bestModel.split("E")[1].split(code)[0].strip("."))
        bestS[ii-1]=np.loadtxt(args.pathToDir+bestModel)
        if args.verbose:
            print "Model number : ",ii, " -> ",bestModel," generated by chain ",chain," at iteration ",idx," (energy "+str(ES[ii-1])+")"
    maxiP=globalMaxP[:,1].max()
    miniP=globalMinP[:,1].min()
    maxiS=globalMaxS[:,1].max()
    miniS=globalMinS[:,1].min()
    dp=(maxiP-miniP)/10
    ds=(maxiS-miniS)/10
    from operator import itemgetter
    idxBestP=min(enumerate(EP), key=itemgetter(1))[0] # index of best model
    idxBestS=min(enumerate(ES), key=itemgetter(1))[0] # index of best S model (it is the same one!)
    plt.figure()
    if analytical:
        plt.plot(realP[:,1],z,color=(0,0,0.5),linewidth=4,label="Real velocity profile")
    if args.verbose:
        print "Models kept after iteration : "+str(args.treshold)+" will be shown"

    for i in np.arange(nBest):
        if iterationBest[i] > args.treshold:
            plt.hold(True)
            plt.grid(True)
            if i == idxBestP:
                plt.plot(bestP[i][:,1],zFilt,linewidth=4,label="Best model")
            else:
                plt.plot(bestP[i][:,1],zFilt)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.xlabel(r'Best P wave velocity models in ($m.s^{-1}$)',fontsize='14')
            plt.ylabel(r'Depth ($m$)',fontsize='14')
            plt.xlim(miniP-dp, maxiP+dp)
    plt.ylim(ymax=z.max())
    plt.gca().invert_yaxis()
    plt.legend()
    if swaves:    
        plt.figure()
        if analytical:
            plt.plot(realS[:,1],z,color=(0,0,0.5),linewidth=4,label="Real velocity profile")
        for i in np.arange(nBest):
            if iterationBest[i] > args.treshold:
                plt.hold(True)
                plt.grid(True)
                if i == idxBestS:
                    plt.plot(bestS[i][:,1],zFilt,linewidth=4,label="Best model")
                else:
                    plt.plot(bestS[i][:,1],zFilt)
                plt.rc('text', usetex=True)
                plt.rc('font', family='serif')
                plt.xlabel(r'Best S wave velocity models in ($m.s^{-1}$)',fontsize='14')
                plt.ylabel(r'Depth ($m$)',fontsize='14')
                plt.xlim(miniS-ds, maxiS+ds)
        plt.ylim(ymax=z.max())
        plt.gca().invert_yaxis()
    plt.legend()


################################ VP/VS ################################
if args.vpvs:
    if not swaves:
        print "Impossible to print Vp/Vs ratio as Swaves had not been calculated during the algorithm"
    else:
        for i in np.arange(nbt):
            plt.figure()
            plt.hold(True)
            plt.grid(True)
            vpFractionalUncertainty = 100*np.sqrt(varPs[i][:,1])/averagesP[i][:,1] # Ex: is vp=2500+/-100m/s -> fractional uncertainty 4% vp = 2500+/-4%
            vsFractionalUncertainty = 100*np.sqrt(varSs[i][:,1])/averagesS[i][:,1] # Ex: is vs=2500+/-100m/s -> fractional uncertainty 4% vs = 2500+/-4%
            ratioFractionalUncertainty = vpFractionalUncertainty + vsFractionalUncertainty
            meanRatio = averagesP[i][:,1]/averagesS[i][:,1]
            numericalUncertainty = meanRatio*(vpFractionalUncertainty + vsFractionalUncertainty)/100
            plt.plot(meanRatio + numericalUncertainty,zFilt,linestyle='--',color=(0.3,0.3,0.7),label="Standard deviation")
            plt.plot(meanRatio - numericalUncertainty,zFilt,linestyle='--',color=(0.3,0.3,0.7))
            plt.plot(meanRatio,zFilt,label="Average Vp/Vs for chain "+str(i))
            plt.plot(meanRatio,zFilt,label="Average Vp/Vs for chain "+str(i))
            # Load best profiles :
            nBest=len(glob.glob1(args.pathToDir,"bestPprofile*"))
            bestP=[0]*nBest
            bestS=[0]*nBest
            EP=[0]*nBest
            ES=[0]*nBest
            ii=0
            for bestModel in glob.glob1(args.pathToDir,"bestPprofile*"):
                ii=ii+1
                EP[ii-1]=float(bestModel.split("E")[1].split(code)[0].strip("."))
                bestP[ii-1]=np.loadtxt(args.pathToDir+bestModel)
            ii=0
            for bestModel in glob.glob1(args.pathToDir,"bestSprofile*"):
                ii=ii+1
                ES[ii-1]=float(bestModel.split("E")[1].split(code)[0].strip("."))
                bestS[ii-1]=np.loadtxt(args.pathToDir+bestModel)
            from operator import itemgetter
            idxBestP=min(enumerate(EP), key=itemgetter(1))[0] # index of best model
            idxBestS=min(enumerate(ES), key=itemgetter(1))[0] # index of best S model (it is the same one!)
            # End of loading best profiles
            plt.plot(bestP[idxBestP][:,1]/bestS[idxBestS][:,1],zFilt,linewidth=4,label="Vp/Vs of the best model")
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.xlabel(r'Ratio Vp/Vs',fontsize='14')
            plt.ylabel(r'Depth ($m$)',fontsize='14')
            plt.ylim(ymax=zFilt.max())
            plt.legend()

################################ SWAPS ################################
if args.swaps:
    print "Not implemented for now"
#    plt.figure()
#   plt.plot(ll[ll[:,2]==3,1],ll[ll[:,2]==3,3])
#plot(exch(exch(:,2)==3,1),exch(exch(:,2)==3,3),'k')
# plot(exch(exch(:,2)==3,1),exch(exch(:,2)==3,3),'y')
# plot(exch(exch(:,2)==2,1),exch(exch(:,2)==2,3),'c')
#plot(exch(exch(:,2)==1,1),exch(exch(:,2)==1,3),'r')
#plot(exch(exch(:,2)==0,1),exch(exch(:,2)==0,3),'b')

plt.show()

    
