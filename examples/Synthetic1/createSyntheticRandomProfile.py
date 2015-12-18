# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 21:31:53 2015

Create random synthetic velocity profile + linear first guesses

@author: alex
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

plt.close('all')
random.seed(2)
X = range(256)
Z = np.linspace(0,1000,256)
pos = 50
posVpVs = 1.7
L = np.array([pos])
vpvs = np.array([posVpVs])
for x in X[1:]:
    pos += random.choice((-0.9,1)) #random.choice((-1,1))
    posVpVs += random.choice((-0.02,0.02))
    L=np.append(L,pos)
    vpvs = np.append(vpvs,posVpVs)

L=70*L
Vp = savitzky_golay(L, 51, 3) # window size 51, polynomial order 3
A = np.array([ Z, np.ones(256)])
Vs = Vp/savitzky_golay(vpvs, 51, 3) # window size 51, polynomial order 3
w = np.linalg.lstsq(A.T,Vp)[0] # obtaining the parameters
# plotting the line
lineP = w[0]*Z+w[1]+500 # regression line
w = np.linalg.lstsq(A.T,Vs)[0] # obtaining the parameters
# plotting the line
lineS = w[0]*Z+w[1]-250 # regression line

plt.figure()
plt.hold(True)
plt.plot(L,Z,label="Random walk")
plt.plot(Vp,Z,linewidth=4,label="P wave velocity from this random walk")
plt.plot(lineP,Z,linewidth=4,label="First guess")
ax = plt.axes()
ax.set_ylim(Z[0],Z[-1])
ax.invert_yaxis()
plt.legend()

plt.figure()
plt.hold(True)
plt.plot(vpvs,Z,linewidth=4,label="Random walk vp/vs")
plt.legend()
ax = plt.axes()
ax.set_ylim(Z[0],Z[-1])
ax.invert_yaxis()

plt.figure()
plt.hold(True)
plt.plot(Vs,Z,linewidth=4,label="S wave velocity from random vp/vs")
plt.plot(lineS,Z,linewidth=4,label="First guess")
plt.legend()
ax = plt.axes()
ax.set_ylim(Z[0],Z[-1])
ax.invert_yaxis()

# Save profiles
np.savetxt("dataExample/realProfileP.txt",np.dstack((Z,Vp))[0])
np.savetxt("dataExample/realProfileS.txt",np.dstack((Z,Vs))[0])
np.savetxt("dataExample/firstGuessP.txt",np.dstack((Z,lineP))[0])
np.savetxt("dataExample/firstGuessS.txt",np.dstack((Z,lineS))[0])

#####################################################################

coordShotsX=[300,500]
coordShotsY=[400]
coordShotsZ=[650]

coordStatsX=[200,300,400,500,600]
coordStatsY=[200,300,400,500,600]
coordStatsZ=[200,300,400,500,600]

Xshots=[]
Yshots=[]
Zshots=[]
Xstats=[]
Ystats=[]
Zstats=[]

#Open a file in write mode:
fo = open("dataExample/coordShots.txt", "w+")
for coordX in coordShotsX:
    for coordY in coordShotsY:
        for coordZ in coordShotsZ:
            Xshots.append(coordX)
            Yshots.append(coordY)
            Zshots.append(coordZ)
            fo.write(str(coordX)+" "+str(coordY)+" "+str(coordZ)+"\n")
# Close opened file
fo.close()

#Open a file in write mode:
fo = open("dataExample/coordStats.txt", "w+")
for coordX in coordStatsX:
    for coordY in coordStatsY:
        for coordZ in coordStatsZ:
            Xstats.append(coordX)
            Ystats.append(coordY)
            Zstats.append(coordZ)
            fo.write(str(coordX)+" "+str(coordY)+" "+str(coordZ)+"\n")
# Close opened file
fo.close()

fig = plt.figure()
ax = fig.gca(projection='3d') #Axes3D(fig)
ax.hold(True)
ax.scatter(Xstats,Ystats,Zstats,zdir='z',s=20,c='b')
if (len(coordShotsX) > 3):
    ax.scatter(Xshots,Yshots,Zshots,zdir='z',s=20,c='r',marker='^')
else:
    ax.scatter(Xshots,Yshots,Zshots,zdir='z',s=200,c='r',marker='^')
ax.set_xlim3d(min(min(Xshots),min(Xstats))-100,max(max(Xshots),max(Xstats))+100)
ax.set_ylim3d(min(min(Yshots),min(Ystats))-100,max(max(Yshots),max(Ystats))+100)
ax.set_zlim3d(min(min(Zshots),min(Zstats))-100,max(max(Zshots),max(Zstats))+100)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Geometry')
ax.invert_zaxis()

