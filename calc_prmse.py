#! /usr/bin/env python

#
# Version History:
#
# Version 1.0:  Initial Release.  17 Oct 2012.
# Version 2.0:  27 Oct 2012.
# Version 2.01: 18 Oct 2016: modified for Neal's basement experiment data.
#
# Usage: python calc_prmse_v2 coordEstFile.txt
#    where coordEstFile.txt is the output of your RTI algorithm.
#    coordEstFile is assumed to be the same number of rows as the 
#    RSS file (the output of listenAllLinks.py).  
# 

import rti
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

    

penalty       = 4.0**2.
noPersonKey   = -99.0
pivotFileName = 'basement/pivot_coords_basement_m.txt'
pivotCoords   = np.loadtxt(pivotFileName)
pathFileName  = 'basement/path_basement_1_f.txt'  # list of pivot numbers in order.
pathInd       = np.loadtxt(pathFileName)
startPathTime = 56000.0  # ms.  I know I started at time 14*4
speed         = 1.0 / 8000.0 # pivot points per millisecond.
eps           = 0.01

coordFileName = 'basement/sensor_coords_basement_m.txt'  
sensorCoords  = np.loadtxt(coordFileName)

if len(sys.argv) >= 2:
    estFileName      = sys.argv[1]
else:
    print("Error: no estimated coordinate filename given.")

dataFileName  = 'basement/basement_listenx_out_1.txt'
dataMat       = np.loadtxt(dataFileName)
[rows,cols]   = dataMat.shape
time_ms       = dataMat[:,-1]

estCoord      = np.loadtxt(estFileName)
if (estCoord.shape[0] != rows):
    sys.exit("Error: estimate coordinate file is not the correct number of rows: should be " + str(rows) + " rows, but is " + str(estCoord.shape[0]))
if (estCoord.shape[1] != 2):
    sys.exit("Error: estimate coordinate file is not the correct number of columns")

# We know the algorithm has some delay (or advance, if it has 
# over-compensated for the delay).  Compute the error for a range of time
# offsets.  Negative offset means the estimate is delayed compared to the actual.
offsetList    = np.arange(-1000.0,1001.0,250)
offsets       = len(offsetList)
actualCoord   = []
err           = np.zeros(len(offsetList))
for j,offset in enumerate(offsetList):
    for i in range(rows):
        actualCoord.append(np.zeros((rows, 2)))
        temp = rti.calcActualPosition(time_ms[i] + offset, pivotCoords, pathInd, startPathTime, speed)
        if len(temp)==0:
            actualCoord[j][i,:] = [noPersonKey, noPersonKey]
        else:
            actualCoord[j][i,:] = temp
    
    #print str(actualCoord-estCoord) 
    err[j]  = rti.prmse(actualCoord[j], estCoord, noPersonKey, penalty)

print str(err)
bestOI = np.argmin(err)
minPRMSE = err[bestOI]

print str(minPRMSE) + " pRMSE at offset = " + str(offsetList[bestOI])

plt.ion()
fig = plt.figure(num=1, figsize=(7, 7))
plt.cla()
rti.plotLocs(sensorCoords)
#imhandle = plt.imshow(image, interpolation='none', origin='lower', 
#                      extent=imageExtent, vmin=0, vmax=8)
plt.ylabel('Y Coordinate (ft)')
plt.xlabel('X Coordinate (ft)')
step = 3
for i in range(0, rows, step):
    ac = actualCoord[bestOI][i,:]
    ec = estCoord[i,:]
    noP = np.abs(ac[0] - noPersonKey) < eps
    estNoP = np.abs(ec[0] - noPersonKey) < eps
    if not(noP) and not(estNoP):
        plt.text(ac[0], ac[1],'X', horizontalalignment='center',
             verticalalignment='center', color='k')
        plt.text(ec[0], ec[1],'o', horizontalalignment='center',
             verticalalignment='center', color='k')
        plt.plot([ac[0], ec[0]], [ac[1], ec[1]],'r')
    elif not(noP) and estNoP:
        plt.text(ac[0], ac[1],'X', horizontalalignment='center',
             verticalalignment='center', color='r')
    elif noP and not(estNoP):
        plt.text(ec[0], ec[1],'o', horizontalalignment='center',
             verticalalignment='center', color='r')

plt.axes().set_aspect('equal', 'datalim')
#axis([-0.5, 25., -0.5, 24])
plt.draw()

# remove the filename extension, add "2.png" and save it.
outfname = estFileName[:-4] + "2.png"
print outfname
fig = fig.savefig(outfname)

raw_input()
