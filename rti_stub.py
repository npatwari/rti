#! /usr/bin/env python

#
# LICENSE:
# Copyright (C) 2016  Neal Patwari
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Author: Neal Patwari, neal.patwari@gmail.com
#
# Version History:
#
# Version 1.0:  Initial Release.  17 Oct 2012
# Version 2.0:  30 Oct 2014.
# Version 2.1:  3 July 2016.
#               Added actualKnown option to turn off location error calculation
#
# Purpose:  
#   This script calls functions in order to calculate radio tomographic images
#   from the following inputs: 
#       1) the locations of radio sensors; 
#       2) received signal strength (RSS) measurements between pairs of those 
#          sensors, on multiple channels.  This RTI code assumes that the RSS 
#          was recorded using listenAllLinks.py or listenx.py, which simply 
#          lists the RSS value on every (Tx, Rx, Channel) combination.  The last 
#          column is the time, in integer milliseconds from the start of the 
#          data collection.  
#   The code plots the RTI estimate.  Currently, RTI methods implemented are:
#       1) multi-channel attenuation-based (traditional) RTI.  This is the method
#          published in MASS 2012: 
#             O. Kaltiokallio, M. Bocca, and N. Patwari, Enhancing the accuracy 
#          of radio tomographic imaging using channel diversity, 9th IEEE Int'l 
#          Conference on Mobile Ad hoc and Sensor Systems, October 8-11, 2012.
#       2) multi-channel variance-based RTI.  This is the straightforward
#          extension of VRTI to multiple channels. Single-channel VRTI is 
#          described in: 
#             J. Wilson and N. Patwari, See Through Walls: Motion Tracking Using 
#          Variance-Based Radio Tomography Networks, IEEE Transactions on Mobile 
#          Computing, vol. 10, no. 5, pp. 612-621, May 2011.
#   The code outputs coordinate estimates into the file, outputFileName.
#
# Usage: 
#   cat filename.txt | python rti_stub.py
#   where filename.txt is the output listenAllLinks.py or listenx.py
#
#   Alternatively, use 
#      python listenx.py -i filename.txt | python rti_stub_v2.py 
#   if the filename.txt was recorded at the gateway node.
#
#   Alternatively, use
#      ssh root@xandem-gateway.local "/opt/xandev/exec/gateway/bin/gateway -l" | python listenx.py | python rti_stub_v2.py 
#   for a real-time operation of the script from the data collected at the Xandem gateway node.  This assumes your gateway is connected to the same local network as the computer running this script.
#   
#   There must be two files in the directory to use actualKnown == True:
#      pivotFileName
#      pathFileName
#
#   There must always be the true sensor coordinate file in the directory:
#      coordFileName
#
#   The output coordinate estimates (one per row of the input file) 
#   is saved to the file named in outputFileName.  A row of -99 -99 indicates 
#   that no person is estimated to be present in the area.


import rti
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 


# Parameters you may change:
#   plotSkip:  Refresh the plot after this many data lines are read
#   buffL:     buffer length, ie, how much recent data to save for each link.
#   startSkip: A serial port seems to have a "memory" of several lines, 
#              which were saved from the previous experiment run.  
#              ** Must be greater than 0.
plotSkip      = 2
startSkip     = 1
buffL         = 4  
calLines      = 50
topChs        = 3  # How many channels to include in multichannel RTI
channels      = 8  # How many channels are measured on each TX,RX combination.
delta_p       = 0.2  # distance units from coordFileName
sigmax2       = 0.5  # image value^2
delta         = 1.0  # distance units from coordFileName
excessPathLen = 0.1 # distance units from coordFileName
units         = 'm' # distance units for plot label
actualKnown   = True
outputFileName= "basement/neals_estimate.txt"

# An image max value above this makes us believe someone is in the area.
personInAreaThreshold = 2.1

# Pivot coords are the position of known "pivot spots" on the floor
# which are the only places the person may turn during their travel
# Pivot spots are numbered (starting from 0)
if actualKnown:
    pivotFileName = 'basement/pivot_coords_basement_m.txt'
    pivotCoords   = np.loadtxt(pivotFileName)

    # The "path" is just a list of the number of the pivot spots, in order
    # that the person traverses them.  The person is assumed to 
    # hit pivot spots at a constant rate.
    pathFileName  = 'basement/path_basement_1_f.txt'  # list of pivot numbers in order.
    pathInd       = np.loadtxt(pathFileName)

    startPathTime = 56000.0  # ms.  I know I started at time 14*4
    speed         = 1.0 / 8000.0 # pivot points per millisecond.

# Load the coordinate file, find the bottom left and top right
coordFileName = 'basement/sensor_coords_basement_m.txt'  
sensorCoords  = np.loadtxt(coordFileName)
sensors       = len(sensorCoords)

# Open a file to write the coordinate estimates to
fout          = open(outputFileName, "w")

# It looks nice to have the sensor coordinates plotted on top of any image.
plt.ion()

# initialize the RTI projection matrix, create a rectangular grid of pixels.
inversion, xVals, yVals = rti.initRTI(sensorCoords, delta_p, sigmax2, delta, excessPathLen)
imageExtent   = (min(xVals) - delta_p/2, max(xVals) + delta_p/2, 
                 min(yVals) - delta_p/2, max(yVals) + delta_p/2)  # the min, max for plot axes
xValsLen      = len(xVals)  # how many pixels along the x-axis.
yValsLen      = len(yVals)  # how many pixels along the y-axis.

# Open the file in argv[1]; if none, use stdin as the file
if len(sys.argv) == 1:
    infile = sys.stdin
    fromStdin = True
else:
    infile = open(sys.argv[1])
    fromStdin = False



# remove junk from start of file.
for i in range(startSkip): 
    line = infile.readline()

# Use the most recent line to determine how many columns (streams) there are.
# The first line is used as the "prevRSS" when reading the following lines.
lineList = [int(i) for i in line.split()]
time_ms  = lineList.pop()  # takes time in ms from the last column.
prevRSS  = np.array(lineList)
numLinks = len(prevRSS)
numPairs = sensors*(sensors-1)
if numLinks != numPairs*channels:
    sys.exit('Error: numLinks = ' + str(numLinks) + '; sensors = ' + str(sensors) + '; channels = ' + str(channels))



# Initialize RSS Buffer, a list of buffers, one for each link.
# For VRTI.
buff = []
for i in range(numLinks):
    buff.append( rti.FixedLenBuffer([0]*buffL))

# Run forever, adding lines as they are available.
counter        = 0
actualCoord    = []
VRTI_err_list  = []
RTI_err_list   = []
sumRSS         = np.zeros(numLinks)
countCalLines  = np.zeros(numLinks)
keepReading    = True
while keepReading:
    print "Counter = " + str(counter)
    line = infile.readline()
    # If at the "end of file", keep reading if reading from stdin.
    if not line:
        keepReading = fromStdin
        continue
    while line[-1] != '\n':   # If line is incomplete, add to it.
        line += infile.readline()
        
    # Get the integers from the line string
    lineList    = [int(i) for i in line.split()]
    time_ms     = lineList.pop(-1)  # remove last element
    rss         = np.array(lineList)
    #print "length of rss = " + str(len(rss))
    if actualKnown:
        actualCoord = rti.calcActualPosition(time_ms, pivotCoords, pathInd, startPathTime, speed)

    # data > -10 means no data measured. Replace with most recently measured value.
    for i in range(numLinks):        
        if (rss[i] > -10):  
            rss[i] = prevRSS[i]
    
        # Store current RSS vector for each link in its FixedLenBuffer
        # For variance-based RTI.
        buff[i].append(rss[i])
    
    # Use first "calLines" data vectors to find average RSS
    if counter < calLines:
        for i in range(numLinks):
            if rss[i] < -10:
                sumRSS[i] += rss[i]
                countCalLines[i] += 1
        fout.write("-99 -99\n")
    
    # At the end of the calLines period, compute the average RSS for each link.
    elif counter == calLines:
        #meanRSS = sumRSS / countCalLines
        # If you have a divide-by-zero problem with the line above, use:
        meanRSS = np.array([s / max(1,countCalLines[i]) for i,s in enumerate(sumRSS)])
        # Sort the meanRSS to decide which channels have the highest average RSS.
        # Make each channel's RSS vector into a separate row
        meanRSS.shape = (channels, numPairs)
        maxInds = meanRSS.transpose().argsort()
        # Sum the highest topChs channels.
        calVec  = rti.sumTopRows(meanRSS, maxInds, topChs)
        # We assume that no person is present during the calibration period.
        fout.write("-99 -99\n")
        
    # When the calibration period is over, calc & display RT images, and write 
    # the current coordinate estimate to the output file.
    if counter >= calLines:
        
        print "RSS on link 1 = " + str(rss[0])
        
        # Compute difference between calVec and current RSS measurement
        rss.shape    = (channels, numPairs)
        curVec       = rti.sumTopRows(rss, maxInds, topChs)
        rss.shape    = numLinks
        scoreVec     = calVec - curVec
        
        # Compute shadowing-based radio tomographic image
        image        = rti.callRTI(scoreVec, inversion, len(xVals), len(yVals))
        RTIMaxCoord  = rti.imageMaxCoord(image, xVals, yVals)
        
        if image.max() > personInAreaThreshold:
            fout.write(str(RTIMaxCoord[0]) + " " + str(RTIMaxCoord[1]) + "\n")
        else:
            fout.write("-99 -99\n")
        #print "RTI Image range:" + str(image.min()) + " - " + str(image.max())
        
        
        # Plot the RTI image each plotSkip data lines.
        if counter % plotSkip == 0:
            if actualKnown:
                rti.plotImage(image, 2, sensorCoords, imageExtent, 8.0, units, time_ms, actualCoord)
            else:
                rti.plotImage(image, 2, sensorCoords, imageExtent, 8.0, units, time_ms)                
            # You must call colorbar() only once, otherwise you get multiple bars.
            if counter==calLines:
                plt.colorbar()

        # Calculate variance on each link and calculate VRTI image
        varVec       = np.array([a.var() for a in buff])
        varVec.shape = (channels, numPairs)
        scoreVec     = rti.sumTopRows(varVec, maxInds, topChs)
    
        # Compute variance-based radio tomographic image
        image        = rti.callRTI(scoreVec, inversion, len(xVals), len(yVals))
        VRTIMaxCoord = rti.imageMaxCoord(image, xVals, yVals)
        # print "VRTI Image range:" + str(image.min()) + " - " + str(image.max())
        if counter % plotSkip == 0:
            rti.plotImage(image, 3, sensorCoords, imageExtent, 16.0, units, time_ms, actualCoord)
            plt.pause(0.001)
        
        # You must call colorbar() only once, otherwise you get multiple bars.
        if counter==calLines:
            plt.colorbar()
  
        # Find the coordinate error
        if (actualKnown):
            if (len(actualCoord) > 0):
                VRTI_err  = np.sqrt(np.sum((actualCoord-VRTIMaxCoord)**2))
                VRTI_err_list.append(VRTI_err)
                RTI_err   = np.sqrt(np.sum((actualCoord-RTIMaxCoord)**2))
                RTI_err_list.append(RTI_err)
            
    # Save RSS in case next line has missing data.
    prevRSS = rss.copy()
    counter += 1

# Overall stats of RMSE over all times
if actualKnown:
    VRTI_RMSE = np.sqrt(np.mean(np.array(VRTI_err_list)**2))
    RTI_RMSE  = np.sqrt(np.mean(np.array(RTI_err_list)**2))
    print "VRTI RMSE = " + str(VRTI_RMSE)
    print "RTI RMSE = " + str(RTI_RMSE)


