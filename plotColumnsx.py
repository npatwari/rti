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
#
# Version History:
#
# Version 1.0:  Initial Release
#
# Version 1.1:  27 Sept 2012.
#   Change: This version won't delete (pop) data unless it is going to enter 
#           in a new piece of data.
#
# Version x.0:  3 July 2016.
#               This version plots data from python listenx.py, which is meant
#               to gather data from the XANDEM HOME gateway.




import sys
import time
from collections import deque
import numpy as np
import matplotlib.pyplot as plt
import getopt, ast
import rss

# Accept command line inputs to list the links to plot
linkCombosToPlot = []
numNodes         = 10
myopts, args = getopt.getopt(sys.argv[1:],"l:n:")
for o, a in myopts:
    if o == "-l":
        inputList = ast.literal_eval( a )
        linkCombosToPlot.append(inputList)
    elif o == "-n":
        numNodes = int(a)

# If no number of nodes given, 
if numNodes <= 1:
    sys.stderr.write('Usage: %s -n numberOfNodes -l "[tx,rx,ch]" ' % sys.argv[0])
    sys.stderr.write('numberOfNodes is required and must be greater than 1')

# If no link combination provided, just plot [tx=1,rx=2,ch=0]
if linkCombosToPlot == []:
    linkCombosToPlot.append([1,2,0])

print "Plotting columns for links:"
print linkCombosToPlot
print "Number of nodes:" + str(numNodes)

# Parameters you may change:
#   plotSkip:  Refresh the plot after this many data lines are read
#   buffL:     buffer length, ie, how much data to plot.
#   startSkip: A serial port seems to have a "memory" of several lines, 
#              which were saved from the previous experiment run.  
#              ** Must be greater than 0.
plotSkip    = 1
buffL       = 40
startSkip   = 1
min_range   = 5

# remove junk from start of file.
for i in range(startSkip): 
  line = sys.stdin.readline()

# Use the most recent line to determine how many columns (streams) there are.
lineInt = [float(i) for i in line.split()]
columns = len(lineInt)
rss_dB  = [int(f) for f in lineInt[:-1]]  # take all columns except for last column
timeSec = lineInt[-1]/1000.0   # Convert from ms to sec
numLinks = len(rss_dB)
numChs  = numLinks / ((numNodes-1)*numNodes)
nodeList = range(1,numNodes+1)
channelList = range(numChs)

streams = len(linkCombosToPlot)
linkNumList = []
for l in linkCombosToPlot:
    linkNumList.append( rss.linkNumForTxRxChLists(l[0], l[1], l[2], nodeList, channelList) )
    

# Init the figure.
plt.ion()
RSSBuffer   = []
TimeBuffer   = []
linePlot    = []
plt.cla()
for i,linkNum in enumerate(linkNumList):
    RSSBuffer.append( deque([rss_dB[i]] * buffL))
    TimeBuffer.append( deque([timeSec] * buffL))
    l, = plt.plot([0]*buffL, RSSBuffer[i], label=str(linkCombosToPlot[i]))
    plt.ylim((-95, -25))
    plt.grid("on")
    plt.ylabel('Received Power (dBm)')
    plt.xlabel('Measurement Time Ago (sec)')
    linePlot.append(l)

# Run forever, adding lines as they are available.
counter = 0
while 1:
    print "counter = " + str(counter)
    line = sys.stdin.readline()
    if not line:
        continue
    while line[-1] != '\n':   # If line is incomplete, add to it.
        line += sys.stdin.readline()
        if not line:
            continue
        
    # Get the integers from the line string
    data = [float(i) for i in line.split()]
    rss_dB  = [int(f) for f in data[:-1]]  # take columns except for last (time) col
    timeSec = data[-1]/1000.0  # Convert from ms to sec
    
    # Append the queue for each stream with the newest RSS and Timestamps
    for i,linkNum in enumerate(linkNumList):
        # data > -10 indicates no data measured.  Don't include a new value.
        if (rss_dB[linkNum] < -10):  
            oldRSS = RSSBuffer[i].popleft()
            oldTS  = TimeBuffer[i].popleft()
            RSSBuffer[i].append(rss_dB[linkNum])
            TimeBuffer[i].append(timeSec)
            
    # Every plotSkip rows, redraw the plot.  Time stamps, relative to the 
    #   maximum time stamp, are on the x-axis.
    counter += 1
    if np.mod(counter, plotSkip) == 0:
        mintime = 0
        maxv    = -200.0  # Lower than I'd expect any RSS to be (dBm)
        minv    = 0.0  # Higher than any RSS could be (dBm)
        for i in range(streams):
            linePlot[i].set_ydata(RSSBuffer[i])
            relTime = np.array(TimeBuffer[i]) - max(TimeBuffer[i])
            linePlot[i].set_xdata(relTime)
            mintime = min(min(relTime), mintime)    # Keep track of lowest reltime
            minv    = min(min(RSSBuffer[i]), minv)  # track lowest rss value
            maxv    = max(max(RSSBuffer[i]), maxv)  # track highest rss value
            print minv, maxv
        
        # Make the minimum and max a multiple of min_range so that axes 
        # don't change as often
        minscale = rss.floor_multiple_of(min_range,minv)
        maxscale = rss.ceil_multiple_of(min_range,maxv)
        print "[" + str(minscale) + ", " + str(maxscale) + "]"
        plt.axis([mintime, 0, minscale, maxscale])
        plt.draw()