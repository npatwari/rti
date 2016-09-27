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
# Version 1.0:  Initial Release 20 Sept 2016
#




import rti
import sys
import numpy as np
import matplotlib.pyplot as plt
import getopt, ast
import rss

# Default inputs if not provided on command line
linkCombosToPlot = []
numNodes         = 10
infile           = sys.stdin
fromStdin        = True

# Accept command line inputs to list the links to use; the input file
myopts, args = getopt.getopt(sys.argv[1:],"l:n:f:")
for o, a in myopts:
    if o == "-l":
        inputList = ast.literal_eval( a )
        linkCombosToPlot.append(inputList)
    elif o == "-n":
        numNodes = int(a)
    elif o == "-f":
        infile = open(a)
        fromStdin = False


# If no number of nodes given, 
if numNodes <= 1:
    sys.stderr.write('Usage: %s -n numberOfNodes -l "[tx,rx,ch]" ' % sys.argv[0])
    sys.stderr.write('numberOfNodes is required and must be greater than 1')

# If no link combination provided, just plot [tx=1,rx=2,ch=0]
if linkCombosToPlot == []:
    linkCombosToPlot.append([1,2,0])
links       = len(linkCombosToPlot)

print "Plotting columns for links:"
print linkCombosToPlot
print "Number of nodes:" + str(numNodes)

# Parameters you may change:
#   plotSkip:  Refresh the plot after this many data lines are read
#   buffL:     buffer length, ie, how much data to plot.
#   startSkip: A serial port seems to have a "memory" of several lines, 
#              which were saved from the previous experiment run.  
#              ** Must be greater than 0.
buffL       = 4
startSkip   = 1

# remove junk from start of file.
for i in range(startSkip): 
  line = infile.readline()

# Use the most recent line to determine how many columns (streams) there are.
lineInt = [float(i) for i in line.split()]
columns = len(lineInt)
rss_dB  = [int(f) for f in lineInt[:-1]]  # take all columns except for last column
numLinks = len(rss_dB)
numChs  = numLinks / ((numNodes-1)*numNodes)
nodeList = range(1,numNodes+1)
channelList = range(numChs)

# Convert the [tx,rx,ch] combos into column numbers
linkNumList = []
for l in linkCombosToPlot:
    linkNumList.append( rss.linkNumForTxRxChLists(l[0], l[1], l[2], nodeList, channelList) )

buff = []
for i in range(links):
    buff.append( rti.FixedLenBuffer([0]*buffL))


# Run forever, adding lines as they are available.
counter = 0
keepReading    = True
while keepReading:
    print "counter = " + str(counter)
    line = infile.readline()
    # If at the "end of file", keep reading if reading from stdin.
    if not line:
        keepReading = fromStdin
        continue
    while line[-1] != '\n':   # If line is incomplete, add to it.
        line += infile.readline()
        
    # Get the integers from the line string
    data = [float(i) for i in line.split()]
    rss_dB  = [int(f) for f in data[:-1]]  # take columns except for last (time) col
    
    # Append the queue for each stream with the newest RSS and Timestamps
    for i,linkNum in enumerate(linkNumList):
        # data > -10 indicates no data measured.  Don't include a new value.
        if (rss_dB[linkNum] < -10):  
            buff[i].append(rss_dB[linkNum])

    # Calculate the variance on all monitored links
    varVec       = np.array([a.var() for a in buff])
    print varVec

    counter += 1
