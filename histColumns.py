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
# Version 1.0:  15 Sept 2014: Initial Release
#
# Version 2.0: 13 Sept 2016: Update to include particular command line options:
#     - plot particular channels rather than all;
#     - change number of channels or number of nodes via command line;
#     - generally to use the getopt, ast packages for command line inputs;
#     - allow up to eight channels.
#    
#
# Usage:
#     python histColumns.py -f file1.txt -f file2.txt ... -f filen.txt
#             -l "[rx1,tx1,ch1]" -l "[rx2,tx2,ch2]" ... -l "[rxn,txn,chn]"
#
#  where file#.txt was created using, for example:
#     python listenSomeColumns_v2.py > file#.txt
#  where rx#, tx#, ch# are the receiver id, transmitter id, and channel number 
#  for a link that you want to plot RSS data for. 
#
# Purpose: Plot a histogram of RSS for particular columns (links)
#     measured.  Compare histograms when using with multiple files,
#     in separate histograms stacked on top of each other.

import sys
import numpy as np
import matplotlib.pyplot as plt
import getopt, ast
import rss

# Accept command line inputs to list the links to plot
linkCombosToPlot = []
fileList         = []
numNodes         = 10   # Xandem kits come with 10 nodes
numChs           = 8    # Xandem kits come programmed with 8 channels
myopts, args = getopt.getopt(sys.argv[1:],"l:n:f:c:")
for o, a in myopts:
    if o == "-l":
        inputList = ast.literal_eval( a )
        linkCombosToPlot.append(inputList)
    elif o == "-n":
        numNodes = int(a)
    elif o == "-c":
        numChs = int(a)
    elif o == "-f":
        fileList.append(a)

files       = len(fileList)
streams     = len(linkCombosToPlot)
if files == 0:
    sys.exit('Usage: histColumns.py -f file1.txt -f file2.txt ... -l "[1,2,3]" ...')

print "Plotting each file in a subplot from top to bottom:" 
print fileList
print "Plotting links:"
print linkCombosToPlot

# Convert the [tx,rx,ch] combos into column numbers
linkNumList = []
nodeList = range(1,numNodes+1)
channelList = range(numChs)
for l in linkCombosToPlot:
    linkNumList.append( rss.linkNumForTxRxChLists(l[0], l[1], l[2], nodeList, channelList) )


# Inputs / Settings
startSkip   = 0
markerlist  = ['-o', '-s', '-v', '-*', '-8','-p','-+','-x']
xbincenters = range(-110, -10)   # lowest to highest POSSIBLE rss

# Create the bin edges.
xbinedges   = [x - 0.5 for x in xbincenters]
xbinedges.append(xbincenters[-1] + 0.5)

# Init the plot.
plt.ion()
plt.cla()
junk        = plt.figure()
fig         = plt.figure()

# Keep track of min and max of data.  Initialize by giving
# absurdly low max value and absurdly high min value.
minrss = 0
maxrss = -110

# Load data from each file.
for i, file in enumerate(fileList):
    data       = np.loadtxt(file, dtype=float, delimiter=' ', skiprows=startSkip)
    cols       = data.shape[1]
    
    # Remove 127 values from, and plot a histogram, for each column
    good_data  = [ data[data[:,k] < 127, k] for k in linkNumList]
    
    # Update the min and max
    maxrss = max(maxrss, max([gd.max() for gd in good_data]))
    minrss = min(minrss, min([gd.min() for gd in good_data]))
    
    # I don't want the histogram "bars", just the count.
    # So I'm calling hist() for a plot I don't want to look at.
    axjunk     = junk.add_subplot(files, 1, i)
    axfig      = fig.add_subplot(files, 1, i)
    n, bins, patches = axjunk.hist(good_data, histtype='bar', bins=xbinedges)
    #print n
    
    # Plot the probability mass function for each channel's
    # (good) data, normalized so the pmf sums to one.
    for j, count in enumerate(n):
        # This commented line caused errors in past versions of matplotlib
        # because "count" was integers and so integer division would give zero.
        # axfig.plot(xbincenters, count/sum(count), markerlist[j])
        total = float(count.sum())
        axfig.plot(xbincenters, count/total, markerlist[j])

# Additional formatting done on each subplot
for i in range(files):
    ax = plt.subplot(files, 1, i)
    ax.set_xlim([minrss-1, maxrss+1])
    ax.set_ylabel('Probability Mass')
    ax.grid()
ax = plt.subplot(files, 1, 0)
ax.set_xlabel('RSS Value (dBm)')

# Wait for user to hit enter so that the plot doesn't disappear immediately
raw_input("Hit enter to end.")