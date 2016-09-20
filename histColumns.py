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
# Version 2.1: 19 Sept 2016: Update to allow a single link to be plotted.
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
links = len(linkNumList)

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
    
    # Print error message if a link has NO good data
    for j in range(links):
       if good_data[j] == []: 
           sys.exit("Error: Link " + str(linkCombosToPlot[j]) + " in file " + file + " has no non-127 values.")
    
    # Update the min and max
    maxrss = max(maxrss, max([gd.max() for gd in good_data]))
    minrss = min(minrss, min([gd.min() for gd in good_data]))
    
    # I don't want the histogram "bars", just the count.
    # So I'm calling hist() for a plot I don't want to look at.
    axjunk     = junk.add_subplot(files, 1, i+1)
    axfig      = fig.add_subplot(files, 1, i+1)
    n, bins, patches = axjunk.hist(good_data, histtype='bar', bins=xbinedges)
    
    # Plot the probability mass function for each channel's
    # (good) data, normalized so the pmf sums to one.
    if links>1:
        for j, count in enumerate(n):
            # This commented line caused errors in past versions of matplotlib
            # because "count" was integers and so integer division would give zero.
            # axfig.plot(xbincenters, count/sum(count), markerlist[j])
            total = float(count.sum())
            axfig.plot(xbincenters, count/total, markerlist[j])
    # if one link, the n variable is only a single list, not a list of lists,
    # so you can't use enumerate(n) and get the same result.
    else:
        total = float(n.sum())
        axfig.plot(xbincenters, n/total, markerlist[0])        

# Additional formatting done on each subplot
for i in range(files):
    ax = plt.subplot(files, 1, i+1)
    ax.set_xlim([minrss-1, maxrss+1])
    ax.set_ylabel('Probability Mass')
    ax.grid()
ax = plt.subplot(files, 1, 1)
ax.set_xlabel('RSS Value (dBm)')

# Wait for user to hit enter so that the plot doesn't disappear immediately
raw_input("Hit enter to end.")