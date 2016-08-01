#! /usr/bin/env python

# This script reads packet data from the listen node through the serial port
# or from a file, with file name given in the -i option on the command line,
# and prints one RSS measurement on each line for all tx, rx, ch combination,
# where the rx != tx, both rx and tx are in the list of sensors, and ch is in 
# the list of channels.  There is an option to beep at a given rate to provide
# feedback to the experimenter about what time it is / where they should be.
# This version reads in the text from the output of the XANDEM HOME kit 
# "gateway" program.
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
# Version 0.1:  Initial Release, 2 July 2016
# Version 0.2:  4 July 2016: outputs time in ms as an integer rather than a float.


import sys
import rss 
import getopt


# If you want beeps and/or second printing, set to True
soundOption     = False
printOption     = False
notMeasuredCode = 127
numNodes        = 10
numChs          = 4
beepRate        = 1.0  # Beeps per second
startSkip       = 24

# Get the file name and any other options from the command line
file_method = "stdin"
fin         = sys.stdin
# Accept command line inputs for file input
myopts, args = getopt.getopt(sys.argv[1:],"i:")
for o, a in myopts:
    if o == "-i":
        file_method = "file"
        fname =  a
        fin = open(a,'r')
    else:
        sys.stderr.write("Usage: %s -i filename.txt" % sys.argv[0])

if printOption:
    print fin

# Constants calculated from the inputs
nodeList        = range(1, numNodes+1)   # Xandem nodes are 1...10
numLinks        = numNodes*(numNodes-1)*numChs
channelList     = range(numChs)
numLines        = numNodes*numChs

# remove junk from start of file.
for i in range(startSkip): 
    line        = fin.readline()

# Use the first line to determine the starting time
lineStrings     = line.split(', ')
time_start      = float(lineStrings[-1])

# Run forever, adding lines from input as they are available.
lineCounter     = 0
beepCounter     = 0
currentLinkRSS  = [127] * numLinks
while 1:
    #  RSS is in columns 11 through 20
    rss_now          = [int(x) for x in lineStrings[11:21]]
    ch_now           = int(lineStrings[6])
    rxid_now         = int(lineStrings[7])
    time_now         = float(lineStrings[-1])
    time_diff        = time_now - time_start
    time_diff_ms     = int(1000.0*time_diff)
    if printOption:
        print "rss_now: " + str(rss_now)
        print "ch_now: " + str(ch_now)
        print "rxid_now: " + str(rxid_now)
        print "time_now: " + str(time_now)

    # If this row is the first row of a new set of RSS data, 
    # then output the old rss data, init a new one 
    # (and enter the new data into the new rss array)
    if (lineCounter == numLines):
        sys.stdout.write(' '.join(map(str,currentLinkRSS)) + ' ' + str(time_diff_ms) + '\n')
        sys.stdout.flush()
        # Restart with a new line by resetting currentLinkRSS
        currentLinkRSS   = [127] * numLinks
        lineCounter      = 0

    # Beeping every beepRate
    if soundOption:
        beepCounter  = rss.beeping(time_diff, beepRate, beepCounter, printOption)

    # Add each RSS to the currentLinkRSS, at the index for the linkId
    for i,rssa in enumerate(rss_now):
        txid         = i+1
        # Link between tx=i to rx=i does not exist
        if txid != rxid_now:   
            linkId   = rss.linkNumForTxRxChLists(txid, rxid_now, ch_now, nodeList, channelList)
            currentLinkRSS[linkId] = rssa
            if printOption:
                print "txid="+str(txid)+", link="+str(linkId)+", rss="+str(rssa)

    # Read a new line, quit if at end of file, or wait for more from stdin
    line = fin.readline()
    if (not line) & (file_method == "file"):
        break
    if (not line) & (file_method == "stdin"):
        continue
    while line[-1] != '\n':   # If line is incomplete, add to it.
        line        += fin.readline()
        if (not line) & (file_method == "file"):
            break

    # Split the line at each comma
    lineStrings      = line.split(', ')
    # Increment the line counter
    lineCounter     += 1


