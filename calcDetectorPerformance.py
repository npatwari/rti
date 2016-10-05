#! /usr/bin/env python

#
# Author: Neal Patwari
#
# License: GPL v3
#
# PURPOSE:
#   Read in a crossings estimate file, the correct crossing times, and the time stamps from the listen_output.txt  
#   and calculate its false alarm rate, detection rate, and "score".
#
# Usage: calcDetectorPerformance.py -c crossingsEstFile.txt -t trueCrossingsFile.txt -r listen_output.txt
#
# Version History:
#
# Version 1.0:  30 Sept 2016: Initial Release
#


import numpy as np
import sys
# Agg is a backend that allows the figure to be saved to a file without displaying it to the screen.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os.path
import getopt

# Increase font size in the figure
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

#  Example call:
#  python calcDetectorPerformance
#
# Accept command line inputs for:
#   crossingDetections file: The i'th row is the output of the detector at the time 
#      that the i'th row of RSS data is collected.  "1" for crossing detected,
#      "0" for no crossing detected.
#   rssDataFile: The output of listenx.py, one row of RSS for each time.
#   trueCrossings file: A list of the true crossing times, in seconds.
myopts, args = getopt.getopt(sys.argv[1:],"c:r:t:")
for o, a in myopts:
    if o == "-c":
        stateEstFile = a
        stateEst    = np.loadtxt(stateEstFile)
        rowsEstFile = len(stateEst)
        print('crossingDetections File ' + stateEstFile + ' contains ' + str(rowsEstFile) + ' records.')
    elif o == "-r":
        # Find the time stamps from the rss data file (the output of listenx.py)
        timeStamps = np.loadtxt(a)[:,-1]
        rowsRSSFile = len(timeStamps)
        print('RSS Data File ' + a + ' contains ' + str(rowsRSSFile) + ' records.')
    # Ground truth file is assumed to have crossing times listed in seconds
    elif o == "-t":
        trueCrossings = np.loadtxt(a)*1000.0  # convert to ms
        numTrueCrossings  = len(trueCrossings)
        print('trueCrossings File ' + a + ' contains ' + str(numTrueCrossings) + ' true crossing times.')
        print(trueCrossings)


# Quit, if we can't determine what time each crossingDetection corresponds to, 
#   or if there are no true crossing times.
if rowsEstFile != rowsRSSFile:
    sys.exit("There must be an identical number of rows in the crossingDetections file as in the RSS data file")
if numTrueCrossings == 0:
    sys.exit("There must be at least one true crossing time")

# milliseconds on each side of the "true" crossing time that a crossing detection is acceptable.
delta             = 1500 

# Process the crossing Time Estimates from the file
detectedInd = np.where(stateEst)[0]
detectedTimes = [timeStamps[i] for i in detectedInd]

# Assume: 
# - The actual crossing times are in order
# - There are no crossings within delta of time 0 or the final time.
pastEndTime = 0
correctDetections = 0
falseAlarms = 0
for tc in trueCrossings:
    
    # Find any detections within +/- delta around an actual crossing time
    # If there are any such detections, it is a correct detection.
    currentCrossInd = np.where(np.abs(detectedTimes-tc) < delta)[0]
    if len(currentCrossInd) > 0:
        correctDetections += 1

    # There shouldn't be detections before the (current crossing time - delta)
    # and after the (past crossing time + delta).  Count such detections as false alarms.
    currentNoCrossInd = np.where((detectedTimes >= pastEndTime) & (detectedTimes < tc - delta))[0]
    falseAlarms += len(currentNoCrossInd)

    # Save the previous tc + delta for use in the next iteration.
    pastEndTime = tc + delta

# Finally, count the false alarms in the period after the last true crossing time
falseAlarms += len(np.where(detectedTimes >= pastEndTime)[0])
    
# Compute the correct detection rate and false alarm rate.  
# Ensure floating point division by adding 0.0.
correctDetectionRate = (correctDetections + 0.0) / len(trueCrossings)
falseAlarmRate       = (falseAlarms + 0.0) / rowsEstFile

# Print to screen the detection and false alarm rates, and score.
print('correctDetectionRate = {:9.7f}'.format(correctDetectionRate))
print('falseAlarmRate = {:9.7f}'.format(falseAlarmRate))
print('score = correctDetectionRate - falseAlarmRate = {:9.7f}'.format(correctDetectionRate - falseAlarmRate))


###########################################################
# PLOT THE RESULTS

# Plot the output, the correct times
plt.figure()
plt.plot(timeStamps/1000., stateEst, '-x')
plt.ylabel('State Est.', fontsize=16)
plt.xlabel('Time (s)', fontsize=16)
plt.xlim(0, timeStamps[-1]/1000.)
plt.ylim(-0.1, 1.1)

# Mark where the training period ends
for tc in trueCrossings:
    
    # Show the +/- delta window around each true crossing time (TCT)
    plt.plot([(tc-delta)/1000., (tc-delta)/1000.], [-0.1, 1.1], 'c-')
    plt.plot([(tc+delta)/1000., (tc+delta)/1000.], [-0.1, 1.1], 'c-')
    # Underline in green a correct detection (Yeah!)
    if np.any(stateEst[np.where((timeStamps >= tc-delta) & (timeStamps <= tc+delta))]):
        plt.plot([(tc-delta)/1000., (tc+delta)/1000.], [-0.05, -0.05], 'g-', linewidth=2)

# Put an X below detections which were false alarms.
for i in detectedInd:
    detectedTime = timeStamps[i]
    if len(np.where(np.abs(trueCrossings - detectedTime) < delta)[0]) == 0:
        plt.plot([detectedTime/1000.], [-0.05], 'rx', linewidth=2)

# Auto save the file with same base file name as the state estimate file.
(basename, ext) = os.path.splitext(stateEstFile)
outfile = basename + '.png'
print("Outfile: " + outfile)
plt.savefig(outfile, bbox_inches='tight')
