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
# Version 1.0:  Initial Release.  27 Oct 2014.
#
# Version 1.1: 26 Jan 2016
#   a) Cleaned up & added comments, 
#   b) Replaced print with sys.stderr.write
#
# Version 1.2: 3 July 2016
#   a) Added floor_multiple_of() and ceil_multiple_of()
#   b) Added beeping()
#

import sys
import platform
from math import floor, ceil 

# Convert Tx, Rx, and Channel numbers to link number
# Possible transmitter and receiver numbers are listed in nodeList
# Possible channels are listed in channelList
def linkNumForTxRxChLists(tx, rx, ch, nodeList, channelList):
    if (nodeList.count(tx) == 0) or (nodeList.count(rx) == 0) or (channelList.count(ch) == 0):
        sys.stderr.write('Error in linkNumForTxRx: tx, rx, or ch number invalid')
    rx_enum = nodeList.index(rx)
    tx_enum = nodeList.index(tx)
    ch_enum = channelList.index(ch)
    nodes = len(nodeList)
    links = nodes*(nodes-1)
    linknum = ch_enum*links + tx_enum*(nodes-1) + rx_enum
    if (rx_enum > tx_enum):
        linknum -= 1
    return linknum

# Convert link number to Tx, Rx, and channel numbers
# Possible transmitter and receiver numbers are listed in nodeList
# Possible channels are listed in channelList
def txRxChForLinkNum(linknum, nodeList, channelList):
    nodes   = len(nodeList)
    # Transmitter / receiver combinations.  A node can't transmit to itself.
    links   = nodes*(nodes-1)  
    ch_enum = linknum / links
    remLN   = linknum % links
    tx_enum = remLN / (nodes-1)
    rx_enum = remLN % (nodes-1)
    if (rx_enum >= tx_enum):
        rx_enum+=1
    if (tx_enum >= nodes) | (ch_enum > len(channelList)):
        sys.stderr.write('Error in txRxForLinkNum: linknum or ch too high for nodes, channels values')
    else:
        ch = channelList[ch_enum]
        tx = nodeList[tx_enum]
        rx = nodeList[rx_enum]
    return (tx, rx, ch)



def hex2signedint(he):
    # Convert from hexidecimal 2's complement to signed 8 bit integer
    return (int(he,16) + 2**7) % 2**8 - 2**7

# We are currently receiving on ch_now.  What was the previous channel txed on?
def prevChannel(channelList, ch_now):
    if (channelList.count(ch_now) > 0):
        i = channelList.index(ch_now)
        rval = channelList[(i-1) % len(channelList)]
    else:
        rval = -1  # Key for bad ch_now input
    return rval

# Return the file name associated with the USB listen node.
#
# This requires some fiddling to get it to automatically return the right file
#       name, depending on the OS and how that OS works.
#
# USER: The following serial "file name" changes depending on your operating
#       system, and what name is assigned to the serial port when your listen
#       node is plugged in.
def serialFileName():    
    system_name = platform.system()
    #
    # LINUX USERS
    if system_name == 'Linux':
        # Automatically grab the USB filename (since the number after /dev/ttyACM may vary)
        usb_file_list = glob.glob('/dev/ttyACM*')
        if len(usb_file_list) > 0:
            serial_filename =  usb_file_list[0]  
        else:
            sys.stderr.write('Error: No Listen node plugged in?\n')
            serial_filename = '0'
    #
    # WINDOWS USERS: Change 'COM#' to match what the system calls your USB port.
    elif system_name == 'Windows':
        serial_filename = 'COM3'
    #
    # MAC USERS
    else:  # 'Darwin' indicates MAC OS X
        # Automatically grab the USB filename (since the number after /dev/tty.usb may vary)
        usb_file_list = glob.glob('/dev/tty.usb*')
        if len(usb_file_list) > 0:
            serial_filename =  usb_file_list[0]  
        else:
            sys.stderr.write('Error: No Listen node plugged in?\n')
    #
    return serial_filename
    
# provide a beep at rate beepRate to provide feedback on what time the computer
# thinks it is, and a printed counter value if printOption==True.
def beeping(time_diff, beepRate, beepCounter, printOption):
    
    curBeepNumber = int(time_diff/beepRate)

    if (curBeepNumber > beepCounter):
        beepCounter = curBeepNumber
        sys.stderr.write('\a')  # BEEP!
        if printOption:
            sys.stderr.write(str(curBeepNumber/4.0 ) + '\n')
        if curBeepNumber % 4 == 0:
            sys.stderr.write('\a')  # Double beep each "measure"
            if printOption:
                sys.stderr.write('---\n')
    return beepCounter

# Finds the next lower multiple of "multiple_of" for "value"
# Useful in setting axes in real time plotting scripts.  
def floor_multiple_of(multiple_of, value):
    return floor(float(value)/multiple_of)*multiple_of

# Finds the next higher multiple of "multiple_of" for "value"
# Useful in setting axes in real time plotting scripts.  
def ceil_multiple_of(multiple_of, value):
    return ceil(float(value)/multiple_of)*multiple_of


