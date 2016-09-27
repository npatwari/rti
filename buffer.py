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
# Version 1.0:  Initial Release.  22 Sept 2016.
#

# ########################################
# Code to provide a fixed-length buffer data type
class MaxLenBuffer:
    
    # frontInd is the index of the most recently appended value.
    # backInd is the index of the oldest value.
    # len is the maximum length of the buffer, older data is overwritten.
    # data is initialized to some value, but that value is not ever used.
    # Without "empty" flag, frontInd=backInd would indicate one value is present in the list.
    # empty==True says that there is no data appended yet.
    def __init__(self, maxlen, initvalue=0):
        self.frontInd = -1   # Key to indicate empty buffer.
        self.backInd  = -1   # Key to indicate empty buffer.
        self.len      = maxlen
        self.data     = [initvalue]*maxlen
    
    
    # Return all of the buffer, oldest to newest.
    # If nothing has yet been appended, return [].
    # If buffer is partially full, just return what has been appended.
    # If more than self.len items have been added, return only the self.len
    #   newest items that remain in the buffer.
    def list(self):
        if (self.frontInd < self.backInd):
            val = self.data[self.backInd:] + self.data[:self.frontInd+1]
        elif self.backInd == -1:
            val = []
        else:
            val = self.data[self.backInd:self.frontInd+1]
        return val
    
    # If backInd == (frontInd+1 mod len) then the append also deletes the oldest item.
    # In that case also increment backInd.
    # Whenever incrementing, check for wrap-around.
    def append(self, newItem):
        
        # Append the new item
        self.frontInd = (self.frontInd + 1) % self.len
        self.data[self.frontInd] = newItem
        
        # If the buffer was empty it is not now. 
        if self.backInd == -1:
            self.backInd = 0 # Now point backInd to the oldest value, same as newest.
        # Otherwise, if the list was full, then the oldest was deleted, 
        # so increment backInd  
        elif self.frontInd == self.backInd:
            self.backInd = (self.backInd + 1) % self.len

    # Clear an existing buffer (without reallocating memory for the data list) 
    def empty(self):
        self.frontInd = -1   # Key to indicate empty buffer.
        self.backInd  = -1   # Key to indicate empty buffer.
        
    # Returns the "front" item, or [] if nothing has yet been appended.
    def mostRecent(self):
        if self.frontInd == -1:
            val = []
        else:
            val = self.data[self.frontInd]
        return val
    
    # How many values are stored?  Must be between 0 and len
    def numStored(self):
        if self.backInd == -1:
            val = 0
        else:
            val = ((self.frontInd-self.backInd) % self.len) + 1
        return val

    # Returns the N items most recently appended, oldest first.
    def mostRecentN(self,N):
        Neff = min(self.numStored(), N)  # "effective" N
        return [self.data[(self.frontInd-i)%self.len] for i in range(Neff-1,-1,-1)]
    
# ########################################
