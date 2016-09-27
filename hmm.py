#! /usr/bin/env python

#
# Author: Neal Patwari
#
# License: GPL v3
#
# Version History:
#
# Version 1.1:  19 Sept 2014: Initial Release
# Version 1.2:  22 Sept 2014: Added Viterbi algorithm
#
# Purpose: Allow forward, backward, and viterbi estimation algorithms
#     to be run for a finite state, finite observation hidden
#     Markov model.  States must be 0, ..., (S-1) for S states.
#     Observations can take any values, but the possible values must be
#     listed in "V".
#
# See: 

import sys
import numpy as np


class nealsHMM:

    # HMM Settings (from Rabiner 1989 tutorial)
    # A: Probability of transition in one step.  Matrix with A[i,j] as the
    #    one-step prob of transition from i to j.
    # B: Observation conditional probability matrix.  Matrix with B[i,j]
    #    as the prob of jth possible observation given current state i.
    # pi: initial state probability.  Vector pi[i] is the probability of
    #    being in state i at time 0.
    # V: All possible observation values (a python list), eg., the alphabet
    #
    # Internal variables:
    # data: list stores the index (zero-starting) of each observation
    # alpha: The forward state probabilities for each time index
    # beta: The backward state probabilities for each time index
    #
    def __init__(self, B, A, pi, V):
        self.B = B
        self.A = A
        self.pi = pi
        self.V = list(V)  # no math is done on V, so it is nice to have a list.
        self.data = []
        self.alpha = []
        self.beta = []
        self.delta = []
        self.phi = []
        self.qstar = []
        self.states = A.shape[0]
    
    # Append new observations in newdata list to the data list.
    def observe(self, newdata):
        for i in newdata:
            # use the index of the value i in V
            # Note that new data NOT in V will cause an error.
            # Make sure to include every possible observation in V.
            s = self.V.index(i)
            self.data.append(s)

    # Compute the forward joint probability alpha.  Compute starting after
    # the previous stoping point in forward algorithm computation.
    def forward(self):
        k = len(self.alpha)
        T = len(self.data)
        if k==0 and k<T:
            temp = self.B[:,self.data[0]] * self.pi
            self.alpha.append(temp/sum(temp))
            k += 1
        
        for t in range(k,T):
            temp = self.B[:,self.data[t]] * np.dot(self.alpha[-1], self.A)
            self.alpha.append(temp/sum(temp))

        return self.alpha

    # Solve the backward algorithm to find beta.  Do for entire data, unless
    # "steps" is specified.
    def backward(self, steps=None):
        T  = len(self.data)
        if steps==None:
            steps = T
        if steps > T:
            sys.exit("Error: backward(steps) being called with steps > number of data points")
        self.beta = [[] for i in range(T)]
        if T > 0:
            self.beta[-1] = [1]*self.states

        for t in range(T-2,-1,-1):
            temp = np.dot(self.A, self.B[:,self.data[t+1]] * self.beta[t+1])
            self.beta[t] = (temp/sum(temp))

        return self.beta

    # Calculate the probability of each state at "steps-1" time steps ago
    #   Specifically, steps==1 is the current time, steps==2 is one time step ago.
    def forwardBackward(self, steps=None):
        self.forward()
        if steps==None:
            steps = len(self.data)
        self.backward(steps)
        temp = self.beta[-steps] * self.alpha[-steps]
        return temp / sum(temp)

    # Find the most likely sequence of states from time T-steps to the current time.
    # Notes:
    # 1. This function assumes that you have added new data using the "observe"
    #    command first. Then you call viterbi to update the result for the new
    #    data points.  The "delta" and "phi" for the old data do not need to be
    #    updated in this new call.  However, "qstar" will be updated (for all time
    #    if steps == T is sent in, or for "steps" recent data points).
    # 2. Use "steps" (steps should be in {1, 2, ..., T} to indicate how far back
    #    in history you want the optimal state "qstar" to be calculated.
    def viterbi(self, steps=None):
        T = len(self.data)
        k = len(self.delta)
        self.qstar = [-1]*T # for i in range(T)]
        if steps == None:
            steps = T
        if steps > T:
            sys.exit("Error: viterbi(steps) being called with steps > number of data points")
        
        
        # Equations 32a and 32b in Rabiner "Initialization"
        if k==0 and k<T:
            self.delta.append( self.B[:,self.data[0]] * self.pi )
            self.phi.append( [0]*self.states )
            k += 1
        
        # Equations 33a and 33b in Rabiner "Recursion"
        for t in range(k,T):
            tempdelta    = [0.0]*self.states
            tempphi      = [0]*self.states
            for j, a_all_j in enumerate(self.A.T):
                prod_a_d     = a_all_j * self.delta[t-1]
                tempdelta[j] = prod_a_d.max() * self.B[j,self.data[t]]
                tempphi[j]   = prod_a_d.argmax()
            self.delta.append( tempdelta / sum(tempdelta) )
            self.phi.append( tempphi )

        # Equations 34a and 34b in Rabiner "Termination"
        # Pstar      = max(self.delta[-1])
        self.qstar[-1] = self.delta[-1].argmax()
        
        # Equation 35 in Rabiner "Path backtracking"
        for t in range(T-2,T-steps-1,-1):
            self.qstar[t] = self.phi[t+1][self.qstar[t+1]]

        return self.qstar


