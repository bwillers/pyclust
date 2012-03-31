# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 04:18:30 2012

@author: Bernard
"""

import numpy as np
import scipy as sp
import struct
from features import *
import random

def loadNtt(filename):
    f = open(filename, 'rb')    
    contents = f.read()   
    assert np.mod(len(contents) - 2**14, 304) == 0     
    # A tetrode spike record is as folows:
    # uint64 - timestamp                    bytes 0:8
    # uint32 - acquisition entity number    bytes 8:12
    # uint32 - classified cel number        bytes 12:16
    # 8 * uint32- params                    bytes 16:48
    # 32 * 4 * int16 - waveform points
    # hence total record size is 2432 bits, 304 bytes
    
    # header is supposed to be 16kbyte, i.e. 16 * 2^10 = 2^14
    
    # Read the header and findthe conversion factors
    header = contents[0:2**14-1]
    a2d_conversion = [1,1,1,1]
    for line in header.split('\n'):
        if line.strip().startswith('-ADBitVolts'):
            a2d_conversion = 1e6 * np.array(map(float, line.split()[1:5]))
            break
    
    # Read the spikes    
    cf = lambda s,p: [ s[i:i+p] for i in range(0,len(s),p) ]
    records = cf(contents[2**14:], 304)
    N = len(records)
    
    spikes = np.zeros([N, 4, 32])
    timestamps= np.zeros([N,1])
    
    for i in range(0,N):
        record = records[i]
        timestamps[i], = struct.unpack('Q', record[0:8])    
        for sample in range(1,33):
            spikes[i-1,:,sample-1] = struct.unpack('h'*4, record[48 + 2 * 4 * (sample-1): 48 + 2*4*sample])
    
    spikes = spikes * np.reshape(a2d_conversion, [1,4,1])
    return Spikeset(spikes, timestamps)

# rough breakdown as follows
# spikeset contains spikes, timestamps for the whole ntt file
# feature contains a name, channel count and data
# boundary contains a feature name, and polygon boundary
# cluster contains a list of boundaries, along with some summary statistics

# convention: N number of spikes, C number of channels, L length of waveforms

# Spike data is N x C x L
class Spikeset:
    def __init__(self, spikes, timestamps):
        self.spikes = spikes
        self.time = timestamps
        self.N = len(timestamps)
        
    def apply_bounds(self, bounds):
        pass

# Clusters have a color, a set of boundaries and some calculation functions
class Cluster:
    def __init__(self, spikeset):
        self.color = (random.randrange(0, 200), random.randrange(0,200), random.randrange(0,200))
        N = np.size(spikeset.spikes,0)
        temp = np.arange(1, N+1)        
        np.random.shuffle(temp)
        self.member = temp < N/10
        
    def __del__(self):
        print "Cluster object being destroyed"
        
    def calculateMembership(self, spikeset):
        pass
        
        #perm = np.random.shuffle(temp)
        #print perm
        #member = perm < N / 2
        #print member
        #print sum(member)


if __name__ == "__main__":
    temp = (255, 180, 30)
    print map(lambda s: s / 255.0, temp)
    print "Loading ntt"
    spikeset = loadNtt('Sample1.ntt')
    print "Creating feature object"
    peaks = Feature_Peak(spikeset)
    print "Constructing cluster object"
    clust = Cluster(spikeset)
    #print clust.member.type
    print type(clust.member)
    #clust.calculateMembership(spikeset)
    w = [True] * spikeset.N
    blah = np.array(w)
    print blah
    