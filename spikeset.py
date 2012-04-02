# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 04:18:30 2012

@author: Bernard
"""

import scipy as sp
import numpy as np
import scipy.io as sio
from features import *
import random
import matplotlib.nxutils as nx
import time

def loadNtt(filename):
    f = open(filename, 'rb')    
    header = f.read(2**14)

    # A tetrode spike record is as folows:
    # uint64 - timestamp                    bytes 0:8
    # uint32 - acquisition entity number    bytes 8:12
    # uint32 - classified cel number        bytes 12:16
    # 8 * uint32- params                    bytes 16:48
    # 32 * 4 * int16 - waveform points
    # hence total record size is 2432 bits, 304 bytes
    
    # header is supposed to be 16kbyte, i.e. 16 * 2^10 = 2^14
    
    # Read the header and findthe conversion factors
    a2d_conversion = [1,1,1,1]
    for line in header.split('\n'):
        if line.strip().startswith('-ADBitVolts'):
            a2d_conversion = 1e6 * np.array(map(float, line.split()[1:5]))
            break
    
    f.seek(2**14)    # start of the spike, records

#    dt = np.dtype([('time', np.uint64, (1)), ('acqe', np.uint32), ('cellno', np.uint32), ('params', np.uint32, (8,1)), ('spikes', np.dtype('>h'), (32,4))])

    dt = np.dtype([('time', np.uint64), ('filer', np.uint8, 2*(2+2+8*2)-1), ('spikes', np.dtype('>h'), (32,4)), ('filler2', np.uint8, 1)])
#    dt = np.dtype([('filler3', np.uint8, 1), ('time', np.uint64), ('filer', np.uint8, 2*(2+2+8*2)-2), ('spikes', np.dtype('>h'), (32,4)),('filler2', np.uint8, 1)])
        
    temp = np.fromfile(f, dt)
    
    # it somehow gets permuted, dont ask me how
    temp['spikes'] = temp['spikes'][:, :, [1,2,3,0]]
    
    #np.permute

    #f.close()  
#    f.seek(0)
#    contents = f.read()
#    
#    # Read the spikes    
#    cf = lambda s,p: [ s[i:i+p] for i in range(0,len(s),p) ]
#    records = cf(contents[2**14:], 304)
#    N = len(records)
#    
#    spikes = np.zeros([N, 4, 32])
#    timestamps= np.zeros([N,1])
#    
#    for i in range(0,N):
#        record = records[i]
#        timestamps[i], = struct.unpack('Q', record[0:8])    
#        for sample in range(1,33):
#            spikes[i-1,:,sample-1] = struct.unpack('h'*4, record[48 + 2 * 4 * (sample-1): 48 + 2*4*sample])
#    
#    spikes = spikes * np.reshape(a2d_conversion, [1,4,1])
#    return Spikeset(spikes, timestamps)
    return Spikeset(temp['spikes'] * np.reshape(a2d_conversion, [1,1,4]), temp['time'])

# rough breakdown as follows
# spikeset contains spikes, timestamps for the whole ntt file
# feature contains a name, channel count and data
# boundary contains a feature name, and polygon boundary
# cluster contains a list of boundaries, along with some summary statistics

# convention: N number of spikes, C number of channels, L length of waveforms

# Spike data is N x L x C
class Spikeset:
    def __init__(self, spikes, timestamps):
        self.spikes = spikes
        self.time = timestamps
        self.N = len(timestamps)
        self.features = [Feature_Peak(self)]
        self.T = (max(self.time) - min(self.time))/1e6
        
    def apply_bounds(self, bounds):
        pass

# Clusters have a color, a set of boundaries and some calculation functions
class Cluster:
    def __init__(self, spikeset):
        self.color = (random.randrange(0, 200), random.randrange(0,200), random.randrange(0,200))
        self.member = np.array([False] * spikeset.N)
        self.bounds = []
        self.isi = []
        self.refractory = np.array([False] * spikeset.N)
        self.stats = {}
        
    def __del__(self):
        print "Cluster object being destroyed"
        
    def addBound(self, bound):
        # bounds should be a 5 or 6 length tuple, 
        # (feature_name_x, feature_chan_x, feature_name_Y, feature_chan_y, polygon_points, extra_data)
        self.removeBound(bound[0], bound[1], bound[2], bound[3])
        self.bounds.append(bound)
        
    def getBoundPolygon(self, feature_name_x, feature_chan_x, feature_name_y, feature_chan_y):        
        f = lambda bound: (bound[0] == feature_name_x) and (bound[1] == feature_chan_x) and (bound[2] == feature_name_y) and (bound[3] == feature_chan_y)
        temp = [bound[4] for bound in self.bounds if f(bound)]
        if temp: 
            return temp[0]
        return None
        
    def removeBound(self, feature_name_x, feature_chan_x, feature_name_y, feature_chan_y):
        f = lambda bound: (bound[0] == feature_name_x) and (bound[1] == feature_chan_x) and (bound[2] == feature_name_y) and (bound[3] == feature_chan_y)        
        self.bounds = [bound for bound in self.bounds if not f(bound)]
        
    def calculateMembership(self, spikeset):
        if not self.bounds:
            self.member = np.array([False] * spikeset.N)
            self.isi = []
            self.refractory = np.array([False] * spikeset.N)
        else:
            self.member = np.array([True] * spikeset.N) # start with everything and cut it down
            for bound in self.bounds:            
                px = [feature.data[:, bound[1]] for feature in spikeset.features if feature.name == bound[0]][0]
                py = [feature.data[:, bound[3]] for feature in spikeset.features if feature.name == bound[2]][0]
    
                data = np.column_stack((px,py))
           
                self.member = np.logical_and(self.member, nx.points_inside_poly(data, bound[4]))
                
            t = spikeset.time[self.member]
            self.refr_period = 2
            self.burst_period = 20
            self.isi = (t[1:] - t[0:-1])/1e3
            
            self.refractory = np.array([False] * spikeset.N)
            ref = np.logical_and(self.isi < self.refr_period, self.isi> 0.8)
            self.refractory[self.member] = np.logical_or(np.append(ref, False), np.append(False, ref))
            
            #if np.sum(self.refractory) == 2:
            #    p = Feature_Peak(spikeset).data
            #    print p[self.refractory,:]
            #    print self.isi[self.isi < self.refr_period]
            #print 'spikes flagged as refractory', np.sum(self.refractory)
            
            #self.refractory[self.member] = np.logical_and(self.isi < self.refr_period, self.isi > 0.9)
            #self.refractory[self.member] = self.isi < self.refr_period
            
            # stats            
            self.stats['num_spikes'] = np.sum(self.member)
            self.stats['mean_rate'] = self.stats['num_spikes'] / spikeset.T            
            
            if self.stats['num_spikes'] <= 1:
                self.stats['burst'] = np.NAN
                self.stats['refr_count'] = np.NAN
                self.stats['csi'] = np.NAN
                self.stats['refr_fp'] = np.NAN
                self.stats['isolation'] = np.NAN
                self.stats['refr_frac'] = np.NAN
            else:
                self.stats['burst'] = 100.0 *  np.sum(self.isi < self.burst_period) / (self.stats['num_spikes'] - 1)
                self.stats['refr_count'] = np.sum(np.logical_and(self.isi < self.refr_period, self.isi> 0.8))
                self.stats['refr_frac'] = 100.0 * self.stats['refr_count'] / (self.stats['num_spikes'] - 1)
                alpha = self.stats['mean_rate'] * self.stats['refr_count'] / (self.stats['num_spikes'] * 2.0 * self.refr_period * 1e-3)
                if alpha > 0.25:
                    self.stats['refr_fp'] = 100.0
                else:
                    self.stats['refr_fp'] = 100.0 * 0.5 * (1 - np.sqrt(1.0 - 4.0 * alpha))
                            
                # csi and isolation needs peaks
                p = Feature_Peak(spikeset).data
                if self.stats['num_spikes']*2 > spikeset.N:
                    self.stats['isolation'] = np.NAN
                else:
                    cvi= np.linalg.inv(np.cov(np.transpose(p[self.member,:])))
                    u = p-np.mean(p[self.member,:], axis=0)
                    m = np.sum(np.dot(u, cvi) * u, axis=1)
                    m = np.sort(m)
                    self.stats['isolation'] = m[self.stats['num_spikes']*2-1]

                chan = np.round(np.mean(np.argmax(p[self.member,:], axis=1)))
                delta = np.diff(p[self.member, chan])
                
                delta= delta[np.logical_and(self.isi < self.burst_period, self.isi > self.refr_period)]
                if np.size(delta,0):
                    print "Neg", np.sum(delta <= 0), "Pos", np.sum(delta > 0), "Total", np.size(delta,0), "Delta", (np.sum(delta < 0) - np.sum(delta > 0))
                    #print delta
                    #print self.isi[self.isi < self.burst_period]
                    self.stats['csi'] = 100.0 * (np.sum(delta <= 0) - np.sum(delta > 0)) / np.size(delta,0)
                else:
                    self.stats['csi'] = np.NAN



if __name__ == "__main__":
    print "Loading ntt"
    spikeset = loadNtt('TT2_neo.ntt')    
    
    print "Constructing cluster object"
    clust = Cluster(spikeset)
    clust.calculateMembership(spikeset)
    poly = np.array([[50,0], [250,0], [50,150], [150,150]], float)
    bound = ('Peak', 0, 'Peak', 1, poly)
    clust.addBound(('Peak', 0, 'Peak', 1, poly))
    clust.addBound(('Peak', 0, 'Peak', 3, poly))
    clust.addBound(('Peak', 1, 'Peak', 2, poly))

    clust.calculateMembership(spikeset)
    
    print clust.stats
    
    clusters = [clust, clust]
    
    temp = np.column_stack([clust.member for clust in clusters])

    sio.savemat('test', {'cluster_id':temp, 'spike_time':spikeset.time}, oned_as='row')
    