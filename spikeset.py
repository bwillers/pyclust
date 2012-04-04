# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 04:18:30 2012

@author: Bernard
"""

import numpy as np
import random
import matplotlib.nxutils as nx
import time
import features

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
    # Neuralynx write little endian for some dumb reason
    dt = np.dtype([('time', '<Q'), ('filer', '<i', 10), ('spikes', np.dtype('<h'), (32,4))])
    temp = np.fromfile(f, dt)

    return Spikeset(temp['spikes'] * np.reshape(a2d_conversion, [1,1,4]), temp['time'], 8)
    
#    else:    
#        # Slow loopy way to read
#        f.seek(2**14)
#        contents = f.read()
#        
#        # Read the spikes    
#        cf = lambda s,p: [ s[i:i+p] for i in range(0,len(s),p) ]
#        records = cf(contents, 304)
#        N = len(records)
#        
#        spikes = np.zeros([N, 32, 4])
#        timestamps= np.zeros([N,1])
#        
#        for i in range(0,N):
#            record = records[i]
#            timestamps[i], = struct.unpack('Q', record[0:8])    
#            for sample in range(1,33):
#                spikes[i,sample-1,:] = struct.unpack('h'*4, record[48 + 2 * 4 * (sample-1): 48 + 2*4*sample])
#        
#        spikes = spikes * np.reshape(a2d_conversion, [1,1,4])
#        return Spikeset(spikes, timestamps)

# rough breakdown as follows
# spikeset contains spikes, timestamps for the whole ntt file
# feature contains a name, channel count and data
# boundary contains a feature name, and polygon boundary
# cluster contains a list of boundaries, along with some summary statistics

# convention: N number of spikes, C number of channels, L length of waveforms

# Spike data is N x L x C
class Spikeset:
    def __init__(self, spikes, timestamps, peak_index = 8):        
        self.spikes = spikes
        self.time = timestamps
        self.N = len(timestamps)
        self.peak_index = 8
        self.features = [features.Feature_Peak(self), features.Feature_Energy(self), features.Feature_Time(self), features.Feature_Valley(self), features.Feature_Aligned_Peak(self)]
        self.T = (max(self.time) - min(self.time))/1e6
        
    def __del__(self):
        print "Spikeset object being destroyed"
        
    def featureNames(self):
        return [feature.name for feature in self.features]
        
    def featureByName(self, name):
        for feature in self.features:
            if feature.name == name: return feature
        return None

# Clusters have a color, a set of boundaries and some calculation functions
class Cluster:
    def __init__(self, spikeset):
        self.color = (random.randrange(100, 200), random.randrange(100,200), random.randrange(100,200))
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
            #ref = np.logical_and(self.isi < self.refr_period, self.isi> 0.8)
            ref = self.isi < self.refr_period
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
                self.stats['refr_count'] = np.sum(self.isi < self.refr_period)
                self.stats['refr_frac'] = float(self.stats['refr_count'] )/ (self.stats['num_spikes'] - 1)
                alpha = self.stats['mean_rate'] * self.stats['refr_count'] / (self.stats['num_spikes'] * 2.0 * self.refr_period * 1e-3)
                if alpha > 0.25:
                    self.stats['refr_fp'] = 100.0
                else:
                    self.stats['refr_fp'] = 100.0 * 0.5 * (1 - np.sqrt(1.0 - 4.0 * alpha))
                            
                # csi and isolation needs peaks
                p = features.Feature_Peak(spikeset).data
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
                    self.stats['csi'] = 100.0 * (np.sum(delta <= 0) - np.sum(delta > 0)) / np.size(delta,0)
                else:
                    self.stats['csi'] = np.NAN



if __name__ == "__main__":

    print "Loading ntt"
    spikeset = loadNtt('Sample2.ntt')
    
    print spikeset.featureNames()
    
    print spikeset.features[0].valid_y_features(1)
    
    print spikeset.featureByName('Peak')
    print spikeset.featureByName('Peak').valid_x_features()
    
    print spikeset.featureByName('Time')
    print spikeset.featureByName('Time').valid_x_features()

    del spikeset