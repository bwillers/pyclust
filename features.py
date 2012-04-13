# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 19:53:55 2012

@author: Bernard
"""

import numpy as np
    
# Feature data is N x C
class Feature:
    def __init__(self, name, spikeset):
        self.name = name;
        self.data = self.calculate(spikeset)
        #shape = np.shape(self.data)
        #self.channels = 1 if len(shape) == 1 else shape[1]
        self.channels = np.size(self.data,1)
        self.valid_y_all_chans = []
        self.valid_y_same_chan = []
        #self.channels = np.size(self.data,1)
        
    def calculate(self, spikeset):
        pass
    
    def valid_y_features(self, feature_chan):
        ret_val = []
        if feature_chan < self.channels - 1:
            ret_val.append((self.name, range(feature_chan+1, self.channels)))
        for item in self.valid_y_all_chans:
            ret_val.append((item, None))
        for item in self.valid_y_same_chan:
            ret_val.append((item, [feature_chan]))
        return ret_val
        
    def valid_x_features(self):
        if self.channels == 1:
            return [0]
        else:
            if self.valid_y_all_chans or self.valid_y_same_chan:
                return range(0, self.channels)
            else:
                return range(0, self.channels-1)
    
class Feature_Peak(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Peak', spikeset)
        self.valid_y_same_chan.append('Energy')
        self.valid_y_same_chan.append('Valley')
        
    def calculate(self, spikeset):
        return np.max(spikeset.spikes, axis=1)
        
class Feature_Valley(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Valley', spikeset)
        
    def calculate(self, spikeset):
        return np.min(spikeset.spikes[:, spikeset.peak_index:,:], axis=1)

class Feature_Energy(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Energy', spikeset)
        self.valid_y_same_chan.append('Peak')
        
    def calculate(self, spikeset):
        return np.sum(np.power(spikeset.spikes,2), axis=1)
        
class Feature_Time(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Time', spikeset)
        self.valid_y_all_chans.append('Peak')
        
    def calculate(self, spikeset):
        temp = np.zeros([spikeset.N,1])
        temp[:,0] = np.array(spikeset.time - spikeset.time[0], dtype=np.float) / 1e6
        return temp

class Feature_Aligned_Peak(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Peak 6-11', spikeset)
        self.valid_y_same_chan.append('Peak')
        self.valid_y_same_chan.append('Energy')
        self.valid_y_same_chan.append('Valley')
        
    def calculate(self, spikeset):
        return np.max(spikeset.spikes[:, range(spikeset.peak_index-3, spikeset.peak_index+4),:], axis=1)        