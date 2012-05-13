# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 19:53:55 2012

@author: Bernard
"""

import numpy as np

# Feature data is N x C


class Feature:

    def __init__(self, name, spikeset):
        self.name = name
        self.data = self.calculate(spikeset)
        shape = np.shape(self.data)
        if len(shape) == 1:     # We need an (N,1) array rather than (N,)
            self.data = np.reshape(self.data, (self.data.size, 1))
        self.channels = np.size(self.data, 1)
        self.valid_y_all_chans = []
        self.valid_y_same_chan = []
        #self.channels = np.size(self.data,1)

    def calculate(self, spikeset):
        pass

    def valid_y_features(self, feature_chan):
        ret_val = []
        if feature_chan < self.channels - 1:
            ret_val.append((self.name, range(feature_chan + 1, self.channels)))
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
                return range(0, self.channels - 1)


#class Feature_Peak(Feature):
#    def __init__(self, spikeset):
#        Feature.__init__(self, 'Peak', spikeset)
#        self.valid_y_same_chan.append('Energy')
#        self.valid_y_same_chan.append('Valley')
#
#    def calculate(self, spikeset):
#        return np.max(spikeset.spikes, axis=1)


class Feature_Valley(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Valley', spikeset)

    def calculate(self, spikeset):
        return np.min(spikeset.spikes[:, spikeset.peak_index:, :], axis=1)


class Feature_Energy(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Energy', spikeset)
        self.valid_y_same_chan.append('Peak')

    def calculate(self, spikeset):
        return np.sum(np.power(spikeset.spikes, 2), axis=1)


class Feature_Time(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Time', spikeset)
        self.valid_y_all_chans.append('Peak')

    def calculate(self, spikeset):
        #temp = np.zeros([spikeset.N, 1])
        #temp[:, 0] =
        return np.array(spikeset.time - spikeset.time[0], dtype=np.float) / 1e6
        #return temp


class Feature_Peak(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Peak', spikeset)
        self.valid_y_same_chan.append('Energy')
        self.valid_y_same_chan.append('Valley')

    def calculate(self, spikeset):
        searchrange = range(spikeset.peak_index - 3, spikeset.peak_index + 4)
        return np.max(spikeset.spikes[:, searchrange, :], axis=1)


class Feature_Barycenter(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Barycenter', spikeset)

    def calculate(self, spikeset):
        p = Feature_Peak(spikeset)
        p = p.data

        p = p - np.min(np.max(p, axis=1))


        x = p[:,0] - p[:,2]
        y = p[:,1] - p[:,3]
        angle = np.arctan2(y,x)
        #amp = p - np.mean(p, axis=0)
        #amp[amp < 0] = 0
        #amp = np.sum(amp, axis=1)
        #amp = np.sum(p, axis=1)
        retval = np.zeros((np.size(p,0),2))
        retval[:,0] = np.cos(angle)# * amp
        w = retval[:,0] > 0
        retval[w,0] = retval[w,0] * p[w,0]
        w = np.logical_not(w)
        retval[w,0] = retval[w,0] * p[w,2]

        retval[:,1] = np.sin(angle)# * amp
        w = retval[:,1] > 0
        retval[w,1] = retval[w,1] * p[w,1]
        w = np.logical_not(w)
        retval[w,1] = retval[w,1] * p[w,3]

        return retval


class Feature_FallArea(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Fall Area/Time', spikeset)
#        self.valid_y_all_chans.append('Fall Time')

    def calculate(self, spikeset):
        wv = np.mean(spikeset.spikes, axis=2)  # Average over channels
        ret = np.sum( wv[:, spikeset.peak_index:], axis=1)
        ind = np.argmin(wv[:, spikeset.peak_index:], axis=1)
        rnd = np.random.rand(np.size(ind,0)) - 0.5
        retval = np.zeros((np.size(wv,0), 2))
        retval[:, 0] = ret
        retval[:, 1] = (ind + rnd) * 1000.0 / spikeset.fs
        print np.shape(retval)
        return retval

