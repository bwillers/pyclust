# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 04:18:30 2012

@author: Bernard
"""

import scipy.special as spspec
import numpy as np
import random
import matplotlib.nxutils as nx
import features
import struct


def loadNtt(filename):
    f = open(filename, 'rb')
    header = f.read(2 ** 14)

    # A tetrode spike record is as folows:
    # uint64 - timestamp                    bytes 0:8
    # uint32 - acquisition entity number    bytes 8:12
    # uint32 - classified cel number        bytes 12:16
    # 8 * uint32- params                    bytes 16:48
    # 32 * 4 * int16 - waveform points
    # hence total record size is 2432 bits, 304 bytes

    # header is supposed to be 16kbyte, i.e. 16 * 2^10 = 2^14

    # Read the header and findthe conversion factors
    a2d_conversion = [1, 1, 1, 1]
    for line in header.split('\n'):
        if line.strip().startswith('-ADBitVolts'):
            a2d_conversion = 1e6 * np.array(map(float, line.split()[1:5]))
            break
    f.seek(2 ** 14)    # start of the spike, records
    # Neuralynx write little endian for some dumb reason
    dt = np.dtype([('time', '<Q'), ('filer', '<i', 10),
        ('spikes', np.dtype('<h'), (32, 4))])
    temp = np.fromfile(f, dt)

    return Spikeset(temp['spikes'] * np.reshape(a2d_conversion, [1, 1, 4]),
            temp['time'], 8, 32556)

def readStringFromBinary(f):
    strlen, = struct.unpack('<I', f.read(4))
    if strlen:
        return f.read(strlen)
    else:
        return ''

def loadDotSpike(filename):
    f = open(filename, 'rb')
    # Everything is little endian
    # A .spike record is as folows:
    # Header:
    # uint16 - version no
    # uint64 - num spikes
    # uint16 - num channels
    # uint16 - num samples per waveform
    # uint32 - sampling frequency
    # uint16 - peak align point
    # 4 * float64 - a2d conversion factor
    # uint32 + n x char - date time string
    # uint32 + n x char - subject string
    # uint32 + n x char - filter description

    version_no, = struct.unpack('<H', f.read(2))
    if version_no != 1:
        f.close()
        return

    print 'Loading', filename
    print 'Format version #', version_no
    num_spikes, = struct.unpack('<Q', f.read(8))
    print 'Num spikes', num_spikes
    num_chans, = struct.unpack('<H', f.read(2))
    print 'Num channels', num_chans
    num_samps, = struct.unpack('<H', f.read(2))
    print 'Num samples', num_samps
    fs, = struct.unpack('<I', f.read(4))
    print 'Sampling frequency', fs
    peak_align, = struct.unpack('<H', f.read(2))
    print 'Peak alignment point', peak_align
    uvolt_conversion = np.array(struct.unpack('<dddd', f.read(8 * 4)))
    print 'Microvolt conversion factor', uvolt_conversion
    datestr = readStringFromBinary(f)
    print 'Date string', datestr
    subjectstr = readStringFromBinary(f)
    print 'Subject string', subjectstr
    filterstr = readStringFromBinary(f)
    print 'Filter string', filterstr

    # Records:
    # uint64 - timestamp                    bytes 0:8
    # numsample x numchannel x int16 - waveform points

    dt = np.dtype([('time', '<Q'),
        ('spikes', np.dtype('<h'), (num_samps, num_chans))])
    temp = np.fromfile(f, dt)

    f.close()

    return Spikeset(temp['spikes'] * np.reshape(uvolt_conversion, [1, 1, 4]),
            temp['time'], peak_align, fs)

def load(filename):
    if filename.endswith('.ntt'):
        return loadNtt(filename)
    if filename.endswith('spike'):
        return loadDotSpike(filename)
    return None

# rough breakdown as follows
# spikeset contains spikes, timestamps for the whole ntt file
# feature contains a name, channel count and data
# boundary contains a feature name, and polygon boundary
# cluster contains a list of boundaries, along with some summary statistics

# convention: N number of spikes, C number of channels, L length of waveforms


# Spike data is N x L x C
class Spikeset:
    def __init__(self, spikes, timestamps, peak_index, sampling_frequency):
        self.spikes = spikes
        self.time = timestamps
        self.N = len(timestamps)
        self.peak_index = peak_index
        self.fs = sampling_frequency
        self.dt_ms = 1000.0 / self.fs
        self.features = [features.Feature_Peak(self), features.Feature_Energy(self),
            features.Feature_Time(self), features.Feature_Valley(self),
            features.Feature_Barycenter(self), features.Feature_FallArea(self),
            features.Feature_PCA(self)]
        self.T = (max(self.time) - min(self.time)) / 1e6

    def __del__(self):
        print "Spikeset object being destroyed"

    def featureNames(self):
        return [feature.name for feature in self.features]

    def featureByName(self, name):
        for feature in self.features:
            if feature.name == name:
                return feature
        return None


# Clusters have a color, a set of boundaries and some calculation functions
class Cluster:
    def __init__(self, spikeset):
        self.color = (random.randrange(100, 200), random.randrange(100, 200),
                random.randrange(100, 200))
        self.member = np.array([False] * spikeset.N)
        self.bounds = []
        self.wave_bounds = []
        self.add_bounds = []
        self.del_bounds = []
        self.isi = []
        self.refractory = np.array([False] * spikeset.N)
        self.stats = {}

    def __del__(self):
        print "Cluster object being destroyed"

    def addBound(self, bound):
        # bounds should be a 5 or 6 length tuple,
        # (name_x, chan_x, name_Y, chan_y, polygon_points, extra_data)
        self.removeBound(bound[0], bound[1], bound[2], bound[3])
        self.bounds.append(bound)

    def getBoundPolygon(self, feature_name_x, feature_chan_x, feature_name_y,
            feature_chan_y, boundtype='limits'):
        f = lambda bound: ((bound[0] == feature_name_x) and
            (bound[1] == feature_chan_x) and (bound[2] == feature_name_y) and
            (bound[3] == feature_chan_y))
        if boundtype == 'add':
            temp = [bound[4] for bound in self.add_bounds if f(bound)]
        elif boundtype == 'del':
            temp = [bound[4] for bound in self.del_bounds if f(bound)]
        elif boundtype == 'limits':
            temp = [bound[4] for bound in self.bounds if f(bound)]

        if temp:
            return temp
        return None

    def removeBound(self, feature_name_x, feature_chan_x, feature_name_y,
            feature_chan_y, boundtype='limits'):
        f = lambda bound: not ((bound[0] == feature_name_x) and
            (bound[1] == feature_chan_x) and (bound[2] == feature_name_y)
            and (bound[3] == feature_chan_y))

        if boundtype == 'add':
            self.add_bounds = [bound for bound in self.add_bounds if f(bound)]
        elif boundtype == 'del':
            self.del_bounds = [bound for bound in self.del_bounds if f(bound)]
        elif boundtype == 'limits':
            self.bounds = [bound for bound in self.bounds if f(bound)]

    def calculateMembership(self, spikeset):
        self.mahal_valid = False
        if (not self.bounds) and (not self.add_bounds):
            self.member = np.array([False] * spikeset.N)
            self.isi = []
            self.refractory = np.array([False] * spikeset.N)
            return

        if self.add_bounds:
            # if we have additive boundaries, use those
            self.member = np.array([False] * spikeset.N)
            for bound in self.add_bounds:
                px = [feature.data[:, bound[1]] for feature in
                    spikeset.features if feature.name == bound[0]][0]
                py = [feature.data[:, bound[3]] for feature in
                    spikeset.features if feature.name == bound[2]][0]

                data = np.column_stack((px, py))

                self.member = np.logical_or(self.member,
                    nx.points_inside_poly(data, bound[4]))
        else:
            self.member = np.array([True] * spikeset.N)

        # start with everything and cut it down
        for bound in self.bounds:
            px = [feature.data[:, bound[1]] for feature in
                spikeset.features if feature.name == bound[0]][0]
            py = [feature.data[:, bound[3]] for feature in
                spikeset.features if feature.name == bound[2]][0]

            data = np.column_stack((px, py))

            self.member = np.logical_and(self.member,
                nx.points_inside_poly(data, bound[4]))

        for (chan, sample, lower_bound, upper_bound) in self.wave_bounds:
            w = np.logical_and(
                spikeset.spikes[:, sample, chan] >= lower_bound,
                spikeset.spikes[:, sample, chan] <= upper_bound)
            self.member = np.logical_and(self.member, w)


        t = spikeset.time[self.member]
        self.refr_period = 2
        self.burst_period = 20
        self.isi = (t[1:] - t[0:-1]) / 1e3

        self.refractory = np.array([False] * spikeset.N)
        #ref = np.logical_and(self.isi < self.refr_period, self.isi> 0.8)
        ref = self.isi < self.refr_period
        self.refractory[self.member] = np.logical_or(np.append(ref, False),
                np.append(False, ref))

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
            self.stats['burst'] = (np.sum(self.isi < self.burst_period) /
                (self.stats['num_spikes'] - 1) * 100)

            self.stats['refr_count'] = np.sum(self.isi < self.refr_period)

            self.stats['refr_frac'] = (float(self.stats['refr_count']) /
                (self.stats['num_spikes'] - 1))

            alpha = (self.stats['mean_rate'] * self.stats['refr_count'] /
                (self.stats['num_spikes'] * 2.0 * self.refr_period * 1e-3))
            if alpha > 0.25:
                self.stats['refr_fp'] = 100.0
            else:
                self.stats['refr_fp'] = 100 * 0.5 * (1 - np.sqrt(1 - 4 * alpha))

            # csi and isolation needs peaks
            p = features.Feature_Peak(spikeset).data
            if self.stats['num_spikes'] * 2 > spikeset.N:
                self.stats['isolation'] = np.NAN
            else:
                try:
                    cvi= np.linalg.inv(np.cov(np.transpose(p[self.member, :])))
                    u = p-np.mean(p[self.member, :], axis=0)
                    m = np.sum(np.dot(u, cvi) * u, axis=1)
                    m = np.sort(m)
                    self.stats['isolation'] = m[self.stats['num_spikes']* 2 - 1]
                except Exception:
                    self.stats['isolation'] = np.NAN

            chan = np.round(np.mean(np.argmax(p[self.member,:], axis=1)))
            delta = np.diff(p[self.member, chan])

            delta= delta[np.logical_and(
                self.isi < self.burst_period,
                self.isi > self.refr_period)]
            if np.size(delta,0):
                self.stats['csi'] = ((100.0 / np.size(delta, 0)) *
                    (np.sum(delta <= 0) - np.sum(delta > 0)))
            else:
                self.stats['csi'] = np.NAN

    def calculateMahal(self, spikeset):
        if np.all(np.logical_not(self.member)):
            return
        # Work on this as a separate tool, too slow every click
        # Compute mahal distance for all waveforms in cluster
        try:
            temp = np.concatenate([spikeset.spikes[self.member, :, i] for i in
                range(np.size(spikeset.spikes, 2))], axis=1)

            cvi= np.linalg.inv(np.cov(np.transpose(temp)))
            u = temp-np.mean(temp, axis=0)
            m = np.sum(np.dot(u, cvi) * u, axis=1)
            self.mahal = m
            self.mahal_valid = True
        except Exception:
            self.mahal = np.NAN


if __name__ == "__main__":
    print "Loading dotspike"
#    ss = loadDotSpike('TT22.spike')
    ss = loadNtt('TT6_neo.ntt')

    d = ss.featureByName('PCA').data

    import pylab

    pylab.ion()

#    pylab.plot(ss.spikes[1,:,])
#    pylab.show()

    from sklearn import mixture
    g = mixture.GMM(n_components=8);
    g.fit(d[:,1:2])
    label = g.predict(d[:, 1:2])


    pylab.scatter(d[0:500,1], d[0:500,2], s=5, c=label[0:500])

    import pdb; pdb.set_trace()
    g = mixture.GMM(n_components=10, cvtype='full')

    g.fit(data)
    Out[171]: GMM(cvtype='full', n_components=10)

    label = g.predict(data)

    scatter(data[:,1], data[:,2], s=5, edgecolor='none', c=label)
    Out[173]: <matplotlib.collections.PathCollection at 0x19e88128>

    scatter(data[:,0], data[:,1], s=5, edgecolor='none', c=label)
    Out[174]: <matplotlib.collections.PathCollection at 0x1a52d3c8>

    scatter(data[:,0], data[:,2], s=5, edgecolor='none', c=label)
    Out[175]: <matplotlib.collections.PathCollection at 0x13fa79b0>

    scatter(data[:,0], data[:,3], s=5, edgecolor='none', c=label)
    Out[176]: <matplotlib.collections.PathCollection at 0x1a313160>

    scatter(data[:,1], data[:,3], s=5, edgecolor='none', c=label)
    Out[177]: <matplotlib.collections.PathCollection at 0x19da7b70>

    #loadDotSpike('TT22.spike')
    #spikeset = loadNtt('Sample2.ntt')
    #clust = Cluster(spikeset)
    # (name_x, chan_x, name_Y, chan_y, polygon_points, extra_data)
    #clust.addBound(('Peak', 0, 'Peak', 1, [[150, 100], [150, 150],
    #    [100, 150], [100, 100]]))

    #clust.calculateMembership(spikeset)

    #pyplot.hist(np.sqrt(clust.mahal))

    #x = np.array(range(1000))

    #plot(chi2f(x,127))

    #m = clust.mahal
    #hold(False)
    #hist(m, normed=True)
    #hold(True)
    #plot(x, chi2f(x,127), 'g', linewidth=3)

    #test = mvtf(x, 1, 128, clust.cvd)
    #plot(x, test, 'r', linewidth=3)

    #k = 20;
    #p = 128;
    #cvd = clust.cvd

    #del spikeset
    pass
