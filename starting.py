# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/bwillers/.spyder2/.temp.py
"""

import matplotlib.pyplot as plt
import numpy as np
import struct

isi_bins = np.logspace(2, 1000, 50)
isi_bin_centers = (isi_bins[0:-1] + isi_bins[1:])
temp = np.histogram([2,5,100,500], isi_bins)
plt.plot(temp[1][0:-1], temp[0])

cf = lambda s,p: [ s[i:i+p] for i in range(0,len(s),p) ]


f = open('Sample3.ntt', 'rb')

contents = f.read()

len(contents)

# A tetrode spike record is as folows:
# uint64 - timestamp                    bytes 0:8
# uint32 - acquisition entity number    bytes 8:12
# uint32 - classified cel number        bytes 12:16
# 8 * uint32- params                    bytes 16:48
# 32 * 4 * int16 - waveform points
# hence total record size is 2432 bits, 304 bytes

# header is supposed to be 16kbyte, i.e. 16 * 2^10 = 2^14

# Read  he headr and findthe conversion factors
header = contents[0:2**14-1]
for line in header.split('\n'):
    if line.strip().startswith('-ADBitVolts'):
        a2d_conversion = 1e6 * np.array(map(float, line.split()[1:5]))
        break

# Read the spikes

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

times = timestamps


isi = np.diff(times, 1, 2) / 1e3
isi = isi / 1e3
isi_bins = np.logspace(np.log10(2), np.log10(1000), 50)
isi_bin_centers = (isi_bins[0:-1] + isi_bins[1:])
isi_bin_count = np.histogram(isi, isi_bins)[0]

# find spike peaks
peaks = spikes[:,:,8]
projections = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]

for i in range(1,7):
    plt.subplot(2,3,i)
    plt.scatter(peaks[:,projections[i-1][0]], peaks[:,projections[i-1][1]], marker='o', s=1)
    plt.xlim([0,300])
    plt.ylim([0,300])
    

plt.show()
print "Done scatter plotting"

np.histogramdd()
test = np.histogramdd(peaks[:,0:2], [np.linspace(0, 250, 100), np.linspace(0, 250, 100)])
plt.figure(2);

np.shape(test[0])
np.shape(test[1])

plt.pcolor(test[1][0], test[1][1], test[0])


bins_x = np.linspace(0, np.max(peaks[:,0]), 100)
bins_y = np.linspace(0, np.max(peaks[:,1]), 100)
count = np.histogram2d(peaks[:,0], peaks[:,1], [bins_x, bins_y])[0]


plot(np.transpose(spikes[:,1,:]), 'b', linewidth=1)


range(1,33)
np.shape(np.std(spikes[:,1,:], axis=0))
np.shape(np.mean(spikes[:,1,:], axis=0))

fig = plt.figure()
ax = fig.add_subplot(2,2,1)
ax.plot(spikes[1,1,:])


print fig.axes
    
    
    
    
    
    
            
            height = self.mp_proj.figure.bbox.height
            width = self.mp_proj.figure.bbox.width
            
            
            width_a = self.mp_proj.axes.bbox.width

            x0 = self.zoom_startpos[0]
            y0 = height_a - self.zoom_startpos[1]
            x1 = event.xdata
            y1 = height_aevent.ydata

            
            x0 = self.zoom_startpos[0]
            y0 = self.zoom_startpos[1]
            x1 = event.xdata
            y1 = event.ydata
            
            rect = [ int(val) for val in min(x0,x1), max(y0,y1), abs(x1-x0), abs(y1-y0)]
            print rect
            
            # now make them relative
            lim_width = float(np.diff(self.mp_proj.axes.get_xlim()));
            lim_height = float(np.diff(self.mp_proj.axes.get_ylim()));
            
            rect = [rect[0]/lim_width, rect[1]/lim_height, rect[2]/lim_width, rect[3]/lim_height]
            print rect            
            
            start_l = self.mp_proj.figure.subplotpars.left * width
            start_b = self.mp_proj.figure.subplotpars.bottom * height
            #width = width - start_l
            #height = height - start_b
            
            rect = [start_l + rect[0]*(width-start_l), height-rect[1] * height, rect[2] * (width-start_l), rect[3] * (height-start_b)]
            self.mp_proj.drawRectangle(rect)
        #print "mouse move", event.xdata, event.ydata