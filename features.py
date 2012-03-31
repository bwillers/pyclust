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
        self.channels = np.size(self.data,1)
        
    def calculate(self, spikeset):
        pass
    
class Feature_Peak(Feature):
    def __init__(self, spikeset):
        Feature.__init__(self, 'Peak', spikeset)
        
    def calculate(self, spikeset):
        return np.max(spikeset.spikes, axis=2)
