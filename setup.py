# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:06:27 2012

@author: Bernard
"""

from distutils.core import setup
import matplotlib
import py2exe


excludes = []
includes = ["scipy.io.matlab.streams"]

opts = {
    "py2exe": {
        "includes":includes,
        "excludes":excludes
    }
}

setup(console=['main.py'], options=opts, data_files=matplotlib.get_py2exe_datafiles())