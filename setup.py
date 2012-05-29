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
data_files = matplotlib.get_py2exe_datafiles()

# something weird with enthought i need to include these manually
data_files.append('mk2_core.dll')
data_files.append('mk2_mc3.dll')
data_files.append('mk2iomp5md.dll')

opts = {
    "py2exe": {
        "includes":includes,
        "excludes":excludes
#        "dll_excludes":["MSVCP90.dll"]
    }
}

setup(console=['main.py'], options=opts, data_files=data_files)
