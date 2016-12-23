#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 17:47:19 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import plotoptics
import datetime
import asi
import pyGps
import pandas
from pandas import read_hdf
import yaml
import plotGps
import plotting
import plotSatellite
import magnetometer

#Type paths to raw data

receiver = {2:'mah22800.h5', 3:'mah32800.h5', 4:'mah42800.h5', 5:'mah52800.h5',
            6:'mah62800.h5', 7:'mah72800.h5', 8:'mah82800.h5', 9:'mah92800.h5',
            13:'ma132800.h5'}
            
yaml = {2:'mah22800.yaml', 3:'mah32800.yaml', 4:'mah42800.yaml', 5:'mah52800.yaml',
            6:'mah62800.yaml', 7:'mah72800.yaml', 8:'mah82800.yaml', 9:'mah92800.yaml',
            13:'ma132800.yaml'}

fn = '/home/smrak/Documents/TheMahali/rinex/'
asi_folder3 = '/home/smrak/Documents/TheMahali/Allsky_multi/'
fname = '/home/smrak/Documents/TheMahali/rinex/mah82800.h5'
navfname = '/home/smrak/Documents/TheMahali/gnss/gps/brdc2800.15n'
timelim = ['10/07/2015', '06:08:00', '10/07/2015', '06:50:00']
sv = 23
rx = 8
ipp_alt = 130
stream = yaml.load(open(fn+'mah82800.yaml', 'r'))
rx_xyz = stream.get('APPROX POSITION XYZ')
data = read_hdf(fname)
obstimes = np.array((data.major_axis))
obstimes = pandas.to_datetime(obstimes) 