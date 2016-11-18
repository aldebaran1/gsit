# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 18:14:50 2016

@author: smrak
"""

import numpy as np
from datetime import date
import mahaliPlot

folder = '/home/smrak/Documents/TheMahali/magnetometer/'
date='2015-10-07'

def readMag(fname):   
    a,b,c,d,e = np.loadtxt(folder+fname,dtype='S',unpack=True,delimiter=',')
    t = np.empty(len(a),dtype='datetime64[s]')
    for i in range(len(a)):
        t[i]=np.datetime64(date+'T'+a[i])
    H = c.astype(float)
    D = d.astype(float)
    Z = e.astype(float)
    #mahaliPlot.plot(H)
    return t, H, D, Z

def magTransform(H, D, Z):
    X = H * np.cos(D)
    Y = H * np.sin(D)
    Z = Z
    return X, Y, Z

def _DtoDeg(d, h):
    D = np.arctan(d/h)
    return D

def getTotalB(x, y, z):
    F = np.sqrt(x**2 + y**2 + z**2)
    return F

def getMagData(fname):
    #start = interval[0]
    #stop = interval[1]
    #fname = 'poker_2015_10_07.csv'
    t,h,d,z = readMag(fname)
    #dt_start = np.datetime64(date+'T'+start)
    #dt_stop = np.datetime64(date+'T'+stop)
    #a = np.where(t == dt_start) [0]
    #b = np.where(t == dt_stop) [0]
    D = _DtoDeg(d, h)    
    
    x, y, z = magTransform(h, D, z)
    F = getTotalB(x, y, z)
    #data = np.vstack((h[a:b], D[a:b], z[a:b]))
    data_xyz = np.vstack((x, y, z, F))

        
    return t, data_xyz
if __name__ == '__main__':
    mag('06:10:00', '06:40:00')    