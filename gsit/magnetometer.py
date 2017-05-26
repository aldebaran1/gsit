#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 18:14:50 2016

@author: Sebastijan Mrak smrak@bu.edu
"""

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

def readMag(fname, date):
    """
    Sebastijan Mrak
    Read CSV file with raw magentometer data. 
    Input: fneme: total path/to/file/file.csv as string
    Output: np array of timestamp and three components [H,D,Z]
    """
    a,b,c,d,e = np.loadtxt(fname,dtype='S',unpack=True,delimiter=',')
    t = np.empty(len(a),dtype='datetime64[s]')
    a = a.astype(str)
    b = b.astype(str)

    for i in range(len(a)):
        t[i]=np.datetime64(date+'T'+a[i])
    H = c.astype(float)
    D = d.astype(float)
    Z = e.astype(float)
    
    return np.array([t, H, D, Z])

def magHDZ2XYZ(F, angle=False):
    """
    Sebastijan Mrak
    Magnetometer components transformation from H, D, Z to X, Y, Z format.
    Inputs: H,D,Z components as np.array. If D component is angle instead of abs. value,
            set optional argument 'angle' to True
    Output: np array with X, Y, Z components
    """
    H = F[1]
    D = F[2]
    Z = F[3]
    if angle:
        D = _DtoDeg(D, H) #From angle to abs. value in [nT]

    X = H * np.cos(D)
    Y = H * np.sin(D)
    Z = Z
    
    return np.array([X, Y, Z])

def _DtoDeg(d, h):
    """
    Sebastijan Mrak
    Transform angle value of D component of geomagnetic field to absolute value.
    Transform follows derivation from 'Introduction to geomagnetic fileds' by 
    Wallace Campbell, Cambridge, June 13, 1997
    """
    D = np.zeros(len(d))
    for i in range(len(d)):
        D[i] = np.arctan(d[1]/h[1])
    return D

def getTotalB(x, y, z):
    """
    Sebastijan Mrak
    Function returns total magnitude of the geomagnetic field. Input parameters 
    can be in either frame of reference, but they must be component's absolute 
    values. 
    """
    F = np.sqrt(x**2 + y**2 + z**2)
    return F

def plotMagnetometer1(t, y1, y2=None, y3=None, y4=None, xlim=None, xlabel=None,
                      ylabel=None, title=None, ylim=None, legend=None, 
                      color1='b', color2='k', color3='g', color4='r'):
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, y1, color1)
    if y2 is not None:
        ax2 = ax1.twinx()
        ax2.plot(t, y2, color2)
        ax2.set_ylabel('By [nT]')
    if y3 is not None:
        ax1.plot(t, y3, color3)
    
    ax1.set_ylabel('Bx [uT]')
    
    plt.xlabel('UT')
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()