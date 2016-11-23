#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:27:15 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

def plotKeogram(t, y, kg, title=None, legend=None, ylim=None, pcolorbar=None):
    """
    """
    t_l = len(t)
    
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.pcolormesh(t, y, np.nan_to_num(kg.T))
    plt.xlabel('UT')
    plt.ylabel('Elevation [deg]')
    plt.title(title)
    if ylim is not None:
        ax.set_ylim(ylim)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        plt.colorbar()
        
    ax.xaxis.set(major_formatter=formatter)
    plt.show()