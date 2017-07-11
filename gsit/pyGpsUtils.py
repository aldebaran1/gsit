#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:49:04 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plotGpsPolarTrajectory(azimuth=[],elevation=[],labels=None,timelim=None):
    """
    Draw a polar plot with trajectories of all given satellites. Azimuth and elevation
    are obligatory inputs in a form of a list!
    """
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='polar')
    for i in range(len(azimuth)):
        if labels is not None:
            ax.plot(np.radians(azimuth[i]), 90-elevation[i], label='sv'+str(labels[i]))
        else:
            ax.plot(np.radians(azimuth[i]), 90-elevation[i])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rmax(80)
    ax.set_rticks([80, 60, 40, 20])
    ax.set_yticklabels([20, 40, 60, 80])
    plt.legend(bbox_to_anchor=(1.1, 1.1))
    if timelim is not None:
        plt.title('Boston: ' + timelim[0].strftime('%m/%d/%y-%H:%M :: ') + timelim[1].strftime('%m/%d/%y-%H:%M'))
    plt.show()
    
def plotGpsMapTrajectory(lat=[], lon=[], rx=None, labels=None, 
                         timelim=None, totalityc=[], totalityu=[], totalityd=[],
                         latlim=[41, 44], ms=10, color='k',
                         lonlim=[-74, -69], center=[42.36, -71.06],
                         parallels=[42,44], meridians = [-73, -70, -67],
                         ax=None, m=None):
    """
    """
    if ax is None:
        (fig,ax) = plt.subplots(1,1,facecolor='w')
        m = Basemap(llcrnrlat=latlim[0],urcrnrlat=latlim[1],
                    llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],
                    projection='merc', resolution='i', ax=ax)
        
        m.drawparallels(parallels,labels=[False, True, True, True], linewidth=1)
        m.drawmeridians(meridians,labels=[True,True,False,True], linewidth=1)
        
        if len(totalityc) > 0:
            x,y = m(totalityc[1], totalityc[0])
            m.plot(x,y, lw=2, color='r')
        if len(totalityu) > 0:
            x,y = m(totalityu[1], totalityu[0])
            m.plot(x,y, lw=2, color='b')
        if len(totalityd) > 0:
            x,y = m(totalityd[1], totalityd[0])
            m.plot(x,y, lw=2, color='b')
            
            
        if len(lat) > 0 and len(lon) > 0:
            for i in range(len(lat)):
                idx = np.where(np.isfinite(lon[i]))[0]
                x,y = m(lon[i][idx], lat[i][idx])
                if labels is not None:
                    m.plot(x,y, lw=2, label='sv'+str(labels[i]))
                else:
                    m.plot(x,y, lw=2)
        if (rx is not None) and isinstance(rx, np.ndarray):
            if len(rx.shape) > 1:
                for i in range(rx.shape[1]):
                    x,y = m(rx[1][i], rx[0][i])
                    m.scatter(x, y, marker='o', color=color, s=ms)
            else:
                x,y = m(rx[1], rx[0])
                m.scatter(x, y, marker='o', color=color, s=ms)
                
        m.drawcoastlines()
        m.drawstates()
        
    else:
        if len(lat) > 0 and len(lon) > 0:
            for i in range(len(lat)):
                idx = np.where(np.isfinite(lon[i]))[0]
                x,y = m(lon[i][idx], lat[i][idx])
                if labels is not None:
                    m.plot(x,y, lw=2, label='sv'+str(labels[i]))
                else:
                    m.plot(x,y, lw=2)
        if (rx is not None) and isinstance(rx, np.ndarray):
            if len(rx.shape) > 1:
                for i in range(rx.shape[1]):
                    x,y = m(rx[1][i], rx[0][i])
                    m.scatter(x, y, marker='o', color=color, s=ms)
            else:
                x,y = m(rx[1], rx[0])
                m.scatter(x, y, marker='o', color=color, s=ms)
    return ax, m
#    plt.tight_layout()
#    plt.legend(bbox_to_anchor=(1.1, 1.1))
#plotGpsMapTrajectory()