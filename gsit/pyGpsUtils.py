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
    
def plotGpsMapTrajectory(azimuth=[], elevation=[], rx=[], labels=None, altitude=120, 
                         timelim=None, latlim=[], lonlim=[], center=[]):
    """
    """
    (fig,ax) = plt.subplots(1,1,figsize=(12,12),facecolor='w')
    m = Basemap(llcrnrlat=latlim[0],urcrnrlat=latlim[1],
                llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],
                projection='merc', resolution='i', ax=ax)
    
    parallels = [42,44]
    m.drawparallels(parallels,labels=[False, True, True, True], linewidth=1)
    meridians = [-73, -70, -67]
    m.drawmeridians(meridians,labels=[True,True,False,True], linewidth=1)
    
    x,y = m(center[1], center[0])
    m.scatter(x, y, marker='o', color='m', s=50)
    m.drawcoastlines()
    m.drawstates()
plotGpsMapTrajectory(latlim=[41, 44], lonlim=[-74, -69], center=[42.36, -71.06])