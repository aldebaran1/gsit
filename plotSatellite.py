# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 11:43:32 2016

@author: smrak
"""

from mpl_toolkits.basemap import Basemap
#import numpy as np
import pyGps
import matplotlib.pyplot as plt

def plotSatelliteGlobalMap(rx, sv, obstimes, color = 'r', lw = 3, center=[-130., 50]):
    #Add optional marks on plot
    all_sky_location = [-147.43, 65.12] 
    
    #Create a Basemap object:
    m = Basemap(width=10000000, height=7000000, projection='aeqd',
                    resolution='c', lat_0=center[1], lon_0=center[0])
    #Plot parallels and meridians
    parallels = [20., 40., 60.]
    m.drawparallels(parallels,labels=[False, True, True, True], linewidth=2)
    meridians = [-50, -90, -120, -150., -180., 150.]
    m.drawmeridians(meridians,labels=[True,True,False,True], linewidth=2)
    m.bluemarble()
    
    if isinstance(sv, int):
        sv_llt = pyGps.getSatellitePosition(rx, sv, obstimes, cs='wsg84')
        
        asi_x, asi_y = m(all_sky_location[0], all_sky_location[1])
        sat_x, sat_y = m(sv_llt[1], sv_llt[0])
        m.plot(sat_x, sat_y, color, linewidth=lw)
        m.scatter(asi_x, asi_y, marker='*', color='m', s=25)  # plot a blue dot there
    else:
        asi_x, asi_y = m(all_sky_location[0], all_sky_location[1])
        m.scatter(asi_x, asi_y, marker='*', color='m', s=25)  # plot a blue dot there
        for sat in sv:
            sv_llt = pyGps.getSatellitePosition(rx, sat, obstimes, cs='wsg84')
            sat_x, sat_y = m(sv_llt[1], sv_llt[0])
            m.plot(sat_x, sat_y, color, linewidth=lw)
    plt.show()