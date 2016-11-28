#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:17:17 2016

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

#Type paths to raw data
fn = '/home/smrak/Documents/TheMahali/rinex/'
asi_folder3 = '/home/smrak/Documents/TheMahali/Allsky_multi/'
fname = '/home/smrak/Documents/TheMahali/rinex/mah82800.h5'
navfname = '/home/smrak/Documents/TheMahali/gnss/gps/brdc2800.15n'
timelim = ['10/07/2015', '06:10:00', '10/07/2015', '06:50:00']
sv = 23
rx = 8
ipp_alt = 130E3
stream = yaml.load(open(fn+'mah82800.yaml', 'r'))
rx_xyz = stream.get('APPROX POSITION XYZ')
data = read_hdf(fname)
obstimes = np.array((data.major_axis))
obstimes = pandas.to_datetime(obstimes) 



# Plot data
# single keogram
def pkg():
    #Get keograms
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t5, el5, i5, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '558', obstimes=obstimes, elevation=ipp[1])
    t4, el4, i4, el_out5 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '428', obstimes=obstimes, elevation=ipp[1])
    t6, el6, i6, el_out6 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                        '630', obstimes=obstimes, elevation=ipp[1])
    plotoptics.plotKeogram(t5, el5, i5/1E3, ylim=[45,85], title='Green line intensity',
                       pcolorbar=True, cmap='viridis', cbartick=[0,2,4,6,8], 
                       cbartitle='kR', ytick=[45,55,65,75,85],
                       xtick=[datetime.datetime(2015, 10, 7, 6, 10, 0), 
                               datetime.datetime(2015, 10, 7, 6, 15, 0),
                               datetime.datetime(2015, 10, 7, 6, 20, 0),
                               datetime.datetime(2015, 10, 7, 6, 25, 0),
                               datetime.datetime(2015, 10, 7, 6, 30, 0),
                               datetime.datetime(2015, 10, 7, 6, 35, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)],
                       xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                              datetime.datetime(2015, 10, 7, 6, 40, 0)])

# Keogram with 2 subplots
def p2kg():
    #Get keograms
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t5, el5, i5, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '558', obstimes=obstimes, elevation=ipp[1])
    t4, el4, i4, el_out5 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '428', obstimes=obstimes, elevation=ipp[1])
    t6, el6, i6, el_out6 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '630', obstimes=obstimes, elevation=ipp[1])
    plotoptics.plot2Keogram(t6, t5, el6, el5, i6/i4, i5/i4, ylim = [45, 85], 
                        title1 = 'Ratio of red and blue line',  
                        title2 = 'Ratio of green and blue line',  
                        pcolorbar=True, ytick=[45, 55, 65, 75, 85],
                        cmap='viridis', cbartick1=[0,0.4,0.8,1.2,1.6],
                        cbartick2=[0,4,8,12,15], cbartitle1='',cbartitle2='',
                        xtick=[datetime.datetime(2015, 10, 7, 6, 10, 0), 
                               datetime.datetime(2015, 10, 7, 6, 15, 0),
                               datetime.datetime(2015, 10, 7, 6, 20, 0),
                               datetime.datetime(2015, 10, 7, 6, 25, 0),
                               datetime.datetime(2015, 10, 7, 6, 30, 0),
                               datetime.datetime(2015, 10, 7, 6, 35, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)],
                        xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                              datetime.datetime(2015, 10, 7, 6, 40, 0)])
#Keogram with 3 subplots

def p3kg():
    #Get keograms
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t5, el5, i5, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '558', obstimes=obstimes, elevation=ipp[1])
    t4, el4, i4, el_out5 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '428', obstimes=obstimes, elevation=ipp[1])
    t6, el6, i6, el_out6 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '630', obstimes=obstimes, elevation=ipp[1])
    plotoptics.plot3Keogram(t5, t6, t4, el5, el6, el4, i5/1E3, i6/1E3, i4/1E3, 
                        ylim=[45,85], title1 = 'Green line, 588nm',
                        title2 = 'Red line, 588nm', title3 = 'Blue line, 588nm',
                        pcolorbar=True, ytick=[45, 55, 65, 75, 85], cmap='viridis',
                        cbartick1=[0,2,4,6,8], cbartick2=[0,0.4,0.6,0.8,1],
                        cbartick3=[0,0.2,0.6,1, 1.4], 
                        cbartitle1='kR', cbartitle2='kR', cbartitle3='kR',
                        xtick=[datetime.datetime(2015, 10, 7, 6, 10, 0), 
                               datetime.datetime(2015, 10, 7, 6, 15, 0),
                               datetime.datetime(2015, 10, 7, 6, 20, 0),
                               datetime.datetime(2015, 10, 7, 6, 25, 0),
                               datetime.datetime(2015, 10, 7, 6, 30, 0),
                               datetime.datetime(2015, 10, 7, 6, 35, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)],
                        xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                              datetime.datetime(2015, 10, 7, 6, 40, 0)])
# Intensity graph
def pInt():
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t, d = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '558', obstimes=obstimes)
    t4, d4 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '428', obstimes=obstimes)
    t6, d6 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                   timelim, '630', obstimes=obstimes)
    plotoptics.plotIntensity(t,d, ylabel='Intensity [R]', xlabel='UT', label1='558nm',
                         t2=t4, y2=d4, t3=t6, y3=d6, label2='428nm', label3='630nm',
                         color1='g', color2 ='b', color3='r', 
                         title='All-sky imager intensity', ylim=[300, 2500], legend=True,
                         xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)])
# TEC
def tec():
    aer = pyGps.getSatellitePosition(np.asarray(rx_xyz), sv, obstimes, navfname, cs='aer')
    L1 = np.array(data['L1', sv, :, 'data'])
    L2 = np.array(data['L2', sv, :, 'data'])
    P1 = np.array(data['C1', sv, :, 'data'])
    P2 = np.array(data['P2', sv, :, 'data'])
    #Phase corrected TEC
    sTECp = pyGps.getPhaseCorrTEC(L1, L2, P1, P2)
    vTECp = pyGps.getVerticalTEC(sTECp, aer[1], 130)
    tec_norm = vTECp - vTECp[np.isfinite(vTECp[0:1000])].min()
    
    plotGps.plotTEC(obstimes, tec_norm, xlabel='UT', ylabel='vTEC [TECU]',
                    ylim=[0,17], ytick=[0,5,10,14,15,16], color = '--m',
                    xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)])

# TEC vs ASI Intensety
def tec_vs_asi():
    aer = pyGps.getSatellitePosition(np.asarray(rx_xyz), sv, obstimes, navfname, cs='aer')
    L1 = np.array(data['L1', sv, :, 'data'])
    L2 = np.array(data['L2', sv, :, 'data'])
    P1 = np.array(data['C1', sv, :, 'data'])
    P2 = np.array(data['P2', sv, :, 'data'])
    l2_lli = np.array(data['L2', sv, :, 'lli'])
    #Phase corrected TEC
    sTECp = pyGps.getPhaseCorrTEC(L1, L2, P1, P2)
    vTECp = pyGps.getVerticalTEC(sTECp, aer[1], 130)
    tec_norm = vTECp - vTECp[np.isfinite(vTECp[0:1000])].min()
    #ASI Intensity data
    # Piercing point for all-sky data 
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    #Asi intensity:
    t5, intensity5 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                              timelim, '558', obstimes=obstimes)
    t4, intensity4 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '428', obstimes=obstimes)
    t6, intensity6 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '630', obstimes=obstimes)
    #Plot 2 subplots
    plotting.plot2subplot(obstimes, t5, tec_norm, (np.array(intensity5)/1E3), 
                          ylabel1='vTEC [TECU]', ylabel2='Intensity [kR]', 
                          xlabel='UT', lli1=l2_lli, ylim2=[0.2, 2.6], lm='line', 
                          color2='g', lli2=l2_lli, t21=t4, y21=(np.array(intensity4)/1E3),
                          t22=t6, y22=(np.array(intensity6)/1E3), label4='630nm', color4='r',
                          color3='b', label2='558nm', label3='428nm', legend2=True,
                          xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                          datetime.datetime(2015, 10, 7, 6, 40, 0)])
    
def tec_keogram_intensity():
    aer = pyGps.getSatellitePosition(np.asarray(rx_xyz), sv, obstimes, navfname, cs='aer')
    L1 = np.array(data['L1', sv, :, 'data'])
    L2 = np.array(data['L2', sv, :, 'data'])
    P1 = np.array(data['C1', sv, :, 'data'])
    P2 = np.array(data['P2', sv, :, 'data'])
    l2_lli = np.array(data['L2', sv, :, 'lli'])
    #Phase corrected TEC
    sTECp = pyGps.getPhaseCorrTEC(L1, L2, P1, P2)
    vTECp = pyGps.getVerticalTEC(sTECp, aer[1], 130)
    tec_norm = vTECp - vTECp[np.isfinite(vTECp[0:1000])].min()
    # Rate of TEC Index
    roti = pyGps.getROTI(vTECp, 10)
    roti_t = obstimes[0:-10]
    #ASI Intensity data
    # Piercing point for all-sky data 
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    #Asi intensity:
    t5, intensity5 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                              timelim, '558', obstimes=obstimes)
    t4, intensity4 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '428', obstimes=obstimes)
    t6, intensity6 = asi.getAllskyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '630', obstimes=obstimes)
    # Asi Keogram
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=ipp[1])
    plotting.plot4subplot(t, obstimes, t5, roti_t, i/1E3, tec_norm, np.array(intensity5)/1E3, 
                      roti, el1=el, ylabel1='Elevation [deg]', ylabel2='vTEC [TECU]', 
                      ylabel3='Intensity [kR]', ylabel4='ROTI[TECU/10s]', 
                      xlabel='UT', ylim1=[45, 85], color5='b', title1='Rx8 - PRN23, Poker Flat 10/07/2015',
                      ytick1=[45,55,65,75,85], pcolorbar=True, cbartick=[0,2,4,6,8],
                      cbartitle='Intensity [kR]', t31=t4, y31=np.array(intensity4)/1E3,
                      t32=t6, y32=np.array(intensity6)/1E3, color2='g', color3='b', color4='r',
                      legend2=True, label2='558nm', label3='428nm',label4='630nm', lli1=l2_lli,
                      ylim3=[0,3.5], ylim4=[0,3.5], ylim2=[0,20], ipp_elevation=aer[1],
                      ytick2=[0,4,8,12,16], ytick3=[0,1,2,3],
                      ytick4=[0,0.5,1,1.5,2,2.5,3],
                      xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                      datetime.datetime(2015, 10, 7, 6, 40, 0)])

if __name__ == '__main__':
    #pkg()
    #p2kg()
    #p3kg()
    #pInt()
    #tec()
    #tec_vs_asi()
    tec_keogram_intensity()