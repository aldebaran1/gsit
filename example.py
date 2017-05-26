#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:17:17 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
#import plotoptics
import datetime
#import asi
#import pyGps
import pandas
from pandas import read_hdf
from gsit import asi, plotGps, pyGps, plotting, plotoptics, plotSatellite, magnetometer
import yaml
#import plotGps
#import plotting
#import plotSatellite
#import magnetometer

#Type paths to raw data

receiver = {2:'mah22800.h5', 3:'mah32800.h5', 4:'mah42800.h5', 5:'mah52800.h5',
            6:'mah62800.h5', 7:'mah72800.h5', 8:'mah82800.h5', 9:'mah92800.h5',
            13:'ma132800.h5'}
            
yml = {2:'mah22800.yaml', 3:'mah32800.yaml', 4:'mah42800.yaml', 5:'mah52800.yaml',
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
    plotoptics.plot2Keogram(t5, t6, el5, el6, i5/1000, i6/i4, ylim = [45, 85], 
                        pcolorbar=True, ytick2=[45, 55, 65, 75],ytick1=[45,55,65,75,85],
                        cmap='viridis', cbartick1=[0,2,4,6,8],
                        cbartick2=[0,0.4,0.8,1.2,1.6,2], cbartitle1='558nm [kR]',cbartitle2='red/blue',
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
    plotoptics.plot3Keogram(t5, t6, t4, el5, el6, el4, i5/1E3, i6/i4, i4/1E3, 
                        ylim=[45,85], title1 = 'Green line, 588nm',
                        title2 = 'Red line, 588nm', title3 = 'Blue line, 588nm',
                        pcolorbar=True, ytick=[45, 55, 65, 75, 85], cmap='viridis',
                        cbartick1=[0,2,4,6,8], cbartick2=[0,0.5,1,1.5,2,3],
                        cbartick3=[0,0.2,0.6,1, 1.4], 
                        cbartitle1='kR', cbartitle2='', cbartitle3='kR',
                        xtick=[datetime.datetime(2015, 10, 7, 6, 10, 0), 
                               datetime.datetime(2015, 10, 7, 6, 15, 0),
                               datetime.datetime(2015, 10, 7, 6, 20, 0),
                               datetime.datetime(2015, 10, 7, 6, 25, 0),
                               datetime.datetime(2015, 10, 7, 6, 30, 0),
                               datetime.datetime(2015, 10, 7, 6, 35, 0),
                               datetime.datetime(2015, 10, 7, 6, 40, 0)],
                        xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                              datetime.datetime(2015, 10, 7, 6, 40, 0)])
    
#Keogram with 4 subplots

def p4kg():
    lli = np.array(data['L2', sv, :, 'lli'])
    #Get keograms
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t5, el5, i5, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '558', obstimes=obstimes, elevation=ipp[1])
    t4, el4, i4, el_out5 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '428', obstimes=obstimes, elevation=ipp[1])
    t6, el6, i6, el_out6 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                       '630', obstimes=obstimes, elevation=ipp[1])
    plotoptics.plot4Keogram(t5, t6, t4, t4, el5, el6, el4, el4, 
                            i5/1E3, i6/1E3, i4/1E3, i6/i4, 
                        ylim=[45,85], pcolorbar=True, ytick1=[55, 65, 75, 85],
                        ytick2=[55, 65, 75, 85],  ytick3=[45, 55, 65, 75, 85],
                        ytick4=[45, 55, 65, 75], cmap='viridis',
                        cbartick1=[0,2,4,6,8], cbartick2=[0,0.4,0.8,1.2,1.6],
                        cbartick3=[0,0.4,0.8,1.2, 1.6], cbartick4=[0, 0.4,0.8,1.2,1.6],
                        cbartitle1='558nm [kR]', cbartitle2='630nm [kR]', 
                        cbartitle3='428nm [kR]', cbartitle4='red/blue',
                        elevation=ipp[1], obstimes=obstimes, lli=lli,
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
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname)
    t, d = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 120,
                                          timelim, '558', obstimes=obstimes)
    t4, d4 = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 90,
                                          timelim, '428', obstimes=obstimes)
    t6, d6 = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 250,
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
    tec_norm = vTECp #- vTECp[np.isfinite(vTECp[0:1000])].min()
    
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
#    aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, 
#                                             navfname, cs='aer')
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, 130E3, 
                                             navfname, cs='wsg84')
    #Asi intensity:
    t5, intensity5 = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 130,
                                              timelim, '558', obstimes=obstimes)
    t4, intensity4 = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 130,
                                          timelim, '428', obstimes=obstimes)
    t6, intensity6 = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], 130,
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
    t5, intensity5 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                              timelim, '558', obstimes=obstimes)
    t4, intensity4 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '428', obstimes=obstimes)
    t6, intensity6 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                          timelim, '630', obstimes=obstimes)
    # Asi Keogram
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=ipp[1])
    plotting.plot4subplot(t, obstimes, t5, roti_t, i/1E3, tec_norm, np.array(intensity5)/1E3, 
                      roti, el1=el, ylabel1='Elevation [deg]', ylabel2='vTEC[TECU]', 
                      ylabel3='Intensity [kR]', ylabel4='ROTI', 
                      xlabel='UT', ylim1=[45, 85], color5='b', title1='Rx8 - PRN23, Poker Flat 10/07/2015',
                      ytick1=[45,55,65,75,85], pcolorbar=True, cbartick=[0,2,4,6,8],
                      cbartitle='558nm [kR]', t31=t4, y31=np.array(intensity4)/1E3,
                      t32=t6, y32=np.array(intensity6)/1E3, color2='g', color3='b', color4='r',
                      legend2=True, label2='558nm', label3='428nm',label4='630nm', lli1=l2_lli,
                      ylim3=[0,3], ylim4=[0,3.5], ylim2=[0,18], ipp_elevation=aer[1],
                      lli_keo=l2_lli, lw=1,
                      ytick3=[0,0.5,1,1.5,2,2.5], ytick2=[0,4,8,12,16], 
                      ytick4=[0,0.5,1,1.5,2,2.5,3],
                      xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                      datetime.datetime(2015, 10, 7, 6, 40, 0)])
    
def tec_keogram_intensity2():
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
    t5, intensity5 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                            timelim, '558', obstimes=obstimes)
    t4, intensity4 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                            timelim, '428', obstimes=obstimes)
    t6, intensity6 = asi.getAllSkyIntensityAER(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                            timelim, '630', obstimes=obstimes)
    # Asi Keogram
    ipp_aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp_aer[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=ipp[1])
    plotting.plot4subplot(t, obstimes, t5, roti_t, i/1E3, tec_norm, np.array(intensity5)/1E3, 
                      roti, el1=el, ylabel1='Elevation [deg]', ylabel2='vTEC [TECU]', 
                      ylabel3='Intensity [kR]', ylabel4='ROTI[TECU/10s]', 
                      xlabel='UT', ylim1=[45, 85], color5='b', title1='Rx8 - PRN23, Poker Flat 10/07/2015',
                      ytick1=[45,55,65,75,85], pcolorbar=True, cbartick=[0,2,4,6,8],
                      cbartitle='Intensity [kR]', t31=t4, y31=np.array(intensity4)/1E3,
                      t32=t6, y32=np.array(intensity6)/1E3, color2='g', color3='b', color4='r',
                      legend2=True, label2='558nm', label3='428nm',label4='630nm', lli1=l2_lli,
                      ylim3=[0,3.5], ylim4=[0,3.5], ylim2=[0,18], ipp_elevation=aer[1],
                      ytick3=[0,1,2,3], ytick2=[0,4,8,12,16], 
                      ytick4=[0,0.5,1,1.5,2,2.5,3],
                      xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                      datetime.datetime(2015, 10, 7, 6, 40, 0)])

# Test lot-lan and aer intensity
def tec_keogram_intensity3():
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
    aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='wsg84')
    #Asi intensity:
    t5, intensity5 = asi.getAllSkyIntensityAER(asi_folder3, aer[0], aer[1], ipp_alt,
                                            timelim, '558', obstimes=obstimes)
    t5a, intensity5a = asi.getAllSkyIntensity(asi_folder3, ipp[0], ipp[1], ipp_alt,
                                            timelim, '558', obstimes=obstimes)
    # Asi Keogram
    ipp_aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp_aer[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=ipp[1])
    plotting.plot4subplot(t, obstimes, t5, roti_t, i/1E3, tec_norm, np.array(intensity5)/1E3, 
                      roti, el1=el, ylabel1='Elevation [deg]', ylabel2='vTEC [TECU]', 
                      ylabel3='Intensity [kR]', ylabel4='ROTI[TECU/10s]', 
                      xlabel='UT', ylim1=[45, 85], color5='b', title1='Rx8 - PRN23, Poker Flat 10/07/2015',
                      ytick1=[45,55,65,75,85], pcolorbar=True, cbartick=[0,2,4,6,8],
                      cbartitle='Intensity [kR]', t31=t5a, y31=np.array(intensity5a)/1E3,
                      color2='g', color3='r',
                      legend2=True, label2='org', label3='int',label4='630nm', lli1=l2_lli,
                      ylim3=[0,3.5], ylim4=[0,3.5], ylim2=[0,18], ipp_elevation=aer[1],
                      ytick3=[0,1,2,3], ytick2=[0,4,8,12,16], 
                      ytick4=[0,0.5,1,1.5,2,2.5,3],
                      xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                      datetime.datetime(2015, 10, 7, 6, 40, 0)])
    
def tec_keogram_intensity5():
    N_roti = 100
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
    roti = pyGps.getROTI(vTECp, N_roti)
    roti_t = obstimes[0:-N_roti]
    #ASI Intensity data
    # Piercing point for all-sky data 
    aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    #ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='wsg84')
    #Asi intensity:
    t5, intensity5 = asi.getAllSkyIntensityAER(asi_folder3, aer[0], aer[1], ipp_alt,
                                            timelim, '558', obstimes=obstimes)
    t6, intensity6 = asi.getAllSkyIntensityAER(asi_folder3, aer[0], aer[1], ipp_alt,
                                            timelim, '630', obstimes=obstimes)
    t4, intensity4 = asi.getAllSkyIntensityAER(asi_folder3, aer[0], aer[1], ipp_alt,
                                            timelim, '428', obstimes=obstimes)
    # Asi Keogram
    ipp_aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp_aer[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=aer[1])
    # Magnetometer
    mag_folder = '/home/smrak/Documents/TheMahali/magnetometer/poker_2015_10_07.csv'
    mag = magnetometer.readMag(mag_folder, '2015-10-07')
    b_t=mag[0][5000:]
    mag_xyz = magnetometer.magHDZ2XYZ(mag, angle=True)
    plotting.plot5subplot(t, t5, obstimes, roti_t, i/1E3, np.array(intensity5)/1E3, tec_norm,
                      roti, el1=el, ylabel1='Elevation [deg]', ylabel3='vTEC [TECU]', 
                      ylabel2='Intensity [kR]', ylabel4='ROTI[TECU/10s]', 
                      xlabel='UT', ylim1=[45, 85], color5='b',
                      ytick1=[45,55,65,75,85], pcolorbar=True, cbartick=[0,2,4,6,8],
                      cbartitle='Intensity [kR]', t31=t6, t32=t4, y31=np.array(intensity6)/1E3,
                      y32=np.array(intensity4)/1E3, color2='g', color3='r', color4='b',
                      legend1=True, label2='558nm', label3='630nm',label4='428nm', lli1=l2_lli,
                      ylim2=[0,3.5], ylim4=[0,3.5], ylim3=[0,18], ipp_elevation=aer[1],
                      ytick2=[0,1,2,3], ytick3=[0,4,8,12,16], ytick4=[0,0.5,1,1.5,2,2.5,3],
                      Bt=b_t, Bx=mag_xyz[0][5000:]/1000, ylabel5='Bx [uT]',cbx='b',
                      ylim5=[12.2,13], ytick5=[12.2,12.4,12.6,12.8],
                      xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                      datetime.datetime(2015, 10, 7, 6, 40, 0)])
    
def phaseScintillation():
    for rr in receiver:
        print (receiver.get(rr))


def plotObservationDay():
    # Green keogram
    ipp_aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='aer')
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='wsg84')
    t, el, i, el_out4 = asi.getASIKeogramIPP(asi_folder3, ipp_aer[0], ipp_alt, timelim, 
                                             '558', obstimes=obstimes, elevation=ipp[1])
    # Magnetic Field
    mag_folder = '/home/smrak/Documents/TheMahali/magnetometer/poker_2015_10_07.csv'
    mag = magnetometer.readMag(mag_folder, '2015-10-07')
    #mag_start = np.datetime64('2015-10-07T00:00:00')
    #mag_stop = np.datetime64('2015-10-07T22:00:00')
    #idx = np.where( (mag[0] > mag_start) & (mag[0] < mag_stop) )
    mag_xyz = magnetometer.magHDZ2XYZ(mag, angle=True)
    
    plotting.plot2subplotK_xy(t, el, i/1000, mag[0][8000:], mag_xyz[0][8000:]/1000, obstimes=obstimes,
                              ipp=ipp_aer[1], xlabel='UT', ylim1=[45,85],
                              ylabel1='Elevation[deg]', ylabel2='Bx [uT]', ytick1 = [45,55,65,75,85],
                              ylim2=[12.2,13.1],ytick2 = [12,12.2,12.4,12.6,12.8,13],
                              pcolorbar=True, cbartick=[0,2,4,6,8], cbartitle='558nm [kR]', 
                              xlim=[datetime.datetime(2015, 10, 7, 6, 10, 0),
                              datetime.datetime(2015, 10, 7, 6, 40, 0)])

def plotGlobalMap():
    start = datetime.datetime(2015, 10, 7, 6, 10, 0)
    stop = datetime.datetime(2015, 10, 7, 6, 40, 0)
    idx = np.where( (obstimes >= start) & (obstimes <= stop) )[0]
    ipp = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfname, cs='wsg84')
    #print (obstimes[idx])
    plotSatellite.plotSatelliteGlobalMap(rx_xyz, sv, obstimes[idx], 
                                         width = 5000000, height = 3000000, center=[-130, 60])
    plotSatellite.plotIPPGlobalMap(ipp[0][idx], ipp[1][idx], 
                                   width = 5000000, height = 3000000, center=[-130, 60])
    plotSatellite.plotSatelliteGlobalMap(rx_xyz, 23, obstimes[0:6000], 
                                         color='green', lw=2, width = 5000000, height = 3000000, 
                                         center=[-130, 60])
    plotSatellite.plotSatelliteGlobalMap(rx_xyz, 23, obstimes[idx], color='red', 
                                         lw=3, width = 5000000, height = 3000000,center=[-130, 60])
if __name__ == '__main__':
    #pkg()
    #p2kg()
    p3kg()
    #p4kg()
    #pInt()
    #tec()
    #tec_vs_asi()
    #tec_keogram_intensity5()
    #plotGlobalMap()
    #plotObservationDay()
    #phaseScintillation()