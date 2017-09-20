#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 14:56:37 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import time, datetime
import pandas
from pandas import read_hdf
import h5py
import yaml
import pyGps


############################# Initials #########################################

receiver = {2:'mah22800.h5', 3:'mah32800.h5', 4:'mah42800.h5', 5:'mah52800.h5',
            6:'mah62800.h5', 7:'mah72800.h5', 8:'mah82800.h5', 9:'mah92800.h5',
            13:'ma132800.h5'}
            
yml = {2:'mah22800.yaml', 3:'mah32800.yaml', 4:'mah42800.yaml', 5:'mah52800.yaml',
            6:'mah62800.yaml', 7:'mah72800.yaml', 8:'mah82800.yaml', 9:'mah92800.yaml',
            13:'ma132800.yaml'}

fn = '/home/smrak/Documents/TheMahali/data/rinex/'
navfname = '/home/smrak/Documents/TheMahali/gnss/gps/brdc2800.15n'
skip = 500
#sv = 9
green_alt = 120
red_alt = 250
blue_alt = 90
#rx = 8
dt = 17
def writePlottingHDF(sv=9):
    h5file = h5py.File('/home/smrak/Documents/TheMahali/plotting1/sv9.h5', 'w')
    data = read_hdf(fn+ str(receiver.get(8)))
    obstimes = np.array((data.major_axis))
    obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=dt)
    obstimes = obstimes[skip:]
    posixtimes = pyGps.datetime2posix(obstimes)
    h5file.create_dataset('time', data=posixtimes)
    posixtimes = pyGps.datetime2posix(obstimes)
    for rx in receiver:
        fname = fn + str(receiver.get(rx))
        stream = yaml.load(open(fn + str(yml.get(rx)), 'r'))
        rx_xyz = stream.get('APPROX POSITION XYZ')
        data = read_hdf(fname)
        lli2 = np.array(data['L2', sv, skip:, 'lli'])
        idx = np.where(lli2%2)[0]
        lli = np.zeros(lli2.shape[0])
        lli[idx] = 1
        
        obstimes = np.array((data.major_axis))
        obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=dt)
        obstimes = obstimes[skip:]
        
    #    posixtimes = pyGps.datetime2posix(obstimes)
        ipp_green = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, green_alt, navfname, cs='wsg84')
        
        ipp_red = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, red_alt, navfname, cs='wsg84')
        ipp_blue = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, blue_alt, navfname, cs='wsg84')
        
        mah = h5file.create_group('mah'+str(rx))
        mah.create_dataset('lli', data=lli)
        green = mah.create_group('120km')
        red = mah.create_group('250km')
        blue = mah.create_group('90km')
        blue.create_dataset('lat', data=ipp_blue[0])
        blue.create_dataset('lon', data=ipp_blue[1])
        red.create_dataset('lat', data=ipp_red[0])
        red.create_dataset('lon', data=ipp_red[1])
        green.create_dataset('lat', data=ipp_green[0])
        green.create_dataset('lon', data=ipp_green[1])
    h5file.close()
    
def writePlottingHDF2():
    h5file = h5py.File('/home/smrak/Documents/TheMahali/plotting1/dub.h5', 'w')
    data = read_hdf(fn+ str(receiver.get(8)))
    obstimes = np.array((data.major_axis))
    obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=dt)
    obstimes = obstimes[skip:]
    posixtimes = pyGps.datetime2posix(obstimes)
    h5file.create_dataset('time', data=posixtimes)
    posixtimes = pyGps.datetime2posix(obstimes)
    gps = [9, 23]
    for sv in gps:
        for rx in receiver:
            fname = fn + str(receiver.get(rx))
            stream = yaml.load(open(fn + str(yml.get(rx)), 'r'))
            rx_xyz = stream.get('APPROX POSITION XYZ')
            data = read_hdf(fname)
            lli2 = np.array(data['L2', sv, skip:, 'lli'])
            idx = np.where(lli2%2)[0]
            lli = np.zeros(lli2.shape[0])
            lli[idx] = 1
            
            obstimes = np.array((data.major_axis))
            obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=dt)
            obstimes = obstimes[skip:]
            
        #    posixtimes = pyGps.datetime2posix(obstimes)
            ipp_green = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, green_alt, navfname, cs='wsg84')
            ipp_red = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, red_alt, navfname, cs='wsg84')
            ipp_blue = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, blue_alt, navfname, cs='wsg84')
            
            mah = h5file.create_group('mah'+str(rx))
            mah.create_dataset('lli', data=lli)
            green = mah.create_group('120km')
            red = mah.create_group('250km')
            blue = mah.create_group('90km')
            blue.create_dataset('lat', data=ipp_blue[0])
            blue.create_dataset('lon', data=ipp_blue[1])
            red.create_dataset('lat', data=ipp_red[0])
            red.create_dataset('lon', data=ipp_red[1])
            green.create_dataset('lat', data=ipp_green[0])
            green.create_dataset('lon', data=ipp_green[1])
    h5file.close()
    
#data = read_hdf(fn+ str(receiver.get(8)))

def wrtiteCollocatedSV2HDF(rx=8):
    stream = yaml.load(open(fn + str(yml.get(rx)), 'r'))
    rx_xyz = stream.get('APPROX POSITION XYZ')
    data = read_hdf('/home/smrak/Documents/TheMahali/data/rinex/mah82800.h5')
    obstimes = np.array((data.major_axis))
    obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=dt)
    posixtimes = pyGps.datetime2posix(obstimes)
    dumb = data['L1', :,1, 'data']
    svlist = dumb.axes
    h5file = h5py.File('/home/smrak/Documents/TheMahali/plotting1/rx8.h5', 'w')
    h5file.create_dataset('time', data=posixtimes)
    for sv in svlist[0]:
        lli2 = data['L2', sv,:, 'lli']
        idx = np.where(lli2%2)[0]
        lli = np.zeros(lli2.shape[0])
        lli[idx] = 1
        aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, obstimes, 0, navfname, cs='aer')
        gr = h5file.create_group(str(sv))
        gr.create_dataset('az', data=aer[0])
        gr.create_dataset('el', data=aer[1])
        gr.create_dataset('lli', data=lli)
    h5file.close()
wrtiteCollocatedSV2HDF()
#writePlottingHDF2()
#writePlottingHDF()
#f = h5py.File('/home/smrak/Documents/TheMahali/plotting1/single.h5', 'r')
#for name in f:
#    print (name)
#f.close()