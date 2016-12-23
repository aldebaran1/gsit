# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 20:21:46 2016

@author: smrak
"""

import numpy as np
import datetime 
import plotting
folder = '/home/smrak/Documents/TheMahali/imf/'

#def openIMF():
fname = 'omniweb20151007.csv'
Y, DinY, H, M, Bx, By, Bz, AE = np.loadtxt(folder+fname,dtype='S',
                                           unpack=True,delimiter=',')
Bx = Bx.astype(float)
By = By.astype(float)
Bz = Bz.astype(float)
AE = AE.astype(int)
H = H.astype(int)
M = M.astype(int)
t = []
for i in range(len(AE)):
    t.append(datetime.datetime(2015,10,7, H[i], M[i], 0))
t = np.array(t)

start = datetime.datetime(2015,10,7, 6, 0, 0)
stop = datetime.datetime(2015,10,7, 6, 40, 0)
idx = np.where( (t>=start) & (t<=stop) )[0]

obstimes = t[idx]

plotting.plotimf(t, Bx, By, Bz, AE/1000, xlabel='UT', ylabel1='B [nT]', 
                 ylabel2='AE Index [uT]', ylim1=[-20,20], ylim2=[0,1.8], 
                 ytick1=[-20,-15,-10,-5,0,5,10,15,20], ytick2=[0,0.4,0.8,1.2,1.6],
                 legend=True, obstimes=obstimes, centerline=True,
                 xlim=[datetime.datetime(2015,10,7,2,0,0),
                      datetime.datetime(2015,10,7,11,0,0)])

#mahaliPlot.plotImfField(np.array(t), Bx, By, Bz, AE)