# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 09:57:11 2016

@author: smrak
"""

import numpy as np
from pandas import read_hdf
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal

#fs = 10
#order = 5
#highcut = 0.025
#nyq = 0.5 * fs
#high = highcut / nyq
#b, a = signal.butter(order, high, btype='highpass', analog=True)
#w, h = signal.freqs(b, a)
#plt.plot(w*nyq, 20 * np.log10(abs(h)))
#plt.ylim([-50, 1])
#plt.xscale('log')


def butter_hp(highcut, fs, order=3):
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high, btype='highpass', analog=False)
    return b, a

rx = 8
sv = 23
skip = 100
rinex = 'C:\\Users\\smrak\\Google Drive\\BU\\software\\gsit\\paper\\'
receiver = {2:'mah22800.h5', 3:'mah32800.h5', 4:'mah42800.h5', 5:'mah52800.h5',
            6:'mah62800.h5', 7:'mah72800.h5', 8:'mah82800.h5', 9:'mah92800.h5',
            13:'ma132800.h5'}
f1 = 1575.42E6
c0 = 3E8
fn = rinex + receiver.get(rx)
data = read_hdf(fn)
obstimes = data.major_axis
L1 = np.array(data['L1', sv, skip:, 'data'])*c0/f1
P1 = np.array(data['C1', sv, skip:, 'data'])
idx = np.where(np.isfinite(L1))[0]
L1 = L1[idx]
P1 = P1[idx]
y = L1-P1
x = np.arange(0, y.shape[0])
f = interpolate.interp1d(x, y)

x_new = np.arange(0, y.shape[0]-1, 0.01)
y_new = f(x_new)

b, a = butter_hp(0.005, 1)
Y1 = signal.lfilter(b, a, y)
Y2 = signal.lfilter(b, a, y_new)
plt.plot(Y2[500:]/c0*f1*2*np.pi)
N = 60
sp = []
Y2 = Y2[500:]
for i in range(len(Y2)-N):
    sp.append(np.std(Y2[i:i+N]))
#plt.plot(Y2)