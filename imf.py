# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 20:21:46 2016

@author: smrak
"""

import numpy as np
import mahaliPlot
import matplotlib.pyplot as plt
folder = '/home/smrak/Documents/TheMahali/imf/'

def openIMF():
    fname = 'space_20151007.csv'
    t, B, Bx, By, Bz, Kp, F, AE = np.loadtxt(folder+fname,dtype='S',unpack=True,delimiter=',')
    t = t.astype(int)    
    mahaliPlot.plotImfField(np.array(t), Bx, By, Bz, AE)
    
    
if __name__ == '__main__':
    openIMF()