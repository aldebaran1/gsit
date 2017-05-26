#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 19:42:25 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

def plotTEC(t, y, xlabel=None, ylabel=None, title=None, xlim=None, ylim=None,
            color='b', ytick=None, xtick=None):
    """
    Sebastijan Mrak
    Plot a single track t-y plot which is meant to be a time dependent graph. Time
    has to be in a datetime.datetime format. It offeres some optional parameters to
    shape the figure.
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(t, y, color)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xtick is not None:
        ax.set_xticks(xtick)
    if ytick is not None:
        ax.set_yticks(ytick)
    ax.xaxis.set(major_formatter=formatter)
    plt.show()
    
def plotTECLossOfLock(t, y, l, lli='x', xlabel=None, ylabel=None, title=None, xlim=None, ylim=None,
            color='b', colorx='xr', ms=10):
    """
    Sebastijan Mrak
    Plot a single track t-y plot which is meant to be a time dependent graph including
    a loss of lock indicators on the plot. Time has to be in a datetime.datetime 
    format. It offeres some optional parameters to shape the figure.
    """
    idx = np.where((l%2)==1)[0]
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(t, y, color)
    if lli == 'x':
        plt.plot(t[idx], y[idx], colorx, ms=ms)
    elif lli == 'line':
        lli_range = ax.get_ylim()
        
        for ix in idx:
            plt.plot([t[ix], t[ix]], [lli_range[0], lli_range[1]], 'r')
    else:
        print ("Enter right parameter for loss of lock index presentation")
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    if xlim is not None:
        ax.set_xlim(xlim)
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.xaxis.set(major_formatter=formatter)
    plt.show()
    
def Plot2Subplot(t1, t2, y1, y2, l=None, lli1='x', lli2='x', xlabel=None, ylabel1=None,
                 ylabel2=None, title1=None, title2=None, xlim=None, ms=10,
                 ylim1=None, ylim2=None, color1='b', colorx='xr', color2='b'):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    plt.plot(t1, y1, color1)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if title1 is not None:
        plt.title(title1)
    if l is not None:
        idx=np.where((l%2)==1)[0]
        if lli1 == 'x':
            plt.plot(t1[idx], y1[idx], colorx, ms=ms)
        elif lli1 == 'line':
            lli_range = ax1.get_ylim()
            
            for ix in idx:
                plt.plot([t1[ix], t1[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    
    ax2=fig.add_subplot(212, sharex=ax1)
    plt.plot(t2, y2, color2)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if title2 is not None:
        plt.title(title2)
    if l is not None:
        idx=np.where((l%2)==1)[0]
        if lli2 == 'x':
            plt.plot(t2[idx], y2[idx], colorx, ms=ms)
        elif lli2 == 'line':
            lli_range = ax2.get_ylim()
            
            for ix in idx:
                plt.plot([t2[ix], t2[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    if ylim2 is not None:
        ax1.set_ylim(ylim2)
    
    if xlabel is not None:
        ax2.set_xlabel(xlabel)
    if xlim is not None:
        ax1.set_xlim(xlim)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()