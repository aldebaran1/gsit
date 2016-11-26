#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 19:42:25 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

def plotTEC(t, y, xlabel=None, ylabel=None, title=None, xlim=None, ylim=None,
            color='b'):
    """
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
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.xaxis.set(major_formatter=formatter)
    plt.show()