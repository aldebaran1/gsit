#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:27:15 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

def plotKeogram(t, y, kg, title=None, legend=None, ylim=None, pcolorbar=None, 
                ytick=None, cmap=None, cbartick=None, cbartitle=None,
                xtick=None, xlim=None):
    """
    Sebastijan Mrak
    Function plotKeogram takes x and y grid values, where x axis is meant to be
    time in datatime.dateimte format. It plots keogram values 'kg' on top of set
    grid with pcolormesh function. There are also many supporting parameters to 
    enrich the figure.
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if cmap is not None:
        plt.pcolormesh(t, y, np.nan_to_num(kg.T), cmap=cmap)
    else:
        plt.pcolormesh(t, y, np.nan_to_num(kg.T))
    plt.xlabel('UT')
    plt.ylabel('Elevation [deg]')
    if title is not None:
        plt.title(title)
    if ylim is not None:
        ax.set_ylim(ylim)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick is not None:
            cbar = plt.colorbar(ticks=cbartick)
            cbar.ax.set_ylabel(cbartitle)
        else:
            plt.colorbar()
    if ytick is not None:
        ax.set_yticks(ytick)
    if xtick is not None:
        ax.set_xticks(xtick)
    if xlim is not None:
        ax.set_xlim(xlim)
        
    ax.xaxis.set(major_formatter=formatter)
    plt.show()
    
def plot2Keogram(t1, t2, y1, y2, kg1, kg2, title1=None, title2=None,
                 legend=None, ylim=None, pcolorbar=None, ytick=None, 
                 cmap=None, cbartick1=None, cbartick2=None, 
                 cbartitle1=None, cbartitle2=None, xtick=None, xlim=None):
    """
    Sebastijan Mrak
    Function plot2Keogram takes x and y grid values, where x axis is meant to be
    time in datatime.dateimte format. It plots keogram values 'kg' on top of set
    grid with pcolormesh function. It takes 2 different sets of input data and plots
    them on separate subplots, where x-axis is shared among subplots.There are 
    also many supporting parameters to enrich the figure.
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    if cmap is not None:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T), cmap=cmap)
    else:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T))
    plt.setp(ax1.get_xticklabels(), visible=False) 
    plt.ylabel('Elevation [deg]')
    if title1 is not None:
        plt.title(title1)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick1 is not None:
            cbar = plt.colorbar(ticks=cbartick1)
            cbar.ax.set_ylabel(cbartitle1)
        else:
            plt.colorbar()
        
    ax2 = fig.add_subplot(212, sharex=ax1)
    if cmap is not None:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T), cmap=cmap)
    else:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    plt.xlabel('UT')
    plt.ylabel('Elevation [deg]')
    if title2 is not None:
        plt.title(title2)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick2 is not None:
            cbar = plt.colorbar(ticks=cbartick2)
            cbar.ax.set_ylabel(cbartitle2)
        else:
            plt.colorbar()
        
    if ylim is not None:
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
        
    if ytick is not None:
        ax1.set_yticks(ytick)
        ax2.set_yticks(ytick)
        
    if xtick is not None:
        ax1.set_xticks(xtick)
    if xlim is not None:
        ax1.set_xlim(xlim)
    
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()
    
def plot3Keogram(t1, t2, t3, y1, y2, y3, kg1, kg2, kg3, title1=None, 
                 title2=None, title3=None, legend=None, ylim=None, 
                 pcolorbar=None, ytick=None, xtick=None, cmap=None,
                 cbartick1=None, cbartick2=None, cbartick3=None, 
                 cbartitle1=None, cbartitle2=None, cbartitle3=None, 
                 xlim=None):
    """
    Sebastijan Mrak
    Function plot2Keogram takes x and y grid values, where x axis is meant to be
    time in datatime.dateimte format. It plots keogram values 'kg' on top of set
    grid with pcolormesh function. It takes 3 different sets of input data and plots
    them on separate subplots, where x-axis is shared among subplots.There are 
    also many supporting parameters to enrich the figure
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    if cmap is not None:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T), cmap=cmap)
    else:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T))
    plt.setp(ax1.get_xticklabels(), visible=False) 
    plt.ylabel('Elevation [deg]')
    if title1 is not None:
        plt.title(title1)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick1 is not None:
            cbar = plt.colorbar(ticks=cbartick1)
            cbar.ax.set_ylabel(cbartitle1)
        else:
            plt.colorbar()
        
    ax2 = fig.add_subplot(312, sharex=ax1)
    plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    if cmap is not None:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T), cmap=cmap)
    else:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    plt.setp(ax2.get_xticklabels(), visible=False) 
    plt.ylabel('Elevation [deg]')
    if title2 is not None:
        plt.title(title2)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick2 is not None:
            cbar2 = plt.colorbar(ticks=cbartick2)
            cbar2.ax.set_ylabel(cbartitle2)
        else:
            plt.colorbar()
        
    ax3 = fig.add_subplot(313, sharex=ax1)
    if cmap is not None:
        plt.pcolormesh(t3, y3, np.nan_to_num(kg3.T), cmap=cmap)
    else:
        plt.pcolormesh(t3, y3, np.nan_to_num(kg3.T))
    plt.xlabel('UT')
    plt.ylabel('Elevation [deg]')
    if title3 is not None:
        plt.title(title3)
    if legend is not None:
        plt.legend()
    if pcolorbar is not None:
        if cbartick2 is not None:
            cbar3 = plt.colorbar(ticks=cbartick3)
            cbar3.ax.set_ylabel(cbartitle3)
        else:
            plt.colorbar()
        
    if ylim is not None:
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
    if ytick is not None:
        ax1.set_yticks(ytick)
        ax2.set_yticks(ytick)
        ax3.set_yticks(ytick)
    if xtick is not None:
        ax1.set_xticks(xtick)
    if xlim is not None:
        ax1.set_xlim(xlim)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()
    
def plotIntensity(t, y, y2=None, y3=None, t2=None, t3=None,
                  ylim=None, xlim=None, ytick=None, xtick=None, 
                  title=None, xlabel=None, ylabel=None, label1=None, label2=None,
                  label3=None, legend=None, color1='b', color2='g', color3='r'):
    """
    Sebastijan Mrak
    function plotIntensity plots up to 3 datasets on signular figure. X-axis is meant
    to be a time variable in datatime.datatime format, y-axis is an intensity of the
    aurora. Function also offers additional paramters to enrich the figure.
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()    
    if (t2 is not None) and (y2 is not None) and (t3 is None):
        plt.plot(t, y, color1, label=label1)
        plt.plot(t2, y2, color2, label=label2)
    elif (t3 is not None) and (y3 is not None):
        plt.plot(t, y, color1, label=label1)
        plt.plot(t2, y2, color2, label=label2)
        plt.plot(t3, y3, color3, label=label3)
    else:
        plt.plot(t, y, color1, label=label1)
    if legend is not None:
        plt.legend()
    #if ylim is not None:
    axes = plt.gca()
    if ylim is not None:
        axes.set_ylim(ylim)
    if xlim is not None:
        axes.set_xlim(xlim)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if title is not None:
        plt.title(title)
    axes.xaxis.set(major_formatter=formatter)
    plt.show()
    
    
    