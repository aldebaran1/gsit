#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 14:21:02 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

#%% Keograms

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
                 legend=None, ylim=None, pcolorbar=None, ytick1=None, ytick2=None,
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
    plt.rc('axes', labelsize=20)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1 = fig.add_subplot(211)
    if cmap is not None:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T), cmap=cmap)
    else:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T))
    plt.setp(ax1.get_xticklabels(), visible=False) 
    if title1 is not None:
        plt.title(title1)
    if legend is not None:
        plt.legend()

    ax2 = fig.add_subplot(212, sharex=ax1)
    if cmap is not None:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T), cmap=cmap)
    else:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    plt.xlabel('UT')
    if title2 is not None:
        plt.title(title2)
    if legend is not None:
        plt.legend()

    if ylim is not None:
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
        
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
        ax2.set_yticks(ytick2)
        
    if xtick is not None:
        ax1.set_xticks(xtick)
    if xlim is not None:
        ax1.set_xlim(xlim)
    
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()
    
    p1 = ax1.get_position()
    pos1 = p1.get_points()
    pos1[0][0] = 0.1
    pos1[1][0] = 0.86
    p1.set_points(pos1)
    ax1.set_position(p1)

    p2 = ax2.get_position()
    pos2 = p2.get_points()
    pos2[0][0] = 0.1
    pos2[1][0] = 0.86
    p2.set_points(pos2)
    ax2.set_position(p2)
    
    if pcolorbar is not None:
        if cbartick1 is not None:
            p1 = ax1.get_position()
            pos1 = p1.get_points()
            cbaxes = fig.add_axes([0.88, pos1[0][1]+0.01, 0.01, pos1[1][1]-pos1[0][1]-0.01])
            cbar = plt.colorbar(cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle1)
        if cbartick2 is not None:
            p2 = ax2.get_position()
            pos2 = p2.get_points()
            cbaxes = fig.add_axes([0.88, pos2[0][1]+0.01, 0.01, pos2[1][1]-pos2[0][1]-0.01])
            cbar = plt.colorbar(ticks=cbartick2, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle2)
    
    fig.text(0.04, 0.5, 'Elevation [deg]', va='center', rotation='vertical', fontsize=20)
    
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
    
def plot4Keogram(t1, t2, t3, t4, y1, y2, y3, y4, kg1, kg2, kg3, kg4, 
                 title1=None, title2=None, title3=None, title4=None, 
                 legend=None, ylim=None, pcolorbar=None, ytick1=None, 
                 ytick2=None, ytick3=None, ytick4=None, xtick=None, cmap=None,
                 cbartick1=None, cbartick2=None, cbartick3=None,cbartick4=None,
                 cbartitle1=None, cbartitle2=None, cbartitle3=None, cbartitle4=None,
                 xlim=None, elevation=None, obstimes=None, lli=None):
    """
    Sebastijan Mrak
    Function plot2Keogram takes x and y grid values, where x axis is meant to be
    time in datatime.dateimte format. It plots keogram values 'kg' on top of set
    grid with pcolormesh function. It takes 3 different sets of input data and plots
    them on separate subplots, where x-axis is shared among subplots.There are 
    also many supporting parameters to enrich the figure
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(8,6), dpi=150)
    
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=8)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)  # fontsize of the tick labels
    ax1 = fig.add_subplot(411)
    if cmap is not None:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T), cmap=cmap)
    else:
        plt.pcolormesh(t1, y1, np.nan_to_num(kg1.T))
    if title1 is not None:
        plt.title(title1)
    if legend is not None:
        plt.legend()
    plt.setp(ax1.get_xticklabels(), visible=False) 
################################################################################
    ax2 = fig.add_subplot(412, sharex=ax1)
    plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    if cmap is not None:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T), cmap=cmap)
    else:
        plt.pcolormesh(t2, y2, np.nan_to_num(kg2.T))
    plt.setp(ax2.get_xticklabels(), visible=False) 
    if title2 is not None:
        plt.title(title2)
    if legend is not None:
        plt.legend()
################################################################################
    ax3 = fig.add_subplot(413, sharex=ax1)
    if cmap is not None:
        plt.pcolormesh(t3, y3, np.nan_to_num(kg3.T), cmap=cmap)
    else:
        plt.pcolormesh(t3, y3, np.nan_to_num(kg3.T))
    if title3 is not None:
        plt.title(title3)
    if legend is not None:
        plt.legend()
    plt.setp(ax3.get_xticklabels(), visible=False) 
################################################################################
    ax4 = fig.add_subplot(414, sharex=ax1)
    if cmap is not None:
        plt.pcolormesh(t4, y4, np.nan_to_num(kg4.T), cmap=cmap)
    else:
        plt.pcolormesh(t4, y4, np.nan_to_num(kg4.T))
#    if elevation is not None:
#        plt.plot(obstimes, elevation, 'r', lw=1)
#    if lli is not None:
#        idx = np.where((lli%2) == 1)[0]
#        plt.scatter(obstimes[idx], elevation[idx], facecolors='none', edgecolors='m', s=10)
    plt.xlabel('UT')
    if title4 is not None:
        plt.title(title4)
    if legend is not None:
        plt.legend()
################################################################################
    if ylim is not None:
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
        ax4.set_ylim(ylim)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
    if ytick2 is not None:
        ax2.set_yticks(ytick2)
    if ytick3 is not None:
        ax3.set_yticks(ytick3)
    if ytick4 is not None:
        ax4.set_yticks(ytick4)
    if xtick is not None:
        ax1.set_xticks(xtick)
    if xlim is not None:
        ax1.set_xlim(xlim)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    
    fig.subplots_adjust(hspace = .01)
    plt.show()
    ############################################################################
    p1 = ax1.get_position()
    pos1 = p1.get_points()
    pos1[0][0] = 0.1
    pos1[1][0] = 0.86
    p1.set_points(pos1)
    ax1.set_position(p1)

    p2 = ax2.get_position()
    pos2 = p2.get_points()
    pos2[0][0] = 0.1
    pos2[1][0] = 0.86
    p2.set_points(pos2)
    ax2.set_position(p2)

    p3 = ax3.get_position()
    pos3 = p3.get_points()
    pos3[0][0] = 0.1
    pos3[1][0] = 0.86
    p3.set_points(pos3)
    ax3.set_position(p3)
    
    p4 = ax4.get_position()
    pos4 = p4.get_points()
    pos4[0][0] = 0.1
    pos4[1][0] = 0.86
    p4.set_points(pos4)
    ax4.set_position(p4)
    if pcolorbar is not None:
        if cbartick1 is not None:
            p1 = ax1.get_position()
            pos1 = p1.get_points()
            cbaxes = fig.add_axes([0.88, pos1[0][1]+0.01, 0.01, pos1[1][1]-pos1[0][1]-0.01])
            cbar = plt.colorbar(cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle1)
        if cbartick2 is not None:
            p2 = ax2.get_position()
            pos2 = p2.get_points()
            cbaxes = fig.add_axes([0.88, pos2[0][1]+0.01, 0.01, pos2[1][1]-pos2[0][1]-0.01])
            cbar = plt.colorbar(ticks=cbartick2, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle2)
        if cbartick3 is not None:
            p3 = ax3.get_position()
            pos3 = p3.get_points()
            cbaxes = fig.add_axes([0.88, pos3[0][1]+0.01, 0.01, pos3[1][1]-pos3[0][1]-0.01])
            cbar = plt.colorbar(ticks=cbartick3, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle3)
        if cbartick4 is not None:
            p4 = ax4.get_position()
            pos4 = p4.get_points()
            cbaxes = fig.add_axes([0.88, pos4[0][1]+0.01, 0.01, pos4[1][1]-pos4[0][1]-0.01])
            cbar = plt.colorbar(ticks=cbartick4, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle4)
        else:
            plt.colorbar()
            
    if elevation is not None:
        ax1.plot(obstimes, elevation, 'r', lw=1)
        ax4.plot(obstimes, elevation, 'r', lw=1)
    if lli is not None:
        idx = np.where((lli%2) == 1)[0]
        ax1.scatter(obstimes[idx], elevation[idx], facecolors='none', edgecolors='m', s=10)
        ax4.scatter(obstimes[idx], elevation[idx], facecolors='none', edgecolors='m', s=10)
    fig.text(0.04, 0.5, 'Elevation [deg]', va='center', rotation='vertical')
    
#%% GPS stuff

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
    
def plot2SubplotGPS(t1, t2, y1, y2, l=None, lli1='x', lli2='x', xlabel=None, ylabel1=None,
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
    
#%% XY Plots

def plot2subplotK_xy(x, y, i, t, B, el=None, xlim=None, ylim1=None, ylim2=None,
                     xlabel=None, ylabel1=None, ylabel2=None, ytick1=None, ytick2=None, 
                     cmap='viridis', pcolorbar=None, cbartick=None, cbartitle=None,
                     obstimes=None, ipp=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    plt.rc('axes', labelsize=20)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1=fig.add_subplot(211)
    plt.pcolormesh(x, y, np.nan_to_num(i.T), cmap=cmap)
    if ipp is not None:
        plt.plot(obstimes, ipp, '-r', lw=1)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    
    ax2 = fig.add_subplot(212, sharex=ax1)
    plt.plot(t, B, lw=2)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
        
    if xlim is not None:
        ax1.set_xlim(xlim)
    if xlabel is not None:
        plt.xlabel(xlabel)
        
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()
    
    #
    p1 = ax1.get_position()
    pos1 = p1.get_points()
    pos1[1][0] = 0.86
    p1.set_points(pos1)
    ax1.set_position(p1)
    # find current position [x,y,width,height]
    p1 = ax1.get_position()
    p2 = ax2.get_position()
    pos1 = p1.get_points()
    pos2 = p2.get_points()
    # set width of second axes equal to first
    pos2[1][0] = pos1[1][0]
    p2.set_points(pos2)
    ax2.set_position(p2)
    
    if pcolorbar is not None:
        if cbartick is not None:
            p1 = ax1.get_position()
            pos1 = p1.get_points()
            cbaxes = fig.add_axes([0.88, pos1[0][1], 0.01, pos1[1][1]-pos1[0][1]])
            cbar = plt.colorbar(ticks=cbartick, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle)
        else:
            plt.colorbar()

def plot2subplot(t1, t2, y1, y2, t21=None, t22=None, y21=None, y22=None,
                 color1='b', color2='b', color3='r', color4='g',
                 xlabel=None, title1=None, title2=None, ylabel1=None, 
                 ylabel2=None, xlim=None, ylim1=None, ylim2=None, 
                 label1=None, label2=None, label3=None, label4=None,
                 legend1=None, legend2=None, lli1=None, lm='x',
                 ms=10, colorx1='xr', lli2=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    if legend1 is not None:
        plt.plot(t1, y1, color1, label=label1)
        plt.legend()
    else:
        plt.plot(t1, y1, color1)
    if lli1 is not None:
        idx=np.where((lli1%2)==1)[0]
        if lm == 'x':
            plt.plot(t1[idx], y1[idx], colorx1, ms=ms)
        elif lm == 'line':
            lli_range = ax1.get_ylim()
            for ix in idx:
                plt.plot([t1[ix], t1[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    plt.setp(ax1.get_xticklabels(), visible=False) 
    if title1 is not None:
        plt.title(title1)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
        
    ax2 = fig.add_subplot(212, sharex=ax1)
    if legend2 is None:
        plt.plot(t2, y2, color2)
        if (t21 is not None) and (y21 is not None):
            plt.plot(t21, y21, color3)
        if (t22 is not None) and (y22 is not None):
            plt.plot(t22, y22, color4)
    else:
        plt.plot(t2, y2, color2, label=label2)
        if (t21 is not None) and (y21 is not None):
            plt.plot(t21, y21, color3, label=label3)
        if (t22 is not None) and (y22 is not None):
            plt.plot(t22, y22, color4, label=label4)
        ax2.legend(loc=2, bbox_to_anchor=(1.001, 0,0, 1), prop={'size':10}, 
                   fancybox=True)
    if lli2 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax2.get_ylim()
        for ix in idx:
            plt.plot([t2[ix], t2[ix]], [lli_range[0], lli_range[1]], 'r')
    if title2 is not None:
        plt.title(title2)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
        
    if xlim is not None:
        ax1.set_xlim(xlim)
    if xlabel is not None:
        plt.xlabel(xlabel)
        
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()
    
def plot3subplot(t1, t2, t3, y1, y2, y3, el1=None, t31=None, t32=None, y31=None, y32=None,
                 color1='b', color2='b', color3='r', color4='g',
                 xlabel=None, title1=None, title2=None, title3=None, ylabel1=None, 
                 ylabel2=None, ylabel3=None, xlim=None, ylim1=None, ylim2=None, ylim3=None,
                 label1=None, label2=None, label3=None, label4=None, ytick1=[],
                 legend1=None, legend2=None, lli1=None, lm='x',
                 ms=10, colorx1='xr', lli2=None, cmap='viridis', pcolorbar=None, 
                 cbartick=None, cbartitle=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(800,600), dpi=150)
    ax1=fig.add_subplot(311)
    plt.pcolormesh(t1, el1, np.nan_to_num(y1.T), cmap=cmap)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if title1 is not None:
        plt.title(title1)
    if pcolorbar is not None:
        if cbartick is not None:
            cbar = plt.colorbar(ticks=cbartick)
            cbar.ax.set_ylabel(cbartitle)
        else:
            plt.colorbar()
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
        
    ax2 = fig.add_subplot(312, sharex=ax1)
    if legend1 is not None:
        plt.plot(t2, y2, color1, label=label1)
        plt.legend()
    else:
        plt.plot(t2, y2, color1)
    if lli1 is not None:
        idx=np.where((lli1%2)==1)[0]
        if lm == 'x':
            plt.plot(t2[idx], y2[idx], colorx1, ms=ms)
        elif lm == 'line':
            lli_range = ax2.get_ylim()
            for ix in idx:
                plt.plot([t2[ix], t2[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    plt.setp(ax1.get_xticklabels(), visible=False) 
    if title2 is not None:
        plt.title(title1)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax3 = fig.add_subplot(313, sharex=ax1)
    if legend2 is None:
        plt.plot(t3, y3, color2)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4)
    else:
        plt.plot(t3, y3, color2, label=label2)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3, label=label3)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4, label=label4)
        ax3.legend(loc=2, bbox_to_anchor=(1.001, 0,0, 1), prop={'size':10}, 
                   fancybox=True)
    if lli2 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax3.get_ylim()
        for ix in idx:
            plt.plot([t1[ix], t1[ix]], [lli_range[0], lli_range[1]], 'r')
    if title3 is not None:
        plt.title(title3)
    if ylabel3 is not None:
        plt.ylabel(ylabel3)
    if ylim3 is not None:
        ax2.set_ylim(ylim3)
        
    if xlim is not None:
        ax1.set_xlim(xlim)
    if xlabel is not None:
        plt.xlabel(xlabel)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    plt.show()
        
    # find current position [x,y,width,height]
    p1 = ax1.get_position()
    p2 = ax2.get_position()
    p3 = ax3.get_position()
    pos1 = p1.get_points()
    pos2 = p2.get_points()
    pos3 = p3.get_points()
    # set width of second axes equal to first
    pos2[1][0] = pos1[1][0]
    p2.set_points(pos2)
    ax2.set_position(p2)
    pos3[1][0] = pos1[1][0]
    p3.set_points(pos3)
    ax3.set_position(p3)    
    
def plot4subplot(t1, t2, t3, t4, y1, y2, y3,y4, el1=None, t31=None, t32=None, 
                 y31=None, y32=None, color1='b', color2='b', color3='r', color4='g',
                 color5='b', xlabel=None, title1=None, title2=None, title3=None, 
                 title4=None, ylabel1=None, ylabel2=None, ylabel3=None, ylabel4=None,
                 xlim=None, ylim1=None, ylim2=None, ylim3=None, ylim4=None,
                 label1=None, label2=None, label3=None, label4=None, ytick1=[],
                 legend1=None, legend2=None, lli1=None, lm='x', lli3=None,
                 ms=10, colorx1='xr', lli2=None, cmap='viridis', pcolorbar=None, 
                 cbartick=None, cbartitle=None, ytick2=None, ytick3=None, ytick4=None,
                 ipp_elevation=None, lli_keo=None, lw=2):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    plt.rc('axes', labelsize=18)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1=fig.add_subplot(411)
    plt.pcolormesh(t1, el1, np.nan_to_num(y1.T), cmap=cmap)
    if ipp_elevation is not None:
        plt.plot(t2, ipp_elevation, '-r', lw=lw)
    if lli_keo is not None:
        idx = np.where((lli_keo%2) == 1)[0]
        ax1.scatter(t2[idx], ipp_elevation[idx], facecolors='none', edgecolors='m', s=15)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if title1 is not None:
        plt.title(title1)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
     
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ############################################################################
    ax2 = fig.add_subplot(412, sharex=ax1)
    if legend1 is not None:
        plt.plot(t2, y2, color1, label=label1)
        plt.legend()
    else:
        plt.plot(t2, y2, color1)
    if lli1 is not None:
        idx=np.where((lli1%2)==1)[0]
        if lm == 'x':
            plt.plot(t2[idx], y2[idx], colorx1, ms=ms)
        elif lm == 'line':
            lli_range = ax2.get_ylim()
            for ix in idx:
                plt.plot([t2[ix], t2[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    if title2 is not None:
        plt.title(title1)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    if ytick2 is not None:
        ax2.set_yticks(ytick2)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ############################################################################
    ax3 = fig.add_subplot(413, sharex=ax1)
    if legend2 is None:
        plt.plot(t3, y3, color2)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4)
    else:
        plt.plot(t3, y3, color2, label=label2)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3, label=label3)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4, label=label4)
        ax3.legend(loc=2, bbox_to_anchor=(0.88, 0,0, 1), prop={'size':12}, 
                   fancybox=True)
    if lli2 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax3.get_ylim()
        for ix in idx:
            plt.plot([t1[ix], t1[ix]], [lli_range[0], lli_range[1]], 'r')
    if title3 is not None:
        plt.title(title3)
    if ylabel3 is not None:
        plt.ylabel(ylabel3)
    if ylim3 is not None:
        ax3.set_ylim(ylim3)
    if ytick3 is not None:
        ax3.set_yticks(ytick3)
    plt.setp(ax3.get_xticklabels(), visible=False) 
    ############################################################################
    ax4 = fig.add_subplot(414, sharex=ax1)
    plt.plot(t4, y4, color5)
    if lli3 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax4.get_ylim()
        for ix in idx:
            plt.plot([t4[ix], t4[ix]], [lli_range[0], lli_range[1]], 'r')
    if title4 is not None:
        plt.title(title4)
    if ylabel4 is not None:
        plt.ylabel(ylabel4)
    if ylim4 is not None:
        ax4.set_ylim(ylim4)
    if ytick4 is not None:
        ax4.set_yticks(ytick4)
    ############################################################################
    if xlim is not None:
        ax1.set_xlim(xlim)
    if xlabel is not None:
        plt.xlabel(xlabel)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()
    
    #
    p1 = ax1.get_position()
    pos1 = p1.get_points()
    pos1[1][0] = 0.86
    p1.set_points(pos1)
    ax1.set_position(p1)
    # find current position [x,y,width,height]
    p1 = ax1.get_position()
    p2 = ax2.get_position()
    p3 = ax3.get_position()
    p4 = ax4.get_position()
    pos1 = p1.get_points()
    pos2 = p2.get_points()
    pos3 = p3.get_points()
    pos4 = p4.get_points()
    # set width of second axes equal to first
    pos2[1][0] = pos1[1][0]
    p2.set_points(pos2)
    ax2.set_position(p2)
    pos3[1][0] = pos1[1][0]
    p3.set_points(pos3)
    ax3.set_position(p3)
    pos4[1][0] = pos1[1][0]
    p4.set_points(pos4)
    ax4.set_position(p4)
    
    if pcolorbar is not None:
        if cbartick is not None:
            p1 = ax1.get_position()
            pos1 = p1.get_points()
            cbaxes = fig.add_axes([0.88, pos1[0][1], 0.01, pos1[1][1]-pos1[0][1]])
            cbar = plt.colorbar(ticks=cbartick, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle)
        else:
            plt.colorbar()
            
def plot5subplot(t1, t2, t3, t4, y1, y2, y3,y4, el1=None, t31=None, t32=None, 
                 y31=None, y32=None, color1='b', color2='b', color3='r', color4='g',
                 color5='b', xlabel=None, title1=None, title2=None, title3=None, 
                 title4=None, ylabel1=None, ylabel2=None, ylabel3=None, ylabel4=None,
                 xlim=None, ylim1=None, ylim2=None, ylim3=None, ylim4=None,
                 label1=None, label2=None, label3=None, label4=None, ytick1=[],
                 legend1=None, legend2=None, lli1=None, lm='x', lli3=None,
                 ms=10, colorx1='xr', lli2=None, cmap='viridis', pcolorbar=None, 
                 cbartick=None, cbartitle=None, ytick2=None, ytick3=None, ytick4=None,
                 ipp_elevation=None, lli_keo=None, lw=2, Bt=None, Bx=None, title5=None,
                 cbx=None, ylabel5=None, ylim5=None, ytick5=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    plt.rc('axes', labelsize=16)
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    ax1=fig.add_subplot(511)
    plt.pcolormesh(t1, el1, np.nan_to_num(y1.T), cmap=cmap)
    if ipp_elevation is not None:
        plt.plot(t3, ipp_elevation, '-r', lw=lw)
    if lli_keo is not None:
        idx = np.where((lli_keo%2) == 1)[0]
        ax1.scatter(t3[idx], ipp_elevation[idx], facecolors='none', edgecolors='m', s=15)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if title1 is not None:
        plt.title(title1)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
     
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ############################################################################
    ax2 = fig.add_subplot(512, sharex=ax1)
    if legend1 is None:
        plt.plot(t2, y2, color2)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4)
    else:
        plt.plot(t2, y2, color2,label=label2)
        plt.plot(t2, y2, color2+'.',)
        if (t31 is not None) and (y31 is not None):
            plt.plot(t31, y31, color3, label=label3)
            plt.plot(t31, y31, color3+'.',)
        if (t32 is not None) and (y32 is not None):
            plt.plot(t32, y32, color4, label=label4)
            plt.plot(t32, y32, color4+'.',)
            ax2.legend(loc=2, bbox_to_anchor=(0.88, 0,0, 1), prop={'size':12}, 
                   fancybox=True)
    if lli2 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax2.get_ylim()
        for ix in idx:
            plt.plot([t1[ix], t1[ix]], [lli_range[0], lli_range[1]], 'r')
    if title2 is not None:
        plt.title(title2)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    if ytick2 is not None:
        ax2.set_yticks(ytick2)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ############################################################################
    ax3 = fig.add_subplot(513, sharex=ax1)
    if legend2 is not None:
        plt.plot(t3, y3, color1, lw=2, label=label1)
        plt.legend()
    else:
        plt.plot(t3, y3, color1, lw=2)
    if lli1 is not None:
        idx=np.where((lli1%2)==1)[0]
        if lm == 'x':
            plt.plot(t3[idx], y3[idx], colorx1, ms=ms)
        elif lm == 'line':
            lli_range = ax3.get_ylim()
            for ix in idx:
                plt.plot([t3[ix], t3[ix]], [lli_range[0], lli_range[1]], 'r')
        else:
            print ("Enter right parameter for loss of lock index presentation")
    if title3 is not None:
        plt.title(title3)
    if ylabel3 is not None:
        plt.ylabel(ylabel3)
    if ylim3 is not None:
        ax3.set_ylim(ylim3)
    if ytick3 is not None:
        ax3.set_yticks(ytick3)
    plt.setp(ax3.get_xticklabels(), visible=False) 
    ############################################################################
    ax4 = fig.add_subplot(514, sharex=ax1)
    plt.plot(t4, y4, color5, lw=2)
    if lli3 is not None:
        idx=np.where((lli2%2)==1)[0]
        lli_range = ax4.get_ylim()
        for ix in idx:
            plt.plot([t4[ix], t4[ix]], [lli_range[0], lli_range[1]], 'r')
    if title4 is not None:
        plt.title(title4)
    if ylabel4 is not None:
        plt.ylabel(ylabel4)
    if ylim4 is not None:
        ax4.set_ylim(ylim4)
    if ytick4 is not None:
        ax4.set_yticks(ytick4)
    plt.setp(ax4.get_xticklabels(), visible=False) 
    ############################################################################
    ax5 = fig.add_subplot(515, sharex=ax1)
    plt.plot(Bt, Bx, color=cbx, lw=2)
    if title5 is not None:
        plt.title(title5)
    if ylabel5 is not None:
        plt.ylabel(ylabel5)
    if ylim5 is not None:
        ax5.set_ylim(ylim5)
    if ytick5 is not None:
        ax5.set_yticks(ytick5)
    ############################################################################
    if xlim is not None:
        ax1.set_xlim(xlim)
    if xlabel is not None:
        plt.xlabel(xlabel)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()
    
    #
    p1 = ax1.get_position()
    pos1 = p1.get_points()
    pos1[1][0] = 0.86
    p1.set_points(pos1)
    ax1.set_position(p1)
    # find current position [x,y,width,height]
    p1 = ax1.get_position()
    p2 = ax2.get_position()
    p3 = ax3.get_position()
    p4 = ax4.get_position()
    p5 = ax5.get_position()
    pos1 = p1.get_points()
    pos2 = p2.get_points()
    pos3 = p3.get_points()
    pos4 = p4.get_points()
    pos5 = p5.get_points()
    # set width of second axes equal to first
    pos2[1][0] = pos1[1][0]
    p2.set_points(pos2)
    ax2.set_position(p2)
    pos3[1][0] = pos1[1][0]
    p3.set_points(pos3)
    ax3.set_position(p3)
    pos4[1][0] = pos1[1][0]
    p4.set_points(pos4)
    ax4.set_position(p4)
    pos5[1][0] = pos1[1][0]
    p5.set_points(pos5)
    ax5.set_position(p5)
    
    if pcolorbar is not None:
        if cbartick is not None:
            p1 = ax1.get_position()
            pos1 = p1.get_points()
            cbaxes = fig.add_axes([0.88, pos1[0][1], 0.01, pos1[1][1]-pos1[0][1]])
            cbar = plt.colorbar(ticks=cbartick, cax = cbaxes)
            cbar.ax.set_ylabel(cbartitle)
        else:
            plt.colorbar()
            
def plotimf(t, Bx, By, Bz, AE, xlabel=None, ylabel1=None, ylabel2=None, xlim=None,
           ylim1=None, ylim2=None, legend=None, lw=1, ytick1=None, ytick2=None, 
           xtick=None, obstimes=None, centerline=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(15,5))
    plt.rc('axes', labelsize=18)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1=fig.add_subplot(211)
    
    plt.plot(t, Bx, 'b', lw=lw, label='Bx')
    plt.plot(t, By, 'm', lw=lw, label='By')
    plt.plot(t, Bz, 'r', lw=lw, label='Bz')
    if centerline is True:
        plt.plot([t[0], t[-1]], [0,0], 'k', lw=lw)
    if obstimes is not None:
        plt.plot
    if legend is not None:
        ax1.legend(loc=2, bbox_to_anchor=(1.01, 0,0, 1), prop={'size':12}, 
                   fancybox=True)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
    if obstimes is not None:
        yrange = ax1.get_ylim()
        plt.plot([obstimes[0],obstimes[0]], [yrange[0], yrange[1]], '--k', lw=1)
        plt.plot([obstimes[-1],obstimes[-1]], [yrange[0], yrange[1]], '--k', lw=1)
        
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ax2=fig.add_subplot(212, sharex=ax1)
    plt.plot(t, AE, color='b', lw=lw)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    if ytick2 is not None:
        ax2.set_yticks(ytick2)
    if obstimes is not None:
        yrange = ax2.get_ylim()
        plt.plot([obstimes[0],obstimes[0]], [yrange[0], yrange[1]], '--k', lw=1)
        plt.plot([obstimes[-1],obstimes[-1]], [yrange[0], yrange[1]], '--k', lw=1)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if xlim is not None:
        ax1.set_xlim(xlim)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()

def plotimfExtended(t, Bx, By, Bz, AE=None, v=None, xlabel=None, ylabel1=None, 
           ylabel2=None, ylabel3=None, xlim=None,
           ylim1=None, ylim2=None, ylim3=None, legend=None, lw=1, ytick1=None, 
           ytick2=None, ytick3=None, xtick=None, obstimes=None, centerline=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    plt.rc('axes', labelsize=18)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1=fig.add_subplot(311)
    
    plt.plot(t, Bx, 'b', lw=lw, label='Bx')
    plt.plot(t, By, 'm', lw=lw, label='By')
    plt.plot(t, Bz, 'r', lw=lw, label='Bz')
    if centerline is True:
        plt.plot([t[0], t[-1]], [0,0], 'k', lw=lw)
    if obstimes is not None:
        plt.plot
    if legend is not None:
        ax1.legend(loc=2, bbox_to_anchor=(1.01, 0,0, 1), prop={'size':12}, 
                   fancybox=True)
    if ylabel1 is not None:
        plt.ylabel(ylabel1)
    if ylim1 is not None:
        ax1.set_ylim(ylim1)
    if ytick1 is not None:
        ax1.set_yticks(ytick1)
    if obstimes is not None:
        yrange = ax1.get_ylim()
        plt.plot([obstimes[0],obstimes[0]], [yrange[0], yrange[1]], '--k', lw=1)
        plt.plot([obstimes[-1],obstimes[-1]], [yrange[0], yrange[1]], '--k', lw=1)
        
    plt.setp(ax1.get_xticklabels(), visible=False) 
    ax2=fig.add_subplot(312, sharex=ax1)
    plt.plot(t, AE, color='b', lw=lw)
    if ylabel2 is not None:
        plt.ylabel(ylabel2)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    if ytick2 is not None:
        ax2.set_yticks(ytick2)
    if obstimes is not None:
        yrange = ax2.get_ylim()
        plt.plot([obstimes[0],obstimes[0]], [yrange[0], yrange[1]], '--k', lw=1)
        plt.plot([obstimes[-1],obstimes[-1]], [yrange[0], yrange[1]], '--k', lw=1)
    
    plt.setp(ax2.get_xticklabels(), visible=False) 
    ax3=fig.add_subplot(313, sharex=ax1)
    plt.plot(t, v, color='b', lw=lw)
    if ylabel3 is not None:
        plt.ylabel(ylabel3)
    if ylim3 is not None:
        ax3.set_ylim(ylim3)
    if ytick3 is not None:
        ax3.set_yticks(ytick3)
        
    if xlabel is not None:
        plt.xlabel(xlabel)
    if xlim is not None:
        ax1.set_xlim(xlim)
    ax1.xaxis.set(major_formatter=formatter)
    fig.tight_layout()
    fig.subplots_adjust(hspace = .01)
    plt.show()
    
def plotMagnetometer(t, y1, y2, ylabel=None, xlabel='UT', xlim=None, ylim=None, 
                     colorx='b', colory='k', colorz='m', zeroline=False):
    """
    Sebastijan Mrak
    Magnetometer, all components
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    plt.rc('axes', labelsize=18)
    plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    ax1=fig.add_subplot(111)
    ax1.plot(t, y1[1][:], colorx)
    ax1.plot(t, y2[0][:], colory)
    
    if xlim is not None:
        ax1.set_xlim(xlim)
    
    ax1.xaxis.set(major_formatter=formatter)
    ax1.set_xlabel(xlabel)
    
    plt.show()
    