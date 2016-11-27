#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 14:21:02 2016

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

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
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),
                   ncol=3, fancybox=True, shadow=True)
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
    fig = plt.figure()
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
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),
                   ncol=3, fancybox=True, shadow=True)
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
                 cbartick=None, cbartitle=None, ytick2=None, ytick3=None, ytick4=None):
    """
    Sebastijan Mrak
    """
    formatter = mdates.DateFormatter('%H:%M')
    fig = plt.figure()
    ax1=fig.add_subplot(411)
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
        ax3.legend(loc=2, bbox_to_anchor=(1.001, 0,0, 1), prop={'size':10}, fancybox=True)
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