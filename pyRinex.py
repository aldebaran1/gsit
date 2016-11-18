#!/usr/bin/env
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 14:15:16 2016
@author of the file: smrak
Autrors of the code: Greg Starr and Michael Hirsch
"""

from __future__ import division,absolute_import,print_function
import numpy as np
from datetime import datetime
from pandas import Panel4D, DataFrame, Series
from pandas.io.pytables import read_hdf
from os.path import splitext,expanduser
from io import BytesIO
from os.path import getsize

def readRinexNav(rinex_nav_filename):
    """
    Michael Hirsch
    It may actually be faster to read the entire file via f.read() and then .split()
    and asarray().reshape() to the final result, but I did it frame by frame.
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
    """

    startcol = 3 #column where numerical data starts
    nfloat = 19  #number of text elements per float data number
    nline = 7    #number of lines per record

    with open(expanduser(rinex_nav_filename),'r') as f:
        #find end of header, which has non-constant length
        while True:
            if 'END OF HEADER' in f.readline(): break
        #handle frame by frame
        sv = []; epoch=[]; raws=''
        while True:
            headln = f.readline()
            if not headln: break
            #handle the header
            sv.append(headln[:2])
            year = int(headln[2:5])
            if (80 <= year <= 99):
                year+=1900
            elif (year < 80): #good till year 2180
                year+=2000
                
            epoch.append(datetime(year = year,
                                  month       = int(headln[5:8]),
                                  day         = int(headln[8:11]),
                                  hour        = int(headln[11:14]),
                                  minute      = int(headln[14:17]),
                                  second      = int(headln[17:20]),
                                  microsecond = int(headln[21])*100000))
            """
            now get the data.
            Use rstrip() to chomp newlines consistently on Windows and Python 2.7/3.4
            Specifically [:-1] doesn't work consistently as .rstrip() does here.
            """
            raw = (headln[22:].rstrip() +
                   ''.join(f.readline()[startcol:].rstrip() for _ in range(nline-1))
                   +f.readline()[startcol:40].rstrip())

            raws += raw + '\n'

    raws = raws.replace('D','E')

    strio = BytesIO(raws.encode())
    darr = np.genfromtxt(strio,delimiter=nfloat)
    nav= DataFrame(darr, epoch,
               ['SVclockBias','SVclockDrift','SVclockDriftRate','IODE',
                'Crs','DeltaN','M0','Cuc','Eccentricity','Cus','sqrtA','TimeEph',
                'Cic','OMEGA','CIS','Io','Crc','omega','OMEGA DOT','IDOT',
                'CodesL2','GPSWeek','L2Pflag','SVacc','SVhealth','TGD','IODC',
                'TransTime','FitIntvl'])
    nav['sv'] = Series(np.asarray(sv,int), index=nav.index)
    #print (type(nav))

    return nav
    
def writeRinexObs2Hdf(rinex_obs_file_name):
    """
    Function writeObs2Hdf takes the rinex obseravation data .15o and writes 
    the observation data into new hdf .h5 file in the same folder as original
    rinex observation data file. 
    Code is resturctured after Greg Starr's rinexObs function
    """
    filename,ext = splitext(expanduser(rinex_obs_file_name))
    with open(rinex_obs_file_name,'r') as f:
        lines = f.read().splitlines(True)
        lines.append('')
        header,version,headlines,obstimes,sats,svset = scan(lines)
            
    data = processBlocks(lines,header,obstimes,svset,headlines,sats)
    h5fn = filename + '.h5'
    data.to_hdf(h5fn,key='data',mode='w',format='table')
    print('Write succesfull. \n {} is a RINEX {} file, {} kB.'.format(
            rinex_obs_file_name,version,getsize(rinex_obs_file_name)/1000.0))    
    
def readRinexObsHdf(hdf5_file_name):
    """
    Function readObsHdf opens the input .h5 file with raw data structured in
    hdf file. Besides restructured observation data in pandas.panel4D, the 
    function finds the original obseravarion rinex data with .15o extension,
    which has to be in the same folder af hdf file, and reads the header of it.
    Function's output is thus, header data structured as dictionary and 
    observarion data structured as pandas.panel4D type.
    Code is resturctured after Greg Starr's rinexObs function
    """
    
    filename, ext = splitext(expanduser(hdf5_file_name))
    obs_data_ext = '.15o'
    obs_header_file_name = filename + obs_data_ext
    
    with open(obs_header_file_name,'r') as f:
        lines = f.read().splitlines(True)
        lines.append('')
        header,version,headlines,obstimes,sats,svset = scan(lines)
    
    data = read_hdf(hdf5_file_name,key='data')
    
    return header, data, list(svset), obstimes           

def scan(lines):
    """
    Written by greg Starr
    This function sets up the rinex file parsing by quickly running through
    the file, looking for the line at which each time block starts, the time 
    of each block, the satellites in view at each time, and overall what 
    satellites are in the rinex file
    inputs:
        lines - list containing each line in the rinex file as a string
    outputs:
        header - all the header info in a dictionary
        verRinex - the rinex file's version
        headlines - a list of ints, the index of lines where each time block
                    starts
        obstimes - list of times corresponding to each block, same length as 
                    headlines
        sats - the satellites in view at each time, should be same length 
               as headlines
        svset - the set of all the satellites in the rinex file
    """
    
    header={}        
    eoh=0
    for i,line in enumerate(lines):
        if "END OF HEADER" in line:
            eoh=i
            break
        if line[60:].strip() not in header:
            header[line[60:].strip()] = line[:60].strip()
        else:
            header[line[60:].strip()] += " "+line[:60].strip()
    verRinex = float(header['RINEX VERSION / TYPE'].split()[0])
    header['APPROX POSITION XYZ'] = [float(i) for i in header[
        'APPROX POSITION XYZ'].split()]
    header['# / TYPES OF OBSERV'] = header['# / TYPES OF OBSERV'].split()
    header['# / TYPES OF OBSERV'][0] = int(header['# / TYPES OF OBSERV'][0])
    header['INTERVAL'] = float(header['INTERVAL'])

    headlines=[]
    obstimes=[]
    sats=[]
    svset=set()
    i = eoh + 1
    while True:
        if not lines[i]: break
        if not int(lines[i][28]):
            #no flag or flag=0
            headlines.append(i)
            obstimes.append(_obstime([lines[i][1:3],lines[i][4:6],
                                   lines[i][7:9],lines[i][10:12],
                                   lines[i][13:15],lines[i][16:26]]))
            numsvs = int(lines[i][30:32])
            if(numsvs > 12):
                sp=[]
                for s in range(numsvs):
                     	if s == 12 : 
				i += 1
			sp.append(int(lines[i][33+(s%12)*3:35+(s%12)*3]))
                   
                sats.append(sp)
            else:
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(numsvs)])
        
            i+=numsvs*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))+1
        else:
            #there was a comment or some header info
            flag=int(lines[i][28])
            if(flag!=4):
                print(flag)
            skip=int(lines[i][30:32])
            i+=skip+1
    for sv in sats:
        svset = svset.union(set(sv))

    return header,verRinex,headlines,obstimes,sats,svset

def _obstime(fol):
    """
    Written by greg Starr
    turns a listed date collected from the rinex file into a datetime, 
    this is just a utility function.
    """
    year = int(fol[0])
    if (80 <= year <= 99):
        year+=1900
    elif (year < 80): #because we might pass in four-digit year
        year+=2000
    return datetime(year = year, month = int(fol[1]), day = int(fol[2]),
                    hour = int(fol[3]), minute = int(fol[4]),
                    second = int(float(fol[5])),
                    microsecond = int(float(fol[5]) % 1) * 100000
                    )

def _block2df(block,obstypes,svnames,svnum):
    """
    input: block of text corresponding to one time increment INTERVAL of 
    RINEX file output: 2-D array of float64 data from block. Future: consider 
    whether best to use Numpy, Pandas, or Xray.
    """
    nobs = len(obstypes)
    stride = 3

    strio = BytesIO(block.encode())
    barr = np.genfromtxt(strio, delimiter=(14,1,1)*5).reshape((svnum,-1), 
                         order='C')

    data = barr[:,0:nobs*stride:stride]
    lli  = barr[:,1:nobs*stride:stride]
    ssi  = barr[:,2:nobs*stride:stride]

    data = np.vstack(([data],[lli],[ssi])).T

    return data

def processBlocks(lines,header,obstimes,svset,headlines,sats):
    """
    turns the rinex file and the info from scan() into a Panel4D
    inputs:
        the info from scan(), see scan() above
    outputs:
        blocks - the Panel4D with all the data, see above for organization
    """
    obstypes = header['# / TYPES OF OBSERV'][1:]
    blocks = np.nan*np.ones((len(obstypes),max(svset)+1,len(obstimes),3))
    
    for i in range(len(headlines)):
        linesinblock = len(sats[i])*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))
        block = ''.join(lines[headlines[i]+1+int(len(sats[i])/13):headlines[i]+linesinblock+1+int(len(sats[i])/13)])
        bdf = _block2df(block,obstypes,sats[i],len(sats[i]))
        blocks[:,np.asarray(sats[i],int),i,:] = bdf
    
    #print (blocks)
    """
    it is way faster to turn a big numpy array into a Panel4D than 
    to make the Panel4D first and assign it one cell at a time,
    Panel4Ds are slow, it is best to use numpy when possible
    """
    blocks = Panel4D(blocks,
                     labels=obstypes,
                     items=np.arange(max(svset)+1),
                     major_axis=obstimes,
                     minor_axis=['data','lli','ssi'])
    blocks = blocks[:,list(svset),:,:]
    
    return blocks
    
if __name__ == '__main__':
    #h, data = rinexobs('/home/smrak/Documents/TheMahali/rinex/ma132800.15o',
    #            '/home/smrak/Documents/TheMahali/rinex/ma132800.h5')
    readRinexObsHdf('/home/smrak/Documents/TheMahali/rinex/ma132800.h5')
    #print (data['C1', 1, :, 0])