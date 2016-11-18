# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 17:28:01 2016

@author: smrak
"""

import numpy as np
import math
import RinexParser
import mahaliPlot
from pandas import DataFrame
from pymap3d.coordconv3d import ecef2geodetic,ecef2aer,aer2geodetic

def getPRNVerticalTEC(sTEC, navdata, obstimes, sv, header, el_filter = False):
    """
    """
    xyz = getSatXYZ(navdata, sv, obstimes)
    receiver_position = np.asarray(header['APPROX POSITION XYZ'], float)[:, None]
    rlat, rlon, ralt = ecef2geodetic(receiver_position)
    az,el,r = ecef2aer(xyz[:,0],xyz[:,1],xyz[:,2],rlat,rlon,ralt)

    Re = 6371.0
    h1 = 400.0
    rc1 = (Re / (Re+h1))
    #print(rc1)
    F = []
    vTEC1 = []

    for i in range(len(el)):
        if (el_filter == True):    
            if ((el[i] < 10) or (math.isnan(sTEC[i]))):
                vTEC1.append(np.nan)   
                F.append(np.nan)
            else:
                F.append(math.cos(math.asin(rc1*math.sin(math.radians(90-el[i])))))
                vTEC1.append((math.cos(math.asin(rc1*math.cos(math.radians(el[i]))))) * sTEC[i]) 
        else:
            if (math.isnan(sTEC[i])):
                vTEC1.append(np.nan)   
                F.append(np.nan)
            else:
                F.append(math.cos(math.asin(rc1*math.sin(math.radians(90-el[i])))))
                vTEC1.append((math.cos(math.asin(rc1*math.cos(math.radians(el[i]))))) * sTEC[i]) 
    return vTEC1, el, F
    
def getVerticalTEC(tec, el, h):
    Re = 6371.0
    rc1 = (Re / (Re + h))
    vTEC =[]
    F = []
    for i in range(len(tec)):
        if np.isnan(tec[i]):
            vTEC.append(np.nan)
            f = math.cos(math.asin(rc1*math.cos(math.radians(el[i]))))
            F.append(f)
        else:
            f = math.cos(math.asin(rc1*math.cos(math.radians(el[i]))))
            vTEC.append(f * tec[i])
            F.append(f)
        
    return np.array(vTEC), np.array(F)
    
def getSatelliteLatLon(sv):
    navdata = RinexParser.readRinexNav('/home/smrak/Documents/TheMahali/gnss/gps/brdc2800.15n')
    header, data, svset, obstimes =  RinexParser.readRinexObsHdf('/home/smrak/Documents/TheMahali/rinex/mah22800.h5')
    xyz = getSatXYZ(navdata, sv, obstimes)
    lat, lon, alt = ecef2geodetic(xyz)
    
    return lat, lon
    
def getIonosphericPP(data, navdata, obstimes, sv, header, ipp_alt):
    """
    """
    xyz = getSatXYZ(navdata, sv, obstimes)
    rec_position_ecef = np.asarray(header['APPROX POSITION XYZ'], float)[:, None]
    rec_lat, rec_lon, rec_alt = ecef2geodetic(rec_position_ecef)
    az,el,r = ecef2aer(xyz[:,0],xyz[:,1],xyz[:,2],rec_lat, rec_lon, rec_alt)
    r_new = []    
    for i in range(len(el)):
        if el[i] > 10:
            r_new.append(ipp_alt / math.sin(math.radians(el[i])))
        else:
            r_new.append(np.nan)
    ipp_lat, ipp_lon, ipp_alt = aer2geodetic(az, el, r_new, rec_lat, rec_lon, rec_alt)    
    
    return ipp_lat, ipp_lon, ipp_alt
    
def getAllIonosphericPP(data, navdata, obstimes, sv_list, header, ipp_alt):
    """
    """
    r_new = []
    lonPP = []
    latPP = []
    for gps in sv_list:    
        xyz = getSatXYZ(navdata, gps, obstimes)
        rec_position_ecef = np.asarray(header['APPROX POSITION XYZ'], float)[:, None]
        rec_lat, rec_lon, rec_alt = ecef2geodetic(rec_position_ecef)
        az,el,r = ecef2aer(xyz[:,0],xyz[:,1],xyz[:,2],rec_lat, rec_lon, rec_alt)
        for i in range(len(el)):
            r_new.append(ipp_alt / math.sin(math.radians(el[i])))
        max_ix = el.tolist().index(max(el))
        ipp_lat, ipp_lon, ipp_alt = aer2geodetic(az[max_ix], el[max_ix], r_new[max_ix],
                                                 rec_lat, rec_lon, rec_alt)
        r_new = []
        lonPP.append(ipp_lon[0])
        latPP.append(ipp_lat[0])
        
    
    return lonPP, latPP    

def getPRNSlantTEC(data, sv):
    """
    """
    
    f1 = 1575420000
    f2 = 1227600000
    
    P1 = np.array(data['C1', sv, :, 0])
    P2 = np.array(data['P2', sv, :, 0])
    sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
        (pow(f2, 2) - pow(f1, 2))) * (P1 - P2)) / pow(10,16)
        
    return sTEC
    
def getPRNSlantTEC2(rx, sv):
    """
    """
    folder = '/home/smrak/Documents/TheMahali/rinex/'
    f1 = 1575420000
    f2 = 1227600000
    data = RinexParser.readRinexObsHdfData(folder+rx)
    P1 = np.array(data['C1', sv, :, 0])
    P2 = np.array(data['P2', sv, :, 0])
    #L1 = np.array(data['L1', sv, :, 0])
    #L2 = np.array(data['L2', sv, :, 0])
    sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
        (pow(f2, 2) - pow(f1, 2))) * (P1 - P2)) / pow(10,16)

        
    return sTEC
    
def getPSlantTEC(rx, sv):
    """
    """
    folder = '/home/smrak/Documents/TheMahali/rinex/'
    f1 = 1575420000
    f2 = 1227600000
    c = 3.0E8
    data = RinexParser.readRinexObsHdfData(folder+rx)
    L1 = np.array(data['L1', sv, :, 'data'])
    L2 = np.array(data['L2', sv, :, 'data'])
    lol1 = np.array(np.nan_to_num(data['L1', sv, :, 'lli']))
    lol2 = np.array(np.nan_to_num(data['L2', sv, :, 'lli']))
    print sum(lol1), sum(np.where(lol2 == 1))
    #L1 = cycleSlipGPS(L1)
    #L2 = cycleSlipGPS(L2)
    #sTEC = 2.85E9/c * (L1/f1 - L2/f2)
    #L1 = np.array(data['L1', sv, :, 0])
    #L2 = np.array(data['L2', sv, :, 0])
    sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
           (pow(f2, 2) - pow(f1, 2))) * (L1/f1 - L2/f2)) / pow(10,16) * c

        
    return sTEC
    
def getDataTime(rx, gps):
    """
    """
    folder = '/home/smrak/Documents/TheMahali/rinex/'
    data = RinexParser.readRinexObsHdfData(folder+rx)
    t = data.major_axis
    
    return t

def solveIter(mu,e):
    """__solvIter returns an iterative solution for Ek
    Mk = Ek - e sin(Ek)
    adapted to accept vectors instead of single values
    from Bill Rideout's tec.py
    """
    thisStart = np.asarray(mu-1.01*e)
    thisEnd = np.asarray(mu + 1.01*e)
    bestGuess = np.zeros(mu.shape)

    for i in range(5): 
        minErr = 10000*np.ones(mu.shape)
        for j in range(5):
            thisGuess = thisStart + j*(thisEnd-thisStart)/10.0
            thisErr = np.asarray(abs(mu - thisGuess + e*np.sin(thisGuess)))
            mask = thisErr<minErr
            minErr[mask] = thisErr[mask]
            bestGuess[mask] = thisGuess[mask]
        
        # reset for next loop
        thisRange = thisEnd - thisStart
        thisStart = bestGuess - thisRange/10.0
        thisEnd = bestGuess + thisRange/10.0
        
    return(bestGuess)

def getSatXYZ(nav,sv,times):
    """
    getSatelliteXYZ returns the satellite XYZ as a tuple at the inputted times
    inputs are rinex navigation data, satellite number, and list of times
    Output: tuple of satellite position in ECEF coordinates (X,Y,Z)
    Algorithm: Based on http://web.ics.purdue.edu/~ecalais/teaching/geodesy/EAS_591T_2003_lab_4.htm
    also based on Bill Rideout's tec.py
    """
    allSvInfo = nav[nav['sv']==sv] 
    timesarray = np.asarray(times,dtype='datetime64[ms]')
    navtimes = np.asarray(allSvInfo.index,dtype='datetime64[ms]')
    bestephind = np.array([np.argmin(abs(navtimes-t)) for t in timesarray])
    info = np.asarray(allSvInfo)[bestephind]
    info = DataFrame(info,index=times,columns=allSvInfo.columns)
    info['sv'] = sv
    info['gpstime'] = np.array([getGpsTime(t) for t in times])
    # constants
    GM = 3986005.0E8 # universal gravational constant
    OeDOT = 7.2921151467E-5
    
    #Basic Parameters
    t = info['gpstime']-info['TimeEph']
    mu = info['M0']+t*(np.sqrt(GM/info['sqrtA']**6)+info['DeltaN'])
    Ek = solveIter(mu,info['Eccentricity'])  
    Vk = np.asarray(np.arctan2(np.sqrt(1.0-info['Eccentricity'])*np.sin(Ek),
                               np.cos(Ek)-info['Eccentricity']),float)
    PhiK = Vk + info['omega']
    #Correct for orbital perturbations
    omega = np.asarray(info['omega']+info['Cus']*np.sin(2.0*PhiK)
             +info['Cuc']*np.cos(2.0*PhiK),float)
    r = np.asarray((info['sqrtA']**2)*(1.0-info['Eccentricity']*np.cos(Ek))
         +info['Crs']*np.sin(2.0*PhiK)+info['Crc']*np.cos(2.0*PhiK),float)
    i = np.asarray(info['Io']+info['IDOT']*t+info['CIS']*np.sin(2.0*PhiK)
         +info['Cic']*np.cos(2.0*PhiK),float)
    
    #Compute the right ascension
    Omega = np.asarray(info['OMEGA']+(info['OMEGA DOT']-OeDOT)*t-
        (OeDOT*info['TimeEph']),float)
    #Convert satellite position from orbital frame to ECEF frame
    cosOmega = np.cos(Omega)
    sinOmega = np.sin(Omega)
    cosomega = np.cos(omega)
    sinomega = np.sin(omega)
    cosi = np.cos(i)
    sini = np.sin(i)
    cosVk = np.cos(Vk)
    sinVk = np.sin(Vk)
    R11 = cosOmega*cosomega - sinOmega*sinomega*cosi
    R12 = -1.0*cosOmega*sinomega - sinOmega*cosomega*cosi
    #R13 = np.sin(Omega)*np.sin(i)
    R21 = sinOmega*cosomega + cosOmega*sinomega*cosi
    R22 = -1.0*sinOmega*sinomega + cosOmega*cosomega*cosi
    #R23 = -1.0*np.cos(Omega)*np.sin(i)
    R31 = sinomega*sini
    R32 = cosomega*sini
    #R33 = np.cos(i)
          
    xyz = np.zeros((len(times),3))
    rv = np.column_stack((r*cosVk,r*sinVk,np.zeros(r.shape)))
    
    R = np.empty((rv.shape[0],3,3))
    R[:,0,0] = R11
    R[:,0,1] = R12
    R[:,0,2] = 0
    R[:,1,0] = R21
    R[:,1,1] = R22
    R[:,1,2] = 0
    R[:,2,0] = R31
    R[:,2,1] = R32
    R[:,2,2] = 0
    

    for i in range(len(times)): #THIS IS THE SLOWEST PART NOW
        xyz[i,:] = (R[i,:,:].dot(rv[i,:]))
        
    return xyz

def getGpsTime(dt):
    """_getGpsTime returns gps time (seconds since midnight Sat/Sun) for a datetime
    """
    total = 0
    days = (dt.weekday()+ 1) % 7 # this makes Sunday = 0, Monday = 1, etc.
    total += days*3600*24
    total += dt.hour * 3600
    total += dt.minute * 60
    total += dt.second
    return(total)
    
def cycleSlipGPS(data):
    """
    Cycle slip detection after GPS SolutionTheory and Practice book (by Hoffman,
    Lichtenegger and Collins, Ch. 9). Mehtod of high order differences.
    """    
    #gps = 5  #satellite number
   
    l_diff_11 = np.diff(data[150:]) # Diff 1st order
    l_diff_11 = np.hstack((np.nan, l_diff_11))
    l_diff_21 = np.diff(l_diff_11) # Diff 2nd order
    l_diff_21 = np.hstack((np.nan, l_diff_21))
    l_diff_31 = np.diff(l_diff_21) # Diff 3rd order 
    l_diff_31 = np.hstack((l_diff_31, np.nan))

    #mahaliPlot.plot(l_diff_31)
    ll_1 = np.round(np.nan_to_num(l_diff_31))

    ###########################################################################        
    # CS detection L1
    cs_ix1= []
    cs_value1 = []
    for i in range(len(ll_1) - 2):
        if abs(ll_1[i+1]) > 1:
            if (ll_1[i+1] % 2) == 0:
                if ((ll_1[i+1] == -2 * ll_1[i+2]) or (ll_1[i+1] == -2 * ll_1[i])):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)
                elif ( abs(ll_1[i+1]) > 150 and ((ll_1[i+1] != -2 * ll_1[i+2]) or 
                       (ll_1[i+1] == -2 * ll_1[i]))):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)
            else: 
                if ((ll_1[i+1] == -2 * ll_1[i+2]+1) or (ll_1[i+1] == -2 * ll_1[i+2]-1)):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)
                elif ( abs(ll_1[i+1]) > 150 and ((ll_1[i+1] != -2 * ll_1[i+2]) or 
                       (ll_1[i+1] == -2 * ll_1[i]))):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)
    #print cs_ix, cs_value
    if len(cs_ix1) > 0:
        cs_ix1 = np.hstack((cs_ix1, 0))
    lll1 = np.nan*np.zeros(len(data))
    count1 = 0
    CS1 = len(cs_ix1)  
    #print cs_ix1
    i = 0
    # cycle slip repair
    if CS1 > 0:
        while i < cs_ix1[0]:
            lll1[i] = data[i]
            i += 1
        while i < len(data):
            if i == cs_ix1[count1]:
                count1 = count1 + 1
            lll1[i] = data[i] - sum(np.round(cs_value1[0 : count1]))
            i += 1
    else:
        lll1 = data
    
    mahaliPlot.plot(l_diff_31)
    
    return lll1