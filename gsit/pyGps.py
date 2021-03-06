# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 17:28:01 2016

@author: Sebasijan Mrak, Greg Starr
smrak@bu.edu

"""

import numpy as np
import math, time, datetime
from gsit import pyRinex
from scipy import signal, interpolate
from pandas import DataFrame
from pymap3d.coordconv3d import ecef2geodetic,ecef2aer,aer2geodetic
import h5py

#constnats for GPS
f1 = 1575420000
f2 = 1227600000
f5 = 1176450000
c0 = 3E8

def getPRNSlantTEC(P1, P2, units='m'):
    """
    Sebsatijan Mrak
    Function returns slant TEC in TECU units. Input data are PRN information
    at two frequences, where f1 and f2 are difined as global variables. It assumes,
    that you use L1 and L2 frquencies. In case of different GNSS constellation or
    use of L5, correct the indexes. 
    Default config. assumes PRN distance in meters [m], otherwise, fulfill the 
    function parameter 'unit' to correct the units.
    Output units are by default in meters.
    """     
    if units == 'm':
        sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (P1 - P2)) / pow(10,16)
    elif units == 'rad':
        sTEC = ((c0/(40.3*2*np.pi)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (P1/f1 - P2/f2)) / pow(10,16)
            
    elif units == 'cycle':
        sTEC = ((c0/(40.3)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (P1/f1 - P2/f2)) / pow(10,16)        
        
    return sTEC
    
def getPSlantTEC(L1, L2, units = 'cycle'):
    """
    Sebsatijan Mrak
    Function returns slant TEC in TECU units. Input data are phase information
    at two frequences, where f1 and f2 are difined as global variables. It assumes,
    that you use L1 and L2 frquencies. In case of different GNSS constellation or
    use of L5, correct the indexes. 
    Default config. assumes phase information in cycles [cycle], otherwise, 
    fulfill the function parameter 'unit' to correct the units. 
    Output units are by default in meters.
    
    Use only if there is no cycle slips in raw phase file!
    """
    if units == 'cycle':
        sTEC = ((c0/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (L1/f1 - L2/f2)) / pow(10,16)
    elif units == 'rad':
        sTEC = ((c0/(40.3*2*np.pi)) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (L1/f1 - L2/f2)) / pow(10,16)
    elif units == 'm':
        sTEC = ((1/40.3) * (( pow(f2, 2) * pow(f1, 2) ) / 
                (pow(f2, 2) - pow(f1, 2))) * (L1 - L2)) / pow(10,16)
        
    return sTEC

def getPhaseCorrTEC(L1, L2, P1, P2, satbias=None, tec_err=False, 
                    intervals=None, fN = None, maxgap=3, maxjump=2):
    """
    Greg Starr
    Function returns a phase corrected TEC, following a paper by Coco el al:
    'Effect of GPS system biases on differential group delay measurements'.
    Imputs are raw data numpy arrays at two frequencies. Algorithm corrects
    the value of a phase TEC with a difference between mean values of a PRN and
    phase TEC, recpectively. This correction is performed on intervals between
    cycle slips. If the length of consequitive interval is 1, then NaN is inserted
    at this place.    
    """
    #Get intervals between nans and/or cycle slips    
    idx, ranges = getIntervals(L1, L2, P1, P2, maxgap=maxgap, maxjump=maxjump)
#    print ('Intervals between cycle slips ', ranges)
    ERR = np.nan * np.zeros(len(L1))
    TEC = np.nan * np.zeros(len(L1))
    for r in ranges:
        if (r[1] - r[0]) > 1:
            if fN is None:
                f1 = 1575420000
                f2 = 1227600000
                range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2[r[0] : r[1]] - 
                                          P1[r[0] : r[1]]) /40.3 / pow(10, 16)
                phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * \
                         (L1[r[0] : r[1]] / f1 - L2[r[0] : r[1]] / f2) / pow(10, 16)
            else: 
                f1 = (1602 + fN*0.5625) * 1000000
                f2 = (1246 + fN*0.4375) * 1000000
                range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2[r[0] : r[1]] - 
                                          P1[r[0] : r[1]]) /40.3 / pow(10, 16)
                phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * \
                         (L1[r[0] : r[1]] / f1 - L2[r[0] : r[1]] / f2) / pow(10, 16)
            tec_difference = np.array(sorted(phase_tec-range_tec))
            
            tec_difference = tec_difference[np.isfinite(tec_difference)]
            median_difference = tec_difference[int(len(tec_difference)/2)]
            difference_width = tec_difference[int(len(tec_difference)*.75)]-tec_difference[int(len(tec_difference)*.25)]
            median_error = difference_width/np.sqrt(len(tec_difference))
            tec = phase_tec - median_difference
            ERR[r[0]:r[1]] = median_error
            TEC[r[0]:r[1]] = tec
    if (tec_err):
        return TEC, ERR
    elif intervals:
        return TEC, ranges
    else:       
        #TEC[idx]=tec_p
        return TEC 

def getPhaseCorrTECGLONASS(L1,L2,P1,P2):
    range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2 - P1) / 40.3 / pow(10, 16)
    phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * (L1 / f1 - L2 / f2) / pow(10, 16)
    
    tec_difference = np.array(sorted(phase_tec-range_tec))
            
    tec_difference = tec_difference[np.isfinite(tec_difference)]
    median_difference = tec_difference[int(len(tec_difference)/2)]
    difference_width = tec_difference[int(len(tec_difference)*.75)]-tec_difference[int(len(tec_difference)*.25)]
    median_error = difference_width/np.sqrt(len(tec_difference))
    tec = phase_tec - median_difference
    
    return tec
                
        
def getIntervals(L1, L2, P1, P2, maxgap=3,maxjump=2):
    """
    Greg Starr
    scans through the phase tec of a satellite and determines where "good"
    intervals begin and end
    inputs:
        L1, L2, P1, P2 - np arrays of the raw data used to calculate phase corrected
        TEC. 
        maxgap - maximum number of nans before starting new interval
        maxjump - maximum jump in phase TEC before starting new interval
    output:
        intervals - list of 2-tuples, beginning and end of each "good" interval
                    as a Pandas/numpy
    """
    
    r = np.array(range(len(P1)))
    idx = np.isfinite(L1) & np.isfinite(L2) & np.isfinite(P1) & np.isfinite(P2)
    r = r[idx]
    intervals=[]
    if len(r)==0:
        return idx, intervals
    phase_tec=2.85E9*(L1/f1-L2/f2)
    beginning=r[0]
    last=r[0]
    for i in r[1:]:
        if (i-last>maxgap) or (abs(phase_tec[i]-phase_tec[last])>maxjump):
            intervals.append((beginning,last))
            beginning=i
        last=i
        if i==r[-1]:
            intervals.append((beginning,last))
    return idx, intervals

def getVerticalTEC(tec, el, h, Fout=False):
    """
    Sebastijan Mrak
    Function takes the slant TEC numpy array, elevation as numpt array and an
    altitude 'h', to calculate the vertival TEC (vTEC). To map the sTEc to vTEC,
    we follow the thin shell approximation, described by Stefan Shaer in his
    desseration and is also illustrated in Brunini and Azpilicueta:
    'GPS slant total electron content accuracy using the single layer model 
    under different geomagnetic regions and ionospheric conditions'.
    Mapping function is defined as:
    
    F(el) = cos (arcsin (Re / (Re+h) * cos(el)))
    vTEC = sTEC * F(el)
    """
    Re = 6371.0
    h = h
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
    
    if Fout:        
        return np.array(vTEC), np.array(F)
    else:
        return np.array(vTEC)
        
def getROTI(tec, length):
    """
    Sebastijan Mrak
    getROTI returns the rate of TEC Index calculated as the standard deviation 
    of the provided TEC on the moving window of the length 'length'. It returns 
    the ROTI as a numpy array data type.
    """
    roti = []    
    for i in range(len(tec)-length):
        roti.append(np.std(tec[i:i+length]))
    
    return np.array(roti)
    
def getSatellitePosition(rx_xyz, sv, obstimes, navfn, cs = 'wsg84'):
    """
    Sebastijan Mrak
    Function returns satellite position in AER and LLA coordinates for a chosen
    satellite vehicle and obseravtion times. Rx_xyz is a receiver position in 
    ECEF CS, as an np. array. Obstimes gat to be in format pandas.to_datetime
    """   
    navdata = pyRinex.readRinexNav(navfn)
    rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])    
    xyz = getSatXYZ(navdata, sv, obstimes)
    az, el, r = ecef2aer(xyz[:,0],xyz[:,1],xyz[:,2],rec_lat, rec_lon, rec_alt)
    lat, lon, alt = ecef2geodetic(xyz[:,0],xyz[:,1],xyz[:,2])
    if cs == 'wsg84':
        return [lat, lon, alt]
    elif cs == 'aer':
        return [az, el, r]
    else:
        print ('Wrong frame of reference. Type "wsg84" or "aer".')
        return 0;
    
def getIonosphericPiercingPoints(rx_xyz, sv, obstimes, ipp_alt, navfn, cs='wsg84', sattype='G'):
    """
    Sebastijan Mrak
    Function returns a list of Ionospheric Piersing Point (IPP) trajectory in WSG84
    coordinate system (CS). Function takes as parameter a receiver location in 
    ECEF CS, satellite number, times ob observation and desired altitude of a IPP
    trajectory. You also have to specify a full path to th navigation data file.
    It returns IPP location in either WSG84 or AER coordinate system.
    """
    
    ipp_alt = ipp_alt * 1000
    rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
    
    if sattype == 'G':
        navdata = pyRinex.readRinexNav(navfn)
        xyz = getSatXYZ(navdata, sv, obstimes)
        az,el,r = ecef2aer(xyz[:,0],xyz[:,1],xyz[:,2],rec_lat, rec_lon, rec_alt)
        y_out = np.array([az, el, r])
        r_new = []    
        for i in range(len(el)):
            if el[i] > 0:
                r_new.append(ipp_alt / math.sin(math.radians(el[i])))
            else:
                r_new.append(np.nan)
        output = np.array(aer2geodetic(az, el, r_new, rec_lat, rec_lon, rec_alt))
        
    elif sattype == 'R':
        f = h5py.File(navfn, 'r')
        sat_xyz = np.array(f[str(sv)+'/data'])[:,0:3]
        sat_time = np.array(f[str(sv)+'/obstimes'])
        az, el, r = ecef2aer(sat_xyz[:,0]*1E3, sat_xyz[:,1]*1E3, sat_xyz[:,2]*1E3, rec_lat, rec_lon, rec_alt)
        t = datetime2posix(obstimes)
        y_out = []
        
        for y in [az, el, r]:
            x_new = t
            f = interpolate.interp1d(sat_time,y,kind='cubic')
            y_new = f(x_new)
            y_out.append(y_new)
        y_out = np.array(y_out)
        
        r_new = []    
        for i in range(y_out.shape[1]):
            if y_out[1, i] > 0:
                r_new.append(ipp_alt / math.sin(math.radians(y_out[1,i])))
            else:
                r_new.append(np.nan)
        ipp_lat, ipp_lon, ipp_alt = aer2geodetic(y_out[0,:], y_out[1,:], r_new, rec_lat, rec_lon, rec_alt)
        output = np.array(aer2geodetic(y_out[0,:], y_out[1,:], r_new, rec_lat, rec_lon, rec_alt))
    else:
        print ('Type in valid sattype initial. "G" for GPS and "R" for glonass')
        
    if (cs == 'wsg84'):
        return output
    elif (cs == 'aer'):
        return y_out
    else:
        print ('Enter either "wsg84" or "aer" as coordinate system. "wsg84" is default one.')
        
    
def getAllIonosphericPP(data, navdata, obstimes, sv_list, header, ipp_alt):
    """
    WTF?
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

def tec2asiCorrelation(tec_ts, tec, asi_ts, intensity):
    """
    """
    
    
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

def getSatXYZ(nav, sv, times):
    """
    Greg Starr
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

def phaseDetrend(y, order):
    """
    Sebastijan Mrak
    Raw phase data detrending using N-th polinom approximation function.
    Detrended output is input data subtracted with polinom approximation.
    Output is of the same length as input data 'y'. 
    """
    x = range(y.shape[0])
    z = np.polyfit(x, y, order)
    f = np.poly1d(z)
    y_new = f(x)
    y_d = y-y_new
    return y_d
    
def butter_hpf(highcut, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high, btype='highpass', analog=False)
    return b, a

def scintHpf(y, fc=0.1, order=5, fs=1):
    """
    Sebastijan Mrak
    Filter the input data 'y' with desired HP filter.  
    """
    b, a = butter_hpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    return y_filt

def phaseScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Phase scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0]-N)
    for i in range(data.shape[0] - N):
        y[i] = np.std(data[i:i+N])
    return y
    
def AmplitudeScintillation(data, N):
    """
    Sebastijan Mrak
    GNSS Amplitude scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0]-N)
    for i in range(data.shape[0] - N):
        y[i] = np.std(data[i:i+N] / np.mean(data[i:i+N]))
    return y
    
def phaseScintillation(data, fc=0.1, filt_order=6, polyfit_order=3, skip=20):
    """
    Sebastijan Mrak
    Standard pocedure to obtain detrended carrier phase scintillation. Input 
    is a list of carrier phase measurements, than this list goes to a process of
    polynominal detrending and finaly detrended data is high-pass filterd (IIR).
    """
    L = np.nan*np.zeros(data.shape[0])
    idx = np.where(np.isfinite(data))[0]
    L1phi = data[idx]
    L1_d = phaseDetrend(L1phi, polyfit_order)
    Y = scintHpf(L1_d, fc, filt_order)
    
    L[idx[skip:]] = Y[skip:]
    return L

def datetime2posix(dtime):
    """
    Convert an input list of datetime format timestamp to posix timestamp
    """
    return [i.replace(tzinfo=datetime.timezone.utc).timestamp() for i in dtime]

def cycleSlipDetect(y, cslim=50, csmargin=0.05):
    """
    Cycle slip detection amgorithm that uses a 3rd order difference methodology to
    supress random flustuations and magnifies the cycle slip. Output is a list of
    indexes of cycle slips and an estimate of cycle slip value at each index
    """
    l_diff_11 = np.diff(y) # Diff 1st order
    l_diff_11 = np.hstack((np.nan, l_diff_11))
    l_diff_21 = np.diff(l_diff_11) # Diff 2nd order
    l_diff_21 = np.hstack((np.nan, l_diff_21))
    l_diff_31 = np.diff(l_diff_21) # Diff 3rd order 
    l_diff_31 = np.hstack((l_diff_31, np.nan))
    
    ll_1 = np.round(np.nan_to_num(l_diff_31))
    
    cs_ix1= []
    cs_value1 = []
    for i in range(len(ll_1) - 2):
        if abs(ll_1[i+1]) > cslim:
            current_value = np.round(ll_1[i+1])
            previous_value = np.round(ll_1[i])
            next_value = np.round(ll_1[i+2])
            #Even valued samples
            if (current_value % 2) == 0:
                
                if ( 
                    (abs(current_value - np.round(-2 * next_value)) < current_value*csmargin) or 
                    (abs(current_value - np.round(-2 * previous_value)) < current_value*csmargin)
                    ):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)
            #Odd values samples
            else: 
                if (
                    (abs(current_value == np.round(-2 * next_value+1)) < current_value*csmargin) or 
                    (abs(current_value == np.round(-2 * next_value-1)) < current_value*csmargin) or 
                    (abs(current_value == np.round(-2 * previous_value+1)) < current_value*csmargin) or 
                    (abs(current_value == np.round(-2 * previous_value-1)) < current_value*csmargin)
                    ):
                    cs_ix1.append(i+1)
                    cs_value1.append(- (l_diff_31[i+1])/2)

    print (cs_ix1, cs_value1)
    print (y[cs_ix1])
    return cs_ix1, cs_value1