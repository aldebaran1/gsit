# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 20:26:09 2016

@author: Sebsatijan Mrak
smrak@bu.edu
"""
from dateutil import parser
import pytz
import numpy as np
from scipy.interpolate import griddata
from pandas.io.pytables import read_hdf
import glob
from datetime import datetime
from GeoData import GeoData
from GeoData.utilityfuncs import readAllskyFITS
import scipy as sp
import math
import matplotlib.pyplot as plt

def getASIKeogram(ASIfolder, azimuth, altitude, interval, cfg_folder=None):
    """
    Function returns Keogram vectors - intensity, time_vector and elevation vector.
    It takes parameters for azimuth, altitude, time interval and folder with AS
    Images. 
    """
    if cfg_folder == None:
        cfg_folder = '/home/smrak/Documents/TheMahali/asi_cfg/'
    
    timelim = str2posix(interval)
    
    start = interval[0]+interval[1]
    stop = interval[2]+interval[3]
    start = datetime(int(interval[0][-4:]),int(interval[0][:-8]), int(interval[0][-7:-5]), 
                     int(interval[1][:-6]),int(interval[1][-5:-3]), int(interval[1][-2:]))
    stop = datetime(int(interval[2][-4:]), int(interval[2][:-8]), int(interval[2][-7:-5]), 
                     int(interval[3][:-6]),int(interval[3][-5:-3]), int(interval[3][-2:]))
    

    wlstr ='*_0'+ '558' +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(datetime(2015, 10, 7, int(img_t[:-4]),int(img_t[-4:-2]), int(img_t[-2:])))
        
    idx = np.where((image_time >= start) & (image_time <= stop) )[0]
    dt = image_time[idx]
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,
                     (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                      cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                      altitude, timelim))
    
    
    [r,az,el] = g1.dataloc.T
    data2D = g1.data['image']
    az = az%360
    #rel = 90-el
    cut = 1
    #az_conj = (azimuth-180)%360    
    
    if (azimuth + cut) <= 360:
        az1 = np.where((az > azimuth) & (az < azimuth+cut))[0]
        #az1_conj = np.where((az > az_conj) & (az < az_conj+cut))[0]
        #az_ix = np.sort(np.hstack((az1, az1_conj)))
        #find_x = -rel[az_ix]*np.sin(np.radians(az[az_ix]))
        #find_y = rel[az_ix]*np.cos(np.radians(az[az_ix]))
    else:
        print ('Use cut between theta = [0, 360]')
    kg = []
    el_array = []
    for i in range(data2D.shape[1]):
        k = data2D[az1,i]
        el_array.append(el[az1])
        kg.append(k)        
    kg = np.array(kg)
    
    
    return kg, np.array(el_array), dt

def getAllSkyIntensity(ASIfolder, IPPlat, IPPlon, obstimes, altitude, wl,  cfg_folder = None):
    """
    Sebastijan Mrak
    function getAllSkyIntensity returns the intensity of light in Rayleighs [R] at
    given time at the location of IPP, given as input parameters IPPlat and IPPlon.
    Function reads all the .FITS images in given folder and interpolates them in 
    WSG84 coordinate system. All computation is performed as GeoData object. You
    also have to specify the folder of calibration pictures for azimuth and elevation
    as 'cfg_folder'. Besides Intensity, function returns also time array with 
    corresponding timestamps of the image.
    """
    if cfg_folder == None:
        cfg_folder = '/home/smrak/Documents/TheMahali/asi_cfg/'
    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(str(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:]))))
    #print image_time
    dt_match = []                                  
    for t in image_time:
        dt_match.append(int(np.where(obstimes == t)[0]))
    ipp_lat = IPPlat[dt_match]
    ipp_lon = IPPlon[dt_match]
    
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,(cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                                                cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), altitude))
    xcoords = g1.__changecoords__('WGS84')
    latlim=[xcoords[:,0].min(), xcoords[:,0].max()]
    lonlim=[xcoords[:,1].min(), xcoords[:,1].max()]
    N = 256
    latvec = sp.linspace(latlim[0], latlim[1], N)
    lonvec = sp.linspace(lonlim[0], lonlim[1], N)
    [LATM,LONM] = sp.meshgrid(latvec,lonvec)

    #Interpolate ASI to certain altitude
    newcoords = sp.column_stack((LATM.flatten(),LONM.flatten(),
                                 altitude*sp.ones(LONM.size)))    
    g1.interpolate(newcoords,'WGS84',method='linear',twodinterp=True)    
    
    lon_vector = np.linspace(lonlim[0], lonlim[1], N)
    lat_vector = np.linspace(latlim[0], latlim[1], N)
    
    intensity = []
    ASI_data_flattern = g1.data['image']
    #print len(ASI_data_flattern.T)
    for i in range(len(ASI_data_flattern.T)):
        aa = np.reshape(ASI_data_flattern[:,i], (N, N))
        lon_lookup = find_nearest(lon_vector, ipp_lon[i])
        lat_lookup = find_nearest(lat_vector, ipp_lat[i])
        lon_lookup_ix = np.where(lon_vector == lon_lookup)[0]
        lat_lookup_ix = np.where(lat_vector == lat_lookup)[0]
        ASI_value = aa[lon_lookup_ix[0]][lat_lookup_ix[0]] 
        intensity.append(ASI_value)        
        
    
    return image_time, intensity  
    
def writeASIntensity2csv(ASIfolder, CSVfname, IPPlat, IPPlon, obstimes, altitude, cfg_folder = None):
    """    
    Sebastijan Mrak
    function getAllSkyIntensity writes as CSV the intensity of light in Rayleighs [R] at
    given time at the location of IPP, given as input parameters IPPlat and IPPlon.
    Function reads all the .FITS images in given folder and interpolates them in 
    WSG84 coordinate system. All computation is performed as GeoData object. You
    also have to specify the folder of calibration pictures for azimuth and elevation
    as 'cfg_folder'. Besides Intensity, function returns also time array with 
    corresponding timestamps of the image.
    """
    
    fwrite = open(CSVfname, 'w')
    
    if cfg_folder == None:
        cfg_folder = '/home/smrak/Documents/TheMahali/asi_cfg/'
    wlstr ='*_0'+ '558' +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(str(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:]))))
    #print image_time
    dt_match = []                                  
    for t in image_time:
        dt_match.append(int(np.where(obstimes == t)[0]))
    ipp_lat = IPPlat[dt_match]
    ipp_lon = IPPlon[dt_match]
    
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,(cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                                                cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), altitude))
    xcoords = g1.__changecoords__('WGS84')
    latlim=[xcoords[:,0].min(), xcoords[:,0].max()]
    lonlim=[xcoords[:,1].min(), xcoords[:,1].max()]
    N = 256
    latvec = sp.linspace(latlim[0], latlim[1], N)
    lonvec = sp.linspace(lonlim[0], lonlim[1], N)
    [LATM,LONM] = sp.meshgrid(latvec,lonvec)

    #Interpolate ASI to certain altitude
    newcoords = sp.column_stack((LATM.flatten(),LONM.flatten(),
                                 altitude*sp.ones(LONM.size)))    
    g1.interpolate(newcoords,'WGS84',method='linear',twodinterp=True)    

    lon_vector = np.linspace(lonlim[0], lonlim[1], N)
    lat_vector = np.linspace(latlim[0], latlim[1], N)
    
    intensity = []
    ASI_data_flattern = g1.data['image']
    #print len(ASI_data_flattern.T)
    for i in range(len(ASI_data_flattern.T)):
        aa = np.reshape(ASI_data_flattern[:,i], (N, N))
        lon_lookup = find_nearest(lon_vector, ipp_lon[i])
        lat_lookup = find_nearest(lat_vector, ipp_lat[i])
        lon_lookup_ix = np.where(lon_vector == lon_lookup)[0]
        lat_lookup_ix = np.where(lat_vector == lat_lookup)[0]
        ASI_value = aa[lon_lookup_ix[0]][lat_lookup_ix[0]] 
        intensity.append(ASI_value) 
        
        fwrite.write(str(image_time[i]) + ', ' + str(ASI_value) + '\n')
        
    fwrite.close
    print ("Document has been created.")

def find_nearest(array,value):
    """
    Find nearest value in array
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def readASIfromCSV(fname):
    """
    Function reads given CSV file, parse it to intensity and time and converts
    types from string to float and np.dateime64 respectively.
    """
    t, asi = np.loadtxt(fname,dtype='S',unpack=True,delimiter=',')
    t = t.astype(np.datetime64)
    asi = asi.astype(float)
    
    return t, asi
    
def str2posix(timelist):
    """ This will take a list of strings with the date along with a start and
        end time and make a list with the posix times.
        Inputs
            timelist - A list of strings with the data followed by two times.
            The date for the second time can also be used, it will be at index
            2 and the second time will be at index 3.
        Outputs
            dtts - A list of posix times from the original inputs"""
    if len(timelist)==3:
        timelist.insert(2,timelist[0])

    (dt1,dt2) = parser.parse(timelist[0]+ ' '+timelist[1]),parser.parse(timelist[2]+ ' '+timelist[3])
    dt1 =dt1.replace(tzinfo=pytz.utc)
    dt2 = dt2.replace(tzinfo=pytz.utc)
    dt1ts = (dt1 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
    dt2ts = (dt2 -datetime(1970,1,1,0,0,0,tzinfo=pytz.utc)).total_seconds()
    return [dt1ts,dt2ts]
    