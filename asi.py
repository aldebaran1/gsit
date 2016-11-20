# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 20:26:09 2016

@author: Sebsatijan Mrak
smrak@bu.edu
"""
import numpy as np
from pandas.io.pytables import read_hdf
import glob
from datetime import datetime
from GeoData import GeoData
from GeoData.utilityfuncs import readAllskyFITS
import scipy as sp
import math
import matplotlib.pyplot as plt

#folder = '/home/smrak/Documents/TheMahali/gnss/mahali/'

def getAllSkyIntensity(ASIfolder, IPPlat, IPPlon, obstimes, altitude, cfg_folder = None):
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
    wlstr ='*_0'+ '558' +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(str(datetime(2015, 10, 07, int(img_t[:-4]), 
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
    
    [lat, lon, alt] = g1.dataloc.T
    lat = np.reshape(lat, (N, N))
    lon = np.reshape(lon, (N, N))
    alt = np.reshape(alt, (N, N))
    
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
        image_time.append(str(datetime(2015, 10, 07, int(img_t[:-4]), 
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
    
    [lat, lon, alt] = g1.dataloc.T
    lat = np.reshape(lat, (N, N))
    lon = np.reshape(lon, (N, N))
    alt = np.reshape(alt, (N, N))
    
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
    print "Document has been created."

def find_nearest(array,value):
    """
    Find nearest value in array
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def openASI():
    fn = 'ASI_rec8_prn23.csv'
    t, data = np.loadtxt(folder+fn,dtype='S',unpack=True,delimiter=',')
    data = data.astype(float)
    #print data
    #mahaliPlot.plot(data)