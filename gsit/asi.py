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
import glob
from datetime import datetime
from GeoData import GeoData
from GeoData.utilityfuncs import readAllskyFITS
import scipy as sp
import math
import matplotlib.pyplot as plt

def getASIKeogramStatic(ASIfolder, azimuth, altitude, interval, wl, cfg_folder=None, obstimes=None):
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

    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(datetime(2015, 10, 7, int(img_t[:-4]),int(img_t[-4:-2]), int(img_t[-2:])))
    image_time = np.array(image_time)

    idx = np.where((image_time >= start) & (image_time <= stop) )[0]
    dt = image_time[idx]
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,
                     (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                      cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                      altitude, timelim))

    [r, az, el] = g1.dataloc.T
    cut = 1
    az1 = np.where((az > azimuth) & (az < azimuth+cut))[0]
    grid_el = np.arange(15, 88.25, 0.25)
    interpolated_int = []
    for i in range(g1.data['image'].shape[1]):
        img_el = el[az1]
        img_int = g1.data['image'][az1,i]
        new_int = griddata(img_el, img_int, grid_el)
        interpolated_int.append(new_int)
    interpolated_int = np.array(interpolated_int)
    return dt, grid_el, interpolated_int

def getASIKeogramIPP(ASIfolder, azimuth, altitude, interval, wl, 
                     cfg_folder=None, obstimes=None, elevation=None):
    """
    Function returns Keogram vectors - intensity, time_vector and elevation vector.
    It takes parameters for azimuth, altitude, time interval and folder with AS
    Images. Here an azimuth is a list! 
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
    

    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(datetime(2015, 10, 7, int(img_t[:-4]),int(img_t[-4:-2]), int(img_t[-2:])))
    image_time = np.array(image_time)

    idx = np.where((image_time >= start) & (image_time <= stop) )[0]
    dt = image_time[idx]
    azimuth = azimuth[idx]
    if elevation is not None:
        el_out = elevation[idx]

    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,
                     (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                      cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                      altitude, timelim))

    [r, az, el] = g1.dataloc.T
    data = g1.data['image']
    cut = 1
    
    grid_el = np.arange(15, 88.25, 0.25)
    interpolated_int = []
    for i in range(data.shape[1]):
        az1 = np.where((az > azimuth[i]) & (az < azimuth[i]+cut))[0]
        img_el = el[az1]
        img_int = data[az1,i]
        new_int = griddata(img_el, img_int, grid_el)
        interpolated_int.append(new_int)
    interpolated_int = np.array(interpolated_int)
    
    if elevation is not None:
        return dt, grid_el, interpolated_int, el_out
    else:
        return dt, grid_el, interpolated_int


def polarProjection(ASIfolder, altitude, wl, ix):
    """
    Plot Lat-Lon interpolated image, interpolated at ceratin altitude.
    """
    cfg_folder = '/home/smrak/Documents/TheMahali/asi_cfg/'
    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(str(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:]))))

    g1 = GeoData.GeoData(readAllskyFITS,(flist558[ix],
                         (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                          cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                          altitude))
    [r, az, el] = g1.dataloc.T
    az = az%360
    rel = 90-el
    x = -rel * np.sin(np.radians(az))
    y = rel * np.cos(np.radians(az))
    grid_x, grid_y = np.mgrid[-90:90:3600j, -90:90:3600j]
    grid_z = griddata(np.array([x,y]).T, 
                      g1.data['image'][:,0].T, 
                      (grid_x, grid_y), method='linear')
    plt.figure()
    plt.imshow(grid_z.T, origin='lower')
    
def latlonProjection(ASIfolder, altitude, wl, ix):
    """
    Plot Lat-Lon interpolated image, interpolated at ceratin altitude.
    """
    cfg_folder = '/home/smrak/Documents/TheMahali/asi_cfg/'
    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(str(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:]))))

    g1 = GeoData.GeoData(readAllskyFITS,(flist558[ix],
                         (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                          cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                          altitude))
    return g1
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
    
    
    im_int = g1.data['image']
    aa = np.reshape(im_int, (N,N))
    plt.figure()
    plt.imshow(aa.T, origin='lower')


def getAllSkyIntensityAER(ASIfolder, IPPaz, IPPel, altitude, interval, wl, obstimes, cfg_folder=None):
    """
    Sebastijan Mrak
    function getAllSkyIntensity returns the intensity of light in Rayleighs [R] at
    given time at the location of IPP, given as input parameters IPPaz and IPPel.
    Function reads all the .FITS images and gets the closes pixel to given azimuth 
    and elevation. This is returned as intensity list. You have to specify the 
    folder of calibration pictures for azimuth and elevation as 'cfg_folder'. 
    Besides Intensity, function returns also time array with corresponding timestamps 
    of the image.
    """
    if cfg_folder == None:
        cfg_folder = '/home/smrak/Documents/TheMahali/data/asi_cfg/'
    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:])))
    image_time = np.array(image_time)
    # Posix timelim format
    timelim = str2posix(interval)
    start = interval[0]+interval[1]
    stop = interval[2]+interval[3]
    start = datetime(int(interval[0][-4:]),int(interval[0][:-8]), int(interval[0][-7:-5]), 
                     int(interval[1][:-6]),int(interval[1][-5:-3]), int(interval[1][-2:]))
    stop = datetime(int(interval[2][-4:]), int(interval[2][:-8]), int(interval[2][-7:-5]), 
                     int(interval[3][:-6]),int(interval[3][-5:-3]), int(interval[3][-2:]))
    #Match image and IPP timestamps
    dt_match = []                                  
    for t in image_time:
        try:
            dt_match.append(int(np.where(obstimes == t)[0]))
        except:
            continue
    idx = np.where((image_time >= start) & (image_time <= stop) )[0]
    dt = image_time[idx]

    #Iopnospheric piercing points
    IPPaz = IPPaz[dt_match]
    IPPel = IPPel[dt_match]
    #OK
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,
                     (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                      cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                      altitude, timelim))

    [r, az, el] = g1.dataloc.T
    data2D = g1.data['image']
    intensity = np.nan * np.zeros(data2D.shape[1])

    for i in range(data2D.shape[1]):
        data = data2D[:,i]
        intensity[i] = findMin(az, el, IPPaz[i], IPPel[i], data)
    return dt, np.array(intensity)
        
def getAllSkyIntensity(ASIfolder, IPPlat, IPPlon, altitude, interval, wl, obstimes, cfg_folder=None):
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
        cfg_folder = '/home/smrak/Documents/TheMahali/data/asi_cfg/'
    wlstr ='*_0'+ wl +'_*.FITS'
    flist558 = sorted(glob.glob(ASIfolder+wlstr))
    image_time = []
    for image in flist558:
        fn_image = str(image)
        img_t = str(int(np.floor(float(fn_image[-15:-5]))))
        image_time.append(datetime(2015, 10, 7, int(img_t[:-4]), 
                                      int(img_t[-4:-2]), int(img_t[-2:])))
    
    image_time = np.array(image_time)
    # Posix timelim format
    timelim = str2posix(interval)
    start = interval[0]+interval[1]
    stop = interval[2]+interval[3]
    start = datetime(int(interval[0][-4:]),int(interval[0][:-8]), int(interval[0][-7:-5]), 
                     int(interval[1][:-6]),int(interval[1][-5:-3]), int(interval[1][-2:]))
    stop = datetime(int(interval[2][-4:]), int(interval[2][:-8]), int(interval[2][-7:-5]), 
                     int(interval[3][:-6]),int(interval[3][-5:-3]), int(interval[3][-2:]))
    #Match image and IPP timestamps
    dt_match = []                                  
    for t in image_time:
        try:
            dt_match.append(int(np.where(obstimes == t)[0]))
        except:
            continue
    idx = np.where((image_time >= start) & (image_time <= stop) )[0]
    dt = image_time[idx]

    #Iopnospheric piercing points
    ipp_lat = IPPlat[dt_match]
    ipp_lon = IPPlon[dt_match]
    
    # Read all-sky data                                       
    g1 = GeoData.GeoData(readAllskyFITS,(flist558,
                         (cfg_folder+'PKR_DASC_20110112_AZ_10deg.FITS',
                         cfg_folder+'PKR_DASC_20110112_EL_10deg.FITS'), 
                         altitude, timelim))
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
        
    
    return dt, np.array(intensity)
    
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
    
def findNearestIdx(array, value):
    """
    """
    idx = (np.abs(array-value)).argmin()
    return idx
    
def findNearestDate(array, value):
    """
    """
    idx = (np.abs(array-value)).argmin()
    return array[idx]

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
    
def findMin(az, el, IPPaz, IPPel, data):
    """
    Sebastijan Mrak
    Find the best approximation of the line of sight azimuth and elevation on the
    all-sky image pixel. Return the value/intensity of the pixel
    """
    cut = 5
    idx = np.where((az > IPPaz-cut) & (az < IPPaz+cut))[0]
    data = data[idx]
    el = el[idx]
    az = az[idx]
    ######################################################
    idx = np.where((el > IPPel-cut) & (el < IPPel+cut))
    data = data[idx]
    el = el[idx]
    az = az[idx]
    ######################################################
    az_diff = abs(az - IPPaz)
    el_diff = abs(el - IPPel)
    #####################################################
    sum_diff = az_diff + el_diff
    I = np.argmin(sum_diff)
    
    return data[I]
    

def plotPizzacut(az, el, data):
    """
    Plot the cut of the original AER all-sky image.
    """
    rel = 90-el
    az = az%360
    x = -rel * np.sin(np.radians(az))
    y = rel * np.cos(np.radians(az))
    grid_x, grid_y = np.mgrid[-90:90:3600j, -90:90:3600j]
    grid_z = griddata(np.array([x,y]).T, data.T, 
                      (grid_x, grid_y), method='linear')
    plt.figure()
    plt.imshow(grid_z.T, origin='lower')
    
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
    