# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 20:26:09 2016

@author: smrak
"""
import numpy as np
from pandas.io.pytables import read_hdf
import glob
import mahaliPlot
from datetime import datetime
from GeoData import GeoData
from GeoData.utilityfuncs import readAllskyFITS
import scipy as sp
import math
import matplotlib.pyplot as plt

folder = '/home/smrak/Documents/TheMahali/gnss/mahali/'
    
def parseASI(fname, gps, altitude, rx):
    data = read_hdf(fname, key='data')
    gps_t = data.major_axis
    IPPlon = data['lon', gps, :, 'data']  
    IPPlat = data['lat', gps, :, 'data']
    ASfn = '/home/smrak/Documents/TheMahali/Allsky_t2/'
    wlstr ='*_0'+ '558' +'_*.FITS'
    flist558 = sorted(glob.glob(ASfn+wlstr))
    AS_list = []
    AS_time = []
    tt = []    
    for i in range(len(IPPlon)):
        tt.append(str(datetime.strptime(str(gps_t[i]),'%Y-%m-%d  %H:%M:%S')))
    gps_time_array_string = np.array(tt)
    wr_fn = 'ASI_rec' + str(rx) +'_prn' + str(gps) + '.csv'
    fwrite = open(wr_fn,'w')    
    
    for image in flist558:
        g1 = GeoData.GeoData(readAllskyFITS,(image,('PKR_DASC_20110112_AZ_10deg.FITS',
                                                    'PKR_DASC_20110112_EL_10deg.FITS'), altitude))
                                                    
        #ASI image transform to WGS84
        xcoords = g1.__changecoords__('WGS84')
        latlim=[xcoords[:,0].min(), xcoords[:,0].max()]
        lonlim=[xcoords[:,1].min(), xcoords[:,1].max()]
        nlat = 256
        nlon = 256
        latvec = sp.linspace(latlim[0],latlim[1],nlat)
        lonvec = sp.linspace(lonlim[0],lonlim[1],nlon)
        [LATM,LONM] = sp.meshgrid(latvec,lonvec)
        #Interpolate ASI to certain altitude
        newcoords = sp.column_stack((LATM.flatten(),LONM.flatten(),150.*sp.ones(LONM.size)))    
        g1.interpolate(newcoords,'WGS84',method='linear',twodinterp=True)    
        #Get a lat lon vectors of the image
        [lat, lon, alt] = g1.dataloc.T
        lat_min = min(lat)
        lat_max = max(lat)
        d_lat = np.diff(lat).max()
        lon_min = min(lon)
        lon_max = max(lon)
        d_lon = np.diff(lon).max()
        lat = np.reshape(lat, (256, 256))
        lon = np.reshape(lon, (256, 256))
        alt = np.reshape(alt, (256, 256))
    
        lon_vector = np.linspace(lon_min, lon_max, 256)
        lat_vector = np.linspace(lat_min, lat_max, 256)
        #Get time of the image!
        fn_str = str(image)
        #print fn_str
        img_t = str(int(np.floor(float(fn_str[-15:-5]))))
        #print img_t
        img_t_dt = str(datetime(2015, 10, 07, int(img_t[:-4]), int(img_t[-4:-2]), int(img_t[-2:])))
        # Get time stamp from GPS where it is equal to timestamp of ASI image
        gps_time_ix = np.where(gps_time_array_string == img_t_dt)[0]
        
        ipp_lon_temp = float(IPPlon[gps_time_ix])
        ipp_lat_temp = float(IPPlat[gps_time_ix])
        lon_lookup = find_nearest(lon_vector, ipp_lon_temp)
        lon_lookup_ix = np.where(lon_vector == lon_lookup)[0]
        lat_lookup = find_nearest(lat_vector, ipp_lat_temp)
        lat_lookup_ix = np.where(lat_vector == lat_lookup)[0]
        
        #print ipp_lat_temp, ipp_lon_temp
        a = g1.data['image']
    
        aa = np.reshape(a, (256, 256))
        #print aa[129][120]
        #plt.imshow(aa) 
        ASI_value = aa[lon_lookup_ix[0]][lat_lookup_ix[0]] 
        AS_list.append(ASI_value)
        AS_time.append(img_t_dt)
        
        print ipp_lat_temp, ipp_lon_temp, ASI_value
        #print ASI_value
        fwrite.write(str(img_t_dt) + ', ' + str(ASI_value) + '\n')
        
    fwrite.close

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
    print data
    mahaliPlot.plot(data)
    
if __name__ == '__main__':
    parseASI('mah72800_alldataIPP130.h5', 9, 130, 7)
    parseASI('mah62800_alldataIPP130.h5', 9, 130, 6)
    parseASI('mah132800_alldataIPP130.h5', 9, 130, 13)
    #parseASI('mah42800_alldataIPP130.h5', 23, 130, 4)
    #parseASI('mah32800_alldataIPP130.h5', 23, 130, 3)
    #parseASI('mah22800_alldataIPP130.h5', 23, 130, 2)
    #parseASI('mah132800_alldataIPP130.h5', 23, 130, 13)