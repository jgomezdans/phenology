#!/usr/bin/env python
"""Some functions to plot data, fit models etc."""

import os
import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np
from osgeo import gdal

def plot_ndvi ( longitude, latitude ):
    """This function plots NDVI for a given longitude and latitude, for 2001 to
    2011"""
    # Check sanity of longitude an latitude values...
    assert ( longitude >= -180 and longitude <= 180 )
    assert ( latitude >= -90 and latitude <= 90 )
    gdal_dataset = gdal.Open ( "/data/geospatial_20/ucfajlg/MODIS/" + \
            "output/NDVI_2001.tif" )
    geoT = gdal_dataset.GetGeoTransform ()
    igeoT = gdal.InvGeoTransform ( geoT)[1]
    [ix, iy] = gdal.ApplyGeoTransform ( igeoT, longitude, latitude)
    data = gdal_dataset.ReadAsArray ()[ :, iy, ix]
    t_range =  [ matplotlib.dates.datestr2num("2001-%d-01"%m) \
            for m in range(1, 13) ]
    for year in xrange ( 2002, 2012 ):
        gdal_dataset = gdal.Open ( "/data/geospatial_20/ucfajlg/MODIS/" + \
            "output/NDVI_%04d.tif" % year )
        data = np.r_[ data, gdal_dataset.ReadAsArray ()[ :, iy, ix] ] 
        t_range +=  [ matplotlib.dates.datestr2num("%04d-%d-01" % (year, m )) \
            for m in range(1, 13) ]
    plt.plot_date ( t_range, data, '-sr', label="NDVI" )
    plt.grid ( True )
    plt.xlabel("Date")
    plt.ylabel("NDVI [-]")
    plt.show()