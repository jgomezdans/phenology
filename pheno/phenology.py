#!/usr/bin/env python
"""Some functions to plot data, fit models etc."""

import os
import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np
from osgeo import gdal

from pheno_utils import *

def calculate_gdd ( year, tbase=10, tmax=40, latitude=None, longitude=None, \
        fname="/data/geospatial_20/ucfajlg/meteo/temp_2m_unscaled.tif" ):
    """This function calculates the Growing Degree Days for a given year from
    the ERA Interim daily mean surface temperature data. The user can select a 
    base temperature in degrees Celsius. By default, the value is 10."""
    a = 0.0020151192442093
    b = 258.72093867714
    # Check that year range is OK
    assert ( year >= 2001 and year <= 2011 )
    assert ( tmax > tbase )
    if year == 2004 or year == 2008:
        n_doys = 366
    else:
        n_doys = 365
    year = year - 2001
    if latitude is None:
        assert ( longitude is None ) 
        g = gdal.Open ( fname )
        temp = g.ReadAsArray()[((year)*n_doys):((year+1)*n_doys), :, :]
        # Scale to degree C
        temp = np.where ( temp!=-32767, temp*a + b - 273.15, -32767)
        b = np.clip ( temp, tbase, tmax )
        c = np.where ( b-tbase<0, 0, b-tbase )
        agdd = c.cumsum (axis=0)
    else:
        assert ( longitude >= -180 and longitude <= 180 )
        assert ( latitude >= -90 and latitude <= 90 )
        (ix, iy) = pixel_loc ( longitude, latitude )
        g = gdal.Open ( fname )
        temp = g.ReadAsArray()[((year)*n_doys):((year+1)*n_doys), iy, ix ]
        # Scale to degree C
        temp = np.where ( temp!=-32767, temp*a + b- 273.15, -32767)
        b = np.clip ( temp, tbase, tmax )
        c = np.where ( b-tbase < 0, 0, b-tbase )
        agdd = c.cumsum (axis=0)
        
    return ( temp, agdd )

def plot_ndvi ( longitude, latitude ):
    """This function plots NDVI for a given longitude and latitude, for 2001 to
    2011"""
    # Check sanity of longitude an latitude values...
    assert ( longitude >= -180 and longitude <= 180 )
    assert ( latitude >= -90 and latitude <= 90 )
    (ix, iy) = pixel_loc ( longitude, latitude )
    data = []
    t_range = []
    for year in xrange ( 2001, 2012 ):
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