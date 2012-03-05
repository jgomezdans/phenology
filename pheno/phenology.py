#!/usr/bin/env python
"""Some functions to plot data, fit models etc."""

import pdb
import os
import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np
from osgeo import gdal

from pheno_utils import *

def quadratic_model ( p, agdd ):
    """A quadratic phenology model. Takes in a lenght 3 vector with parameters
    for a quadratic function of AGDD ``agdd``"""
    return p[0]*agdd**2 + p[1]*agdd + p[2]
    
def fourier_model ( p, agdd, n_harm):
    """A simple Fourier model Takes in the Fourier coefficients and phase and 
    uses AGDD ``agdd`` as its time axis. By default ``n_harm`` harmonics are 
    used."""
    
    integration_time = len ( agdd )
    t = np.arange ( 1, integration_time + 1)
    result = t*.0 + p[0]
    for i in xrange ( 1, n_harm*2, 2 ):
        result += p[i]*np.cos ( 2*np.pi*t/integration_time + p[i+1] ) 
    return result
    
def dbl_logistic_model ( p, agdd ):
    """A double logistic model, as in Sobrino and Juliean, or Zhang et al"""
    return p[0] + p[1]* ( 1./(1+np.exp(-p[2]*(agdd-p[3]))) + \
                          1./(1+np.exp(-p[4]*(agdd-p[5])))  - 1 )

def mismatch_function ( p, pheno_model, ndvi, years, n_harm=3 ):
    """The NDVI/Phenology model mismatch function. This can be a multi-year
    function that will be minimised wrt to the VI observations"""
    n_extra_params = 0
    if pheno_model == "quadratic":
        pheno_func = quadratic_model
        n_params = 3
    elif pheno_model == "fourier"
        pheno_func = fourier_model
        n_params = 1 + 2*n_harm
        n_extra_params = 1
    elif pheno_model == "dbl_logistic":
        pheno_func = dbl_logistic_model
        n_params = 6
    # output stores the predictions    
    output = []
    # Initial solution. This will get updated with each year
    xinit = [0,]*n_params
    for year in years:
        ndvi_i = ndvi [ (year-2001)*12:( year - 2001 + 1)*12 ]
        ( temp, agdd ) = calculate_gdd( year, latitude=latitude, \
            longitude=longitude )
        t = (year-2001)*12 + np.array( [ 16,  44,  75, 105, 136, 166, 197, 228,\
            258, 289, 319, 350 ] )
        # We will interpolate NDVI to be daily. For this we need the following array
        ti = (year-2001)*12 + np.arange ( 1, 366 )
        # This is a simple linear interpolator. Strictly, *NOT* needed, but makes
        # everything else easier.
        ndvid = np.interp ( ti, t, ndvi_i )
        if n_extra_params == 0:
            fitness = lambda p, ndvi_in, agdd: ndvi_in - pheno_func ( p, agdd )
            ( xsol, msg ) = leastsq ( fitness, xinit, args=(ndvid, agdd) )
            output.append ( pheno_func ( xsol, agdd ) )
        else:
            fitness = lambda p, ndvi_in, agdd, n_harm: \
                    ndvi_in - pheno_func ( p, agdd, n_harm=n_harm )
            ( xsol, msg ) = leastsq ( fitness, xinit, args=(ndvid, agdd, \
                            n_harm ) )
            output.append ( pheno_func ( xsol, agdd, n_harm=n_harm ) )
        # We can use the solution as the next year's initialisation value
        # This doesn't force this solution, but should be quicker than starting
        # from an arbitrary point
        xinit = xsol.copy()
        
        
def fit_phenology_model ( longitude, latitude, year, pheno_model="quadratic", \
            tbase=10, tmax=40, n_harm=3 ):
    from scipy.optimize import leastsq
    # These next few lines retrieve the mean daily temperature and
    # AGDD, but with
    # This test is to see whether we get a list of years, or a single year
    if isinstance ( year, list ):
        years = year
    elif isinstance ( year, int ) or isinstance ( year, float ):
        years = [ year]
    else:
        raise TypeError, "year has to be a scalar or  list"
    ndvi_all = get_ndvi (  longitude, latitude )
    pdb.set_trace()
    for year in years:
        ndvi = ndvi_all [ (year-2001)*12:( year - 2001 + 1)*12 ]
        ( temp, agdd ) = calculate_gdd( year, latitude=latitude, \
            longitude=longitude )
        # Grab NDVI. Only the first year
        # Need to clear plot
        plt.clf()
        # The following array are the mid-month DoY dates to which NDVI could relate
        # to
        t = (year-2001)*12 + np.array( [ 16,  44,  75, 105, 136, 166, 197, 228,\
                258, 289, 319, 350 ] )
        # We will interpolate NDVI to be daily. For this we need the following array
        ti = (year-2001)*12 + np.arange ( 1, 366 )
        # This is a simple linear interpolator. Strictly, *NOT* needed, but makes
        # everything else easier.
        ndvid = np.interp ( ti, t, ndvi )
        # The fitness function is defined as a lambda function for simplicity
        plt.plot ( ti, ndvid, '-r', label="Obs NDVI" )
        if pheno_model == "quadratic":
            fitf = lambda p, ndvid, agdd: \
                    ndvid - quadratic_model ( p, agdd )
            # Fit fitf using leastsq, with an initial guess of 0, 0, 0
            ( xsol, msg ) = leastsq ( fitf, [0, 0,0], args=(ndvid, agdd) )
            plt.plot ( ti, quadratic_model ( xsol, agdd), \
                '-g', label="Quadratic Fit" )
        elif pheno_model == "fourier":
            fitf = lambda p, ndvi, agdd, n_harm : \
                    ndvi - fourier_model ( p, agdd, n_harm=n_harm )
            ( xsol, msg ) = leastsq ( fitf, [0,]*(2*n_harm+1), \
                    args=(ndvid, agdd, n_harm) )
            plt.plot ( ti, fourier_model ( xsol, agdd, n=_harmn_harm), \
                    '-g', label="Fourier Fit" )
        elif pheno_model == "dbl_logistic":
            fitf = lambda p, ndvi, agdd: \
                    ndvi - dbl_logistic_model ( p, agdd )
            ( xsol, msg ) = leastsq ( fitf, [0,]*6, args=( ndvid, agdd ) )
            plt.plot ( ti, fourier_model ( xsol, agdd ), \
                '-g', label="Logistic Fit" )
        
    plt.rcParams['legend.fontsize'] = 9 # Otherwise too big
    plt.legend(loc='best', fancybox=True, shadow=True ) # Legend
    plt.grid ( True )
    return ( agdd, ndvi, xsol, msg )
    
    
def calculate_gdd ( year=None, , tbase=10, tmax=40, \
        latitude=None, longitude=None, \
        fname="/data/geospatial_20/ucfajlg/meteo/temp_2m_unscaled.tif" ):
    """This function calculates the Growing Degree Days for a given year from
    the ERA Interim daily mean surface temperature data. The user can select a 
    base temperature in degrees Celsius. By default, the value is 10. If no
    year is specified, the whole time series is retrieved. Note that if you 
    don't specify time and location, the operation can be quite slow and will
    return lots of data."""
    a = 0.0020151192442093
    b = 258.72093867714
    # Check that year range is OK
    assert ( year >= 2001 and year <= 2011 or year is None )
    assert ( tmax > tbase )
    if year == 2004 or year == 2008:
        n_doys = 366
    else:
        n_doys = 365
    year = year - 2001
    if latitude is None:
        assert ( longitude is None ) 
        g = gdal.Open ( fname )
        if year is not None:
            temp = g.ReadAsArray()[((year)*n_doys):((year+1)*n_doys), :, :]
        else:
            temp = g.ReadAsArray()
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
        if year is not None:
            temp = g.ReadAsArray()[((year)*n_doys):((year+1)*n_doys), iy, ix ]
        else:
            temp = g.ReadAsArray()[:, iy, ix ]
        # Scale to degree C
        temp = np.where ( temp!=-32767, temp*a + b- 273.15, -32767)
        b = np.clip ( temp, tbase, tmax )
        c = np.where ( b-tbase < 0, 0, b-tbase )
        agdd = c.cumsum (axis=0)
        
    return ( temp, agdd )

def get_ndvi ( longitude, latitude, plot=False ):
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
    if plot:
        plt.plot_date ( t_range, data, '-sr', label="NDVI" )
        plt.grid ( True )
        plt.xlabel("Date")
        plt.ylabel("NDVI [-]")
        plt.show()
    return ( data )