================================================================
Fitting a simple phenology model to satellite observations
================================================================

Phenology refers to the study of timing of recurrent events, such as
leaf emergence, leaf fall, etc. In the context of carbon science, the timing of
emergence of leaves and the duration of the leaf-on period is crucial to 
understand the impact of decidious trees on atmospheric carbon levels. The 
timing of vegetation is also important from the point of view of other ecosystem 
components: vegetation flowering dates and fruit production are an important 
constraint on the temporal behaviour of fruit-eating wildlife, for example. 
Phenology is also important in agricultural practices, as certain management 
practices are only applied at particular phenological stages of crop development.

Monitoring phenology is important also because plant phenology reacts to a 
warmer climate, where leaf emergence happens earlier, and leaf fall happens 
later due to milder winters. Historically, researchers have recorded bud burst,
date of flowering and leaf fall, breeding times, etc. for individual species. 
Phenology networks such as the `USA's phenology network <http://www.usanpn.org/home>`_
or efforts such as the `UK's NatureCalendar <http://www.naturescalendar.org.uk/>'_
rely on meticulous observation of phenology events by dedicated researchers and
members of the public. While these efforts are invaluable, and often allow the
extension of data series way back in the past, they tend to only monitor a handful
of species over a limited geographical extent. Vegetation indices, such as NDVI 
or EVI, broadly respond to the amount of vegetation within the scene, and 
typical temporal trajectories can be used to track the temporal development of 
vegetation and can be related to phenology measurements in the ground. 

.. figure:: ndvi_annual_trend.png
   :scale: 25%
   :alt: A typical NDVI profile
   
   The monthly NDVI trajectory for a site located in the UK for 2002.

   
These measurements do not track individual species, as the footprint of typical sensors
such as MODIS, AVHRR or MERIS is of a few hundred meters, and will usually 
result in a mixture of species. A further complication is that the VI is a
combination of soil (and snow) reflectance and vegetation biochemical composition
and structure, so it is hard to associate any point in a typical trajectory
with particular phenological events such as budburst or leaf fall. In these
cases, simple phenology models can be used to interpret the signal.

Plotting NDVI temporal trajectories
-------------------------------------

Some data and code has been provided for you to visualise VI trajectories. We
shall be using the MODIS monthly NDVI product, gridded to a global grid of 
:math:`1.5\times 1.5` degrees. The function ``plot_ndvi`` will plot the complete
time series for a given logitude and latitude. For example, this plots the
NDVI for the gridcell that covers the `Hainich <http://www.bgc-jena.mpg.de/public/carboeur/sites/hainich.html>`_
FluxNet site, an area in Easter Germany where trees are mostly decidious (beech).

.. plot::
   :include-source: 
    
   from phenology import *
   get_ndvi ( 10.25, 51.05, plot=True ) # NDVI around Hainich FluxNET site
   plt.show()

The following plot shows some location near the city of Tomsk, in Siberia.

.. plot::
   :include-source:  
   
   from phenology import *
   get_ndvi ( 86, 57, plot=True) # NDVI around Siberia
   plt.show()

This is the NDVI temporal trajectory for the Tumbarumba Fluxnet site in Australia

.. plot::
   :include-source:  
   
   from phenology import *
   get_ndvi ( 148, -35, plot=True) # Tumbarumba, Oz
   plt.show()     

Try a few sites with different latitude and longitudes, and discuss what the
differences are.


Simple phenology models
=========================

Phenology has been found to be very responsive to temperature. *Growing degree
days* are a useful way to define phenological events, rather than calendar days, 
as the latter are rather arbitrary. Growing degree days accumulate anytime the 
average temperature for the day is more than a certain base temperature. This 
allows for some normalisation of timing of events across latitudinal ranges, for
example. 

.. todo::
    
   Put a few references on climate and phenology and how it impacts C cycle.

Accessing the AGDD data
------------------------

Data from the `ERA interim archive <http://data-portal.ecmwf.int/data/d/interim_daily/>`_
have been prepared. The data available for this practial is the mean daily 
temperature. Some functions have been provided for you to access the data easily:

.. plot::
   :include-source: 
  
   # Import some libraries, in case you haven't yet imported them
   import matplotlib.pyplot as plt
   import numpy as np
   from phenology import *
   # These next few lines retrieve the mean daily temperature and
   # AGDD for the three sites mentioned above
   ( temp_hainich, agdd_hainich ) = calculate_gdd( 2005, \
            latitude=51, longitude=10 )
   ( temp_tomsk, agdd_tomsk ) = calculate_gdd( 2005, \
            latitude=57, longitude=86 )
   ( temp_tumbarumba, agdd_tumbarumba ) = calculate_gdd( 2005, \
            latitude=-35, longitude=148 )
   # Temporal range for plots
   t_range =  np.arange ( 1, 366 )
   # First subplot is Hainich (DE)
   plt.subplot ( 3, 1, 1)
   # Put a grey area for the AGDD calculation bounds
   plt.axhspan ( 10, 40, xmin=0, xmax=366, color='0.9' )
   # Plot temperature
   plt.plot ( t_range, temp_hainich, '-r', label="Tm" )
   plt.ylabel("Mean Temp [degC]")
   plt.grid ( True )
   plt.twinx()
   plt.plot ( t_range, agdd_hainich, '-g', label="AGDD" )
   plt.ylabel ( "AGDD [degC]")
   plt.subplot ( 3, 1, 2)
   # Second subplot is Tomsk. Everything as before
   plt.axhspan ( 10, 40, xmin=0, xmax=366, color='0.9' )
   plt.plot ( t_range, temp_tomsk, '-r', label="Tm" )
   plt.ylabel("Mean Temp [degC]")
   plt.grid ( True )
   plt.twinx()
   plt.plot ( t_range, agdd_tomsk, '-g', label="AGDD" )
   plt.ylabel ( "AGDD [degC]")
   plt.subplot ( 3, 1, 3)
   # Third subplot is Tumbarumba. Everything as before
   plt.axhspan ( 10, 40, xmin=0, xmax=366, color='0.9' )
   plt.plot ( t_range, temp_tumbarumba, '-r', label="Tm" )
   plt.ylabel("Mean Temp [degC]")
   plt.grid ( True )
   plt.rcParams['legend.fontsize'] = 9 # Otherwise too big
   plt.legend(loc='best', fancybox=True, shadow=True ) # Legend
   plt.twinx()
   plt.plot ( t_range, agdd_tumbarumba, '-g', label="AGDD" )
   plt.ylabel ( "AGDD [degC]")
   plt.xlabel("DoY/2005")
   plt.rcParams['legend.fontsize'] = 9 # Otherwise too big
   plt.legend(loc='best', fancybox=True, shadow=True ) # Legend
   plt.show()
            
Examine the previous plots, noting particularly the inflexion points in the 
AGDD curve, and how they relate to the base and maximum mean daily temperatures
(shown in the grey area). Also not how for the Tumbarumba site, there is a 
seasonality with respect to the Northern Hemisphere sites. 

Phenology models
-----------------

Inspection of typical evolution of vegetation indices like the one carried out
above suggest that a simple phenology model that goes from a minimum to a maximum
and then decreases again may be suitable, at least for the Northern Hemisphere.
One such method used successfully in `De Beurs and Henebry (2008)`_ by using
a simple quadratic function of AGDD:

.. math::
    
    NDVI(t) = a\cdot AGDD^{2} + b\cdot AGDD + c 

.. plot::
   :include-source:
       
   # Import some libraries, in case you haven't yet imported them
   import matplotlib.pyplot as plt
   import numpy as np
   from scipy.optimize import leastsq
   from phenology import *
   # These next few lines retrieve the mean daily temperature and
   # AGDD, but with 
   ( temp_tomsk, agdd_tomsk ) = calculate_gdd( 2005, latitude=57, longitude=86 )
   # Grab NDVI. Only the first year
   ndvi = plot_ndvi (  86, 57 )[:12]
   # Need to clear plot
   plt.clf()
   # The following array are the mid-month DoY dates to which NDVI could relate
   # to
   t = np.array([ 16,  44,  75, 105, 136, 166, 197, 228, 258, 289, 319, 350])
   # We will interpolate NDVI to be daily. For this we need the following array
   ti = np.arange ( 1, 366 )
   # This is a simple linear interpolator
   ndvid = np.interp ( ti, t, ndvi )
   # The fitness function is defined as a lambda function for simplicity
   fitf = lambda p, ndvid, agdd: \
        ndvid -( p[0]*agdd**2 +p[1]*agdd+p[2])
   # Fit fitf using leastsq, with an initial guess of 0, 0, 0
   ( xsol, msg ) = leastsq ( fitf, [0, 0,0], args=(ndvid, agdd_tomsk) )
   plt.plot ( ti, ndvid, '-r', label="Obs NDVI" )
   plt.plot ( ti, xsol[0]*agdd_tomsk*agdd_tomsk +xsol[1]*agdd_tomsk+xsol[2], \
       '-g', label="Fit" )
   # Now, try a different base temperature
   ( temp_tomsk, agdd_tomsk ) = calculate_gdd( 2005, tbase=-5,latitude=57, \
        longitude=86 )
   ( xsol, msg ) = leastsq ( fitf, [0, 0,0], args=(ndvid, agdd_tomsk) )
   plt.plot ( ti, xsol[0]*agdd_tomsk*agdd_tomsk +xsol[1]*agdd_tomsk+xsol[2], \
      '-b', label="Fit (-5)" )
   plt.rcParams['legend.fontsize'] = 9 # Otherwise too big
   plt.legend(loc='best', fancybox=True, shadow=True ) # Legend
   plt.grid ( True )
   plt.show()     
    
We can see that the quadratic model has some complications even fitting a simple
NDVI profile like that of Siberia. The model, as introduced above, will also 
struggle to cope with NDVI patterns typical of the Southern Hemisphere (unless
a temporal shift is introduced). These limitations have lead to the development
of more complex and robust methods. 
    
Other more complex models have been developed in the literature, that make use
of different temporal template shapes (such as asymetric Gaussian functions, or
the double logistic function). A double logistic model (after e.g. 
`Zhang et al. (2003)`_ or `Sobrino and Julien (2011)`_ ) is given by 

.. math::
    
    NDVI(t) =   NDVI_{0} + (NDVI_{M} - NDVI_{0} )\cdot
    \left[\frac{1}{1+\exp(-m_{s}(AGDD-S))} + 
    \frac{1}{1+\exp(m_{A}(AGDD-A))} - 1\right]
    
In comparison with the quadratic model, the double logistic model will provide a
more flexible fit to the data by virtue of having 6 parameters (the :math:`NDVI_{0}`
and :math:`NDVI_{M}` terms are maximum and minimum NDVI, and can be readily 
estimated from the time series).

Other methods to fit a model to observations of NDVI rely on Fourier analysis
ideas. Fourier analysis states that within a closed interval, any periodic
function can be expressed as a sum of increasing frequency sine waves:
    
.. math::
    
    NDVI(t) = \overline{NDVI}(t) + \sum_{i=1}^{N/2}A_{i}\cos(2\pi i t/N) + \phi_{i}
    
where :math:`\overline{NDVI}(t)` is the mean NDVI value within the period of 
interest :math:`(0,N)`. :math:`A_{i}` and :math:`\phi_{i}` are the magnitude and
phase of the :math:`i`-th harmonic, respectively. Usually, only a few terms of the
summation are required to produce a reasonable fit to the observations. An added
benefit is that the different harmonics allow for a more detailed exploration of
the temporal dynamics observed by the sensor: the first term (the 0-th harmonic)
can be related to the mean biome amount of vegetation. The first and second
harmonics relate to the dynamics of annual and biannual evoluation of vegetation.
Finally, frequency-domain analysis is fairly robust against noise. However, there
are some shortcomings: data gaps need to be filled in or "padded", and the
frequency at which one can extract information is governed by the periodicity of
the data, wihch in our case is monthly. Also, fast events might be blurred. For
a more in-depht analysis, see e.g. `Moody and Johnson (2001)`_


    
    
.. _De Beurs and Henebry (2008): http://geography.vt.edu/deBeurs_Henebry_JClimate.pdf

.. _Sobrino and Julien (2011): http://www.uv.es/juy/Doc/Sobrino_GIMMS-global-trends_IJRS_2011.pdf

.. _Zhang et al. (2003): http://www.sciencedirect.com/science/article/pii/S0034425702001359
.. _Moody and Johnson (2001): ftp://ftp.ccrs.nrcan.gc.ca/ftp/ad/Phenology/PhenologyPapers/Moody_2001_AVHRR_DFourierTransPhenology_USA.pdf