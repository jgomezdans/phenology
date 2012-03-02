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
   :scale: 50%
   :alt: A typical NDVI profile
   
   The monthly NDVI trajectory for a site located in the UK for 2002.

   
These measurements do not track individual species, as the footprint of typical sensors
such as MODIS, AVHRR or MERIS is of a few hundred meters, and will usually 
result in a mixture of species. A further complication is that the VI is a
combination of soil (and snow) reflectance and vegetation biochemical composition
and structure, so it is hard to associate any point in a typical trajectory
with particular phenological events such as budburst or leaf fall. In these
cases, simple phenology models are used.


Simple phenology models
=========================

Inspection of typical evolution of vegetation indices indicates that a simple 
phenology model that assumes a quadratic relationship between AGDD and the index
might be appropriate. 

.. math::
    
    NDVI(t) = a\cdot AGDD^{2} + b\cdot AGDD + c 
    
Other more complex models have been developed in the literature, that make use
of different temporal template shapes (such as asymetric Gaussian functions, or
the double logistic function). A double logistic model (after e.g. 
`Zhang et al. (2003)`_ or `Sobrino and Julien (2011)`_ ) is given by 

.. math::
    
    NDVI(t) =   NDVI_{0} + (NDVI_{M} - NDVI_{0} )\cdot
    \left[\frac{1}{1+\exp(-m_{s}(AGDD-S))} + 
    \frac{1}{1+\exp(m_{A}(AGDD-A))} - 1\right]
    



    
    
.. _Sobrino and Julien (2011): http://www.uv.es/juy/Doc/Sobrino_GIMMS-global-trends_IJRS_2011.pdf

.. _Zhang et al. (2003): http://www.sciencedirect.com/science/article/pii/S0034425702001359
    