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

Monitoring phenology is important also because plant phenology reacts to a 
warmer climate, where leaf emergence happens earlier, and leaf fall happens 
later due to milder winters. 


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
    \left[\frac{1}{1+\exp(-m_{s}(t-S))} + 
    \frac{1}{1+\exp(m_{A}(t-A))} - 1\right]
    



    
    
.. _Sobrino and Julien (2011): http://www.uv.es/juy/Doc/Sobrino_GIMMS-global-trends_IJRS_2011.pdf

.. _Zhang et al. (2003): http://www.sciencedirect.com/science/article/pii/S0034425702001359
    