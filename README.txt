%% Model description %%
This Snow Cover Runoff (SCR) model uses monthly calculations of snow cover
frequency (SCF) from the Modis instrument along with previous monthly streamflow to predict
monthly streamflow with a one month lead time.
 
This code requires MATLAB. Octave is an open-source option, though this code has only been tested in Matlab.

The SCR model is based upon the paper:
Sproles, E.A., Kerr, T., Orrego Nelson, C. Water Resources Management
(2016) 30: 2581. doi:10.1007/s11269-016-1271-4


SCF is calculated using Google Earth Engine,
and the code can be found on this same github site. 

This sample code is set up to run for the La Laguna sub-watershed of the r?o Elqui
in northern central Chile. The model has also been tested on the John Day
(easten Oregon, USA) and the r?o Aragon in northern Spain.

%% Model structure %%

The DOI for this code is:
10.5281/zenodo.582652

The model has seven parameters that are calibrated using this model based
upon Monte Carlo simulations and Dotty Plots.

The data that accompanies this code has four columns (Year, Month, Q , N)    

The inputs are
Q is a monthly time series of stream flow
N is the snowcover frequency for the basin, and should be a decimal value.
nr is the number of realizations for the monte carlo simulation. A good
number to start with is 1000. 

The model has seven parameters based upon the formula:
Qpredicted = a*((tsSnow(m-1)).^b) + c*tsQshort(m-1) + d*(tsQlong(m-1));

tsSnow is the moving timeseries for SCF. The model optimizes the # of
    months for the time series.
   
tsQshort is the moving timeseries for Q in the short term. The model optimizes the # of
    months for the time series.
   
tsQlong is the moving timeseries for Q in the long term, and conceptually reresents baseflowconditions
    The model optimizes the # of months for the time series.
   
a, c, and d are scaling coefficients specific to each model forcing.
b is an exponential scaling parameter and represents the tapering effect of snowpack contributions to streamflow as it melts (Leibowitz et al. 2012). 
