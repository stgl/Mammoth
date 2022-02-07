Mammoth

Repository that documents processing of Eddy Covariance, GNSS Time Series, Strain Calculation, and Figure Production for Hilley, Lewicki, and Baden (2022).

The key files in this repository are:

1) Analysis.ipynb - This notebook was used for initial development of the GPS processing module.  The results from this notebook appear nowhere in the published study.
2) AverageSGSSimulations.ipynb - This notebook was used to compute the source-weight function for the eddy covariance measurements.  The notebook reads a set of sequential gaussian simulations that are based on measured accumulation-chamber flux measurements.  The averaging gives an estimate of the spatial distribution of CO2 fluxes coming from the ground.
3) CalculateHSLStresses.ipynb - This notebook uses the GNSS time-series (that is detrended) to calculate strains across the Long Valley Caldera area.  These strains are then used with Hooke's Law to calculate stresses, given the plane-stress assumption used by GPSGridder.  The notebook should reproduce Figure 4 and Figure S1.
4) CalculateHSLStressesUnfiltered.ipynb - This notebook is identical to CalculateHSLStresses.ipynb, excepting that the GNSS time series is not smoothed using a moving-window Gaussian weighting function.
5) 
