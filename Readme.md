Mammoth

Repository that documents processing of Eddy Covariance, GNSS Time Series, Strain Calculation, and Figure Production for Hilley, Lewicki, and Baden (2022).

The key files in this repository are:

1) Analysis.ipynb - This notebook was used for initial development of the GPS processing module.  The results from this notebook appear nowhere in the published study.
2) AverageSGSSimulations.ipynb - This notebook was used to compute the source-weight function for the eddy covariance measurements.  The notebook reads a set of sequential gaussian simulations that are based on measured accumulation-chamber flux measurements.  The averaging gives an estimate of the spatial distribution of CO2 fluxes coming from the ground.
3) CalculateHSLStresses.ipynb - This notebook uses the GNSS time-series (that is detrended) to calculate strains across the Long Valley Caldera area.  These strains are then used with Hooke's Law to calculate stresses, given the plane-stress assumption used by GPSGridder.  The notebook should reproduce Figure 4 and Figure S1.
4) CalculateHSLStressesUnfiltered.ipynb - This notebook is identical to CalculateHSLStresses.ipynb, excepting that the GNSS time series is not smoothed using a moving-window Gaussian weighting function.
5) ECFilteredTS.ipynb - Notebook to read, plot, and fit sinusoidal variations (2x daily, daily, yearly) and multi-year (erfc function) to eddy covariance flux measurements.
6) EddyData.py - Python module to load and manipulate eddy covariance data.
7) EddyParserDevelopment.ipynb - Jupyter notebook that was used to develop EddyData.py functionality.
8) FootprintAnalysis.ipynb - Jupyter Notebook to read EC time-series and footprint models (per half-hour metereological data) and analyze whether or not changing footprint source areas could account for changes in EC time series.  Should reproduce Figure 3 in text.
9) GNSSTimeSeriesPlotting.ipynb - Jupyter notebook to read and display example data from GNSS stations in the Long Valley area.
10) GPSDevel.ipynb - Notebook used for development of gpsutils.py (results not used in text).
11) HotCreek.ipynb - Jupyter notebook used to perform preliminary analysis of Hot Creek gauging station data (in Long Valley).  Analysis not used in article.
12) LoadAndProcessGPSData.ipynb - Notebook used to create pickled data files containing Long Valley (and surrounding) GNSS station and network data.
13) LoadAndProcessRegionalGPSData.ipynb - Notebook used to created pickled data files containing greater Sierra Nevada GNSS station and network information.
14) MMFGeometry.ipynb - Jupyter Notebook documenting three-dimensional trace of Mammoth Mountain Fault, as extracted based on mapping cross-referenced with digital topographic data.
15) MakeAnimations.ipynb - Jupyter Notebook to make Movie S1. 
16) MakeAnimationsSierras.ipynb - Jupyter Notebook to create animation of strains throughout the entire Sierra Nevada during the time series analyzed (not presented in contribution).
17) MakeAnimationsSierrasVertical.ipynb - Jupyter notebook to create animation of vertical motions throughout the Sierra Nevada based on GNSS time-series.
18) MakeAnimationsVertical.ipynb - Jupyter notebook to create animation of vertical motions for the Long Valley Area based on GNSS time-series.
19) NoiseAnalysis.ipynb - Exploration of the effect of smoothing window on eddy convariance flux time-series (reproduces Figure contained in Lewicki 2021)
20) VerticalAnalysis.ipynb - Exploration of the correlations between vertical motions and inferred stress changes along Mammoth Mountain Fault (not reported in manuscript).
21) filters.py - module that contains filters that can be applied to time-series.
22) footprints.py - Footprint source-area functions.
23) gpsgridderDevel.ipynb - Notebook used for development of python GPS gridder implementation.
24) gpsutils.py - module for reading and manipulating GPS time series and networks.
25) plotting.py - module containing plotting routines.
26) strains.ipynb - early work on calculating strains.  Approach abandoned in favor of Strain2D modification.
27) utils.py - misc manipulation utilities.
