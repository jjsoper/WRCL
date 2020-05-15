# Overview
This repository contains the Wachusett Reservoir chloride loading (WRCL) models developed using the USGS [LOADEST]('https://water.usgs.gov/software/loadest/') program. The models are implemented with the R package [rloadest]('https://github.com/USGS-R/rloadest').

#### Data_Init.R
Generates calibration data (discharge and chloride) for input into LOADEST. Discharge may be separated into components of baseflow and runoff.  
#### LOADEST_Init.R
Develops loading models and generates various figures, diagnostic reports, and statistics.  
Key functions:
1. generate.reports
2. model.stats
#### Loading_Analysis.R
Function to set up monthly chloride loading for direct data analysis in R.  
Key functions:
1. load.setup
