Scripts for processing IGETS data: **igetsDataTools**
=====================================================
This repository contains **Matlab/Octave** scripts related to loading/writing/processing of gravimeter data stored in **International Geodynamics and
Earth Tide Service** ([IGETS](http://igets.u-strasbg.fr/data_products.php)) database.  
### Dependency
* To run any script, download the [hydroGravityLib](https://github.com/emenems/hydroGravityLib) (Matlab/Octave HydroGravity Library) and adjust the paths in the script accordingly
* All script should work with both, Matlab and Octave, although the loading of one second monthly data in Octave takes lot of time :-(  
* Scripts tested using iGrav006 gravimeter data

### Scripts
* `igets_convert_tsf_to_1sec.m`: convert 1 second iGrav [TSoft](http://seismologie.oma.be/en/downloads/tsoft) files to ggp/igets data format/database, i.e., convert daily data to monthly files  
* `igets_filter_1sec_to_1min.m`: filter and re-sample/decimate the created or downloaded monthly ggp/igets data to 1 minute data
* `igets_create_residual_data_1hour.m`: filter and re-sample the 1 minute data to hourly values
* `igets_load_data.m` load + stack one second/minute/hour data and save the whole time series into one file (ggp/tsf/mat format)  

Use settings part in head of each script to adjust the processing for your gravimeter.

### Auxiliary files
Besides the IGETS data, some auxiliary files are needed. The `data` folder contains examples of such files:
* modified filter files: nonetheless, use (modified) filters recommended by IGETS!
* tides time series: used in `igets_create_residual_data_1hour.m`
* correction file used in `igets_create_residual_data_1hour.m`
