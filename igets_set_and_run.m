%% Set all parameters and run this code
clear;close all;clc

%% Gravimeter and enviroment settings
% Path to hydro-gravimetry library
hydroGravityLib = 'f:\mikolaj\code\libraries\matlab_octave_library';
% Gravimeter/site settings:
station =    'Wettzell';
instrument = 'iGrav006';
latitude =   '49.1449    0.0001 measured';
longitude =  '12.8769    0.0001 measured';
height =     '609.76     0.3000 measured';
author =     'M. Mikolaj (mikolaj@gfz-potsdam.de)';
% Set what steps should be carried out:
%   Load raw TSoft gravimeter data and export the result to monthly ggp/igets
%   files
raw_to_ggp = 1; % 1 = yes, 0 = no
%   Convert 1 second monthly ggp/igets data to 1 minute monthly data 
sec_to_min = 1; % 1 = yes, 0 = no
%   Convert 1 minute monthly ggp/igets data to 1 hour gravity RESIDUALS data 
min_to_res = 1; % 1 = yes, 0 = no
%   Export monthly data in ggp/igets format to one file with whole time
%   series (to tsf, ggp/igets, or mat format)
export_ggp = 1;
% Inport data in tsf or mat format to ggp/igets 
inport_data = 0;
% Set which channels should be exported (e.g., gravity and pressure
% channels in INPUT tsf files)
input_channels = [1,2];% e.g, [15] for PCB temperature;
channel_names = {'gravity','pressure'}; % e.g., {'PCB-Temp'} ;
channel_units = {'V','hPa'};% e.g., {'degC'};
% Process time interval. If filter applied, the time series will be cut to
% remove NaNs at the edges (filter effect)
start_time = [2015 04 01 00 00 00];
end_time   = [2016 08 31 23 59 59];
% Set OUTPUT file format (more or less fixed)
header_offset = 21; % fixed for ggp/igets file format
file_format = 'preterna';
nanval = 99999.999; % Flagged NaN values
out_precision = {'%10.6f','%10.4f'}; % for 1sec and 1min values. 1 hour precision is fixed
    
%% One SECOND data
% Path to raw iGrav data stored in TSoft format. Following data
% structure is required:
%  iGrav gravimeter:
%	'path_igrav'\'prefix_raw'+_YYYY\MMDD\Data_+'prefix_raw'+_MMDD+'prefix_raw'
path_raw = 'y:\iGrav\iGrav006 Data';
% Raw gravimeter data prefix
prefix_raw = 'iGrav006';
% Raw gravimeter data suffix
suffix_raw = '.tsf';% '.tsf' for iGrav and '.030' for SG030
% Resulting 1 second monthly time series folder structure:
%   'path_1sec'\YYYY\'prefix_1sec'+YYYYMM+'suffix_1sec'
path_1sec = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1';
prefix_1sec = 'IGETS-IGRAV-SEC-we006-';
suffix_1sec = '00.aux';
% Besided data, also logfile will be created in the same folder as data:
prefix_1sec_log = 'IGETS-IGRAV-LOG-we006-'; % e.g., 'IGETS-IGRAV-AUXLOG-we006-'
suffix_1sec_log = '00.log';
% Maximum time interval to be automatically interpolated in case of missing 
% or NaN data. Set to 0 for not interpolation. Will be applied before filtering 
% and noted in logfile
fill_missing_1sec = 10; % in seconds 
% Filter from 1 second to 1 minute. This file must be in special format: all 
% header lines start '%' and only half of the impulse response is given (second 
% will be created via mirroring).
filter_file_1sec = fullfile('data','g1s1md.gwr');
% Text to be added to the data file header
header_add_1sec = {'To get (geoid) height of the sensor add 1.05 m (0.03 measured)';...
    'Always check the STATLOG files for details on time series quality'};

%% One MINUTE data
% Following data structure will be used for minute data:
%   'path_1min'\YYYY\'prefix_1min'+YYYYMM+'suffix_1min'
path_1min = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1';
prefix_1min = 'IGETS-IGRAV-MIN-we006-';
suffix_1min = '00.ggp';
% Text to be added to the data file header
header_add_1min = {'To get (geoid) height of the sensor add 1.05 m (0.03 measured)';... 
    'Always check the STATLOG files for details on time series quality'};

%% One HOUR RESIDUAL data
% Following data structure will be used for hourly data:
%   'path_1hour'\YYYY\'prefix_1hour'+YYYYMM+'suffix_1hour'
path_1hour = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3';
prefix_1hour = 'IGETS-IGRAV-RESHOUR-we006-';
suffix_1hour = 'r1.ggp';
% Set output logfile prefix (will be created for each month). Same output folder 
% as for data will be used!
prefix_1hour_log = 'IGETS-IGRAV-RESHOUR-we006-';
suffix_1hour_log = 'r1.log';
% Corrections/Calibration used for hourly time series.
% Calibration parameters:
calib_delay = -11.7;        % seconds (use negative value for SGs)
calib_factor = -914.416;    % nm/s^2/V
% Drift parameters (in polyval matlab/octave format)
drift_fit = [0.2581 -189960];
% Tides file contain time series of tidal effect 
input_tides = 'y:\iGrav\Corrections\Tides\WE_wet2009_TideEffect_CurrentFile_60sec.tsf';
% Set URLs for Atmacs correction (computed using getAtmacs.m):
atmacs_url_loc = {'http://atmacs.bkg.bund.de/data/results/lm/we_lm2_12km_19deg.grav',...
                  'http://atmacs.bkg.bund.de/data/results/iconeu/we_iconeu_70km.grav'};
atmacs_url_glo = {'http://atmacs.bkg.bund.de/data/results/icon4lm/we_icon384_19deg.grav',...
                  'http://atmacs.bkg.bund.de/data/results/icon/we_icon384_20deg.grav'};
% Admittance used for residual pressure effect
atmacs_admittance = -3; % hPa/nm/s^2
% Correction file (with steps and gaps, correctTimeInterval.m will be used)
corr_file = fullfile('data','correction_file_example.txt');
% Maximum time interval to be automatically interpolated in case of missing 
% or NaN data. Set to 0 for not interpolation. Will be applied before filtering 
% and interpolation      
fill_missing_1min = 60; % in mintutes   
% Filter from 1 minute to 1 hour (same as for previous filter)
filter_file_1min = fullfile('data','g1m1h.nlf');
% Text to be added to the data file header   
header_add_1hour = {'Gravity residuals = gravity - tides - atmosphere - polar motion - drift - steps';...
                  'Tides: local model based on SG029 time series (provided by H. Wziontek, BKG)';...
                  'Atmospheric effect: Atmcas + 3*(observed-model pressure)';...
                  'Polar motion: IERS pol coordinates (1.16 amplitude factor)';...
                  'Drift: degree one polynomial fit (empirically estimated)';...
                  'Steps: determined via visual inspection. STATLOG files for details on time series quality';...
                  'All reductions applied before filtering/decimation';...
                  'Gravity is not corrected for PCB-Temp effects (see Level1/AUX time series)';...
                  'Calibration parameters: see IGETS-IGRAV-CAL-we006-20160700.cal';};

%% Export data to other formats 
% Output file (format detected automatically)
file_export = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3\IGETS-IGRAV-RESHOUR-we006-ALL.tsf';
% Inputs = what should be exported
path_export = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3';
prefix_export = 'IGETS-IGRAV-RESHOUR-we006-';
suffix_export = 'r1.ggp';
% Flagged values in output file (used only for ggp/igets file format)
nanval_export = 99999.99;
out_precision_export = '%10.2f';

%% Import data to IGETS (e.g., convert tsf to ggp)
file_inport = 'f:\mikolaj\documents\manuscript\Guentner_et_al_HESS_2017\Data\Input\Meteo\Output_for_IGETS.tsf';%fullfile('data','tides_file_example.tsf');
input_channels_inport = 1:6;
path_inport = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006';
suffix_inport = 'r1.aux';
prefix_inport = 'IGETS-IGRAV-AUX-we006-';
filter_file_inport = []; % e.g., fullfile('data','g1m1h.nlf');
resample_inport = [];%'hour';
fill_missing_inport = 0; % in hours if 'resample_inport' = 'hour'
filter_after_inport = 0;
header_add_inport = {'Precip_gauge: precipitation measured by in-situ tipping bucked and corrected for under-catch following Richter (1995).';...
	'Precip_lysi: precipitation as derived from in-situ lysimeter measurements after AVAT filtering (Peters et al., 2014)';...
	'Evap_lysi: actual evapotranspiration as derived from in-situ lysimeter measurements after AVAT filtering (Peters et al., 2014)';...
	'Evap0: reference evapotranspiration calculated from the in-situ meteorological data with the Penman-Monteith approach following the FAO-56 standard (Allen et al., 1998)';...
	'Discharge: streamflow measured in Chamerau and provided by Bayerisches Landesamt für Umwelt, (2016, http://www.hnd.bayern.de)';... 
	'PCB-corr: PCB correction estimated using PCB-Temperature measurements (add in order to correct gravity)';...
	'For all details on data processing, see Guentner et al., (2017, HESS)'};
out_precision_inport = '%10.2f';
channel_names_inport = {'Precip_gauge','Precip_lysi','Evap_lysi',...
                        'Evap0','Discharge','PCB-corr'};% e.g. {'PCB'};
channel_units_inport = {'mm','mm','mm','mm','mm','nm/s^2'};% e.g. {'degC'};

%%%%%%%%%%%%%%%%%%%%% Call Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(hydroGravityLib)
%% raw_to_ggp
if raw_to_ggp == 1
    igets_raw_to_ggp('start',start_time,'stop',end_time,'input_path',path_raw,...
          'input_prefix',prefix_raw,'input_suffix',suffix_raw,...
          'input_channels',input_channels,'channel_names',channel_names,...
          'channel_units',channel_units,'output_path',path_1sec,'output_prefix',prefix_1sec,...
          'output_suffix',suffix_1sec,'logfile_prefix',prefix_1sec_log,...
          'logfile_suffix',suffix_1sec_log,'file_format',file_format,...
          'station',station,'instrument',instrument,'author',author,...
          'latitude',latitude,'longitude',longitude,'height',height,...
          'header_add',header_add_1sec,'out_precision',out_precision,...
          'fill_missing',fill_missing_1sec,'nanval',nanval,...
          'header_offset',header_offset);
end

%% 1sec_to_1min
if sec_to_min == 1
    igets_sec_to_min('start',start_time,'stop',end_time,'input_path',path_1sec,...
          'input_prefix',prefix_1sec,'input_suffix',suffix_1sec,...
          'input_channels',input_channels,'channel_names',channel_names,...
          'channel_units',channel_units,'output_path',path_1min,...
          'output_prefix',prefix_1min,'output_suffix',suffix_1min,...
          'station',station,'instrument',instrument,'author',author,...
          'latitude',latitude,'longitude',longitude,'height',height,...
          'header_add',header_add_1min,'out_precision',out_precision,...
          'fill_missing',0,'nanval',nanval,'file_format',file_format,...
          'header_offset',header_offset,'filter_file',filter_file_1sec);
end

%% min_to_res
if min_to_res == 1
    igets_min_to_hour_res('start',start_time,'stop',end_time,'input_path',path_1min,...
          'input_prefix',prefix_1min,'input_suffix',suffix_1min,...
          'output_path',path_1hour,'output_prefix',prefix_1hour,...
          'output_suffix',suffix_1hour,'file_format',file_format,...,
          'logfile_prefix',prefix_1hour_log,'logfile_suffix',suffix_1hour_log,...
          'station',station,'instrument',instrument,'author',author,...
          'latitude',latitude,'longitude',longitude,'height',height,...
          'header_add',header_add_1hour,'fill_missing',fill_missing_1min,...
          'nanval',nanval,'header_offset',header_offset,...
          'filter_file',filter_file_1min,'calib_delay',calib_delay,...
          'calib_factor',calib_factor,'drift',drift_fit,...
          'input_tides',input_tides,'atmacs_loc',atmacs_url_loc,...
          'atmacs_glo',atmacs_url_glo,'corr_file',corr_file,...
          'admittance',atmacs_admittance);
end

%% export_ggp
if export_ggp == 1
    igets_export_data('start',start_time,'stop',end_time,'input_path',path_export,...
          'input_prefix',prefix_export,'input_suffix',suffix_export,...
          'output_file',file_export,'out_precision',out_precision_export,...
          'station',station,'instrument',instrument,'author',author,...
          'latitude',latitude,'longitude',longitude,'height',height,...
          'nanval_in',nanval,'nanval_out',nanval_export);
end

%% inport_data
if inport_data == 1
    igets_import_data('start',start_time,'stop',end_time,'input_file',file_inport,...
          'output_path',path_inport,'output_suffix',suffix_inport,...
          'output_prefix',prefix_inport,'out_precision',out_precision_inport,...
          'station',station,'instrument',instrument,'author',author,...
          'latitude',latitude,'longitude',longitude,'height',height,...
          'nanval_in',nanval,'header_add',header_add_inport,...
          'channel_names',channel_names_inport,'channel_units',channel_units_inport,...
          'resample',resample_inport,'filter_after',filter_after_inport,...
          'fill_missing',fill_missing_inport,'input_channels',input_channels_inport,...
          'filter_file',filter_file_inport);
end

rmpath(hydroGravityLib)