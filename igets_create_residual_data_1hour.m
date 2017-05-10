%% Convert 1 minute data to 1 hour residuals
% This script will convert one minute monthly ggp/igets files to one hour 
% ggp/igets gravity residuals.
% Use the 'igets_convert_tsf_to_1sec.m' first to get the monthly 1 second 
% data. Then use 'igets_filter_1sec_to_1min.m' for conversion to 1 minute
% data.
% Following data structure of INPUT files is required
%     'input_path\YYYY\input_prefix+YYYYMM+input_suffix'
% Following OUTPUT structure will be used
%     'output_path\YYYY\output_prefix+YYYYMM+output_suffix'
% Following steps are carried out:
%   1. Load data for all months and stack to one time series
%   2. Load/compute correction time series:
%       - tides (given file)
%       - atmcas (requires Internet connection)
%       - polar motion (requires Internet connection)
%       - drift (given polynomial coef.)
%   3 If required, the input time series will be automatically 
%      inspected for missing/NaN data and specified interval can be filled 
%      by interpolated values. This should be normally done while creating
%      1 second data.
%   4. Filter the loaded time series (if filter selected). 
%   5. Re-sample/decimate the time series to 1 hour resolution
%   6. Write the result
% Optional steps (before writing results): 
%   * Interpolate user specified (not automatic detection) time intervals,
%      or insertion of NaNs for unreliable time intervals (see 'corr_file')
%   * Remove user specified steps (see 'corr_file')
% 
% Script tested on Matlab R2015b (preferred) and Octave 4.2.1 (rather slow)
%
%                                                    M.Mikolaj
%                                                    mikolaj@gfz-potsdam.de
clear
close all
clc
% Add path containing hydroGravity library (loadtsf.m, findTimeStep.m, ... functions)
% Download from: https://github.com/emenems/hydroGravityLib
addpath('f:\mikolaj\code\libraries\matlab_octave_library')


%% Main INPUT settings
% Process time interval (to be loaded)
start_time = [2015 03 05 14 00 00];
end_time   = [2017 02 27 15 00 00];
% INPUT File path/name settings (to tsf files)
input_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1'; % year/month/day will be generated automatically
input_prefix = 'IGETS-IGRAV-MIN-we006-'; % file name prefix
input_suffix = '00.ggp';
% Set which channels should be loaded (e.g., gravity and pressure))
input_channels = [1,2]; % 
% Set filter from 1 min to 1 hour. This file must be modified: all header 
% lines start  '%' and only half of the impulse response is given (second 
% will be created via mirroring).
filter_file = fullfile('data','g1m1h.nlf');

% Set tides/atmo/polar motion/drift corrections.
% See \data\tides_example_file.tsf how the 'input_tides' file should look
% like.
input_tides = 'y:\iGrav\Corrections\Tides\WE_wet2009_TideEffect_CurrentFile_60sec.tsf';
% Set URLs for Atmacs correction (computed using getAtmacs.m):
atmacs_url_link_loc = {'http://atmacs.bkg.bund.de/data/results/lm/we_lm2_12km_19deg.grav',...
                       'http://atmacs.bkg.bund.de/data/results/iconeu/we_iconeu_70km.grav'};
atmacs_url_link_glo = {'http://atmacs.bkg.bund.de/data/results/icon4lm/we_icon384_19deg.grav',...
                       'http://atmacs.bkg.bund.de/data/results/icon/we_icon384_20deg.grav'};
atmcas_admittance = -3; % nm/s^2/hPa
% Polar motion effect (computed using getEOPeffect.m):
pol_lat = 49.14490;
pol_lon = 12.87687;
% Drift parameters (in polyval format)
drift_fit = [0.2581 -189960];
% Calibration parameters:
calib_delay = -11.7; % seconds (use negative value for SGs)
calib_factor = -914.416; % nm/s^2/V
% Correction file (with steps and gaps, correctTimeInterval.m will be used)
corr_file = 'y:\iGrav\Corrections\iGrav006_correction_file_plotGrav.txt';
% Maximum time interval to be automatically interpolated in case of missing 
% or NaN data. Set to 0 for not interpolation. 
fill_missing = 60; % minutes (=in input units)

%% Main OUTPUT settings
% Output time series (to be written)
start_time_write = [2015 04 02 12 00 00];
end_time_write   = [2017 02 27 12 00 00];
% Set OUTPUT file naming
output_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3';
output_prefix = 'IGETS-IGRAV-RESHOUR-we006-';
output_suffix = 'r1.ggp';
% Set output logfile prefix (will be created for each month). To avoid
% writing logfile, set to '' (=empty); Same output folder as for data will
% be used!
logfile_prefix = 'IGETS-IGRAV-RESHOUR-we006-';
logfile_suffix = 'r1.log';

% Set header for OUTPUT file
header = {'Filename','';... % file name will be appended automatically
            'Station','Wettzell';...
            'Instrument','iGrav006';...
            'N. Latitude (deg)', '49.1449    0.0001 measured';...
            'E. Longitude (deg)','12.8769    0.0001 measured';...
            'Elevation MSL (m)', '609.76     0.3000 measured';...
            'Author','M. Mikolaj (mikolaj@gfz-potsdam.de)'};      
% Header footer (OUTPUT)
header_add = {'Gravity residuals = gravity - tides - atmosphere - polar motion - drift - steps';...
              'Tides: local model based on SG029 time series';...
              'Atmospheric effect: Atmcas + 3*(observed-model pressure)';...
              'Polar motion: IERS pol coordinates (1.16 amplitude factor)';...
              'Drift: degree one polynomial fit';...
              'Steps: determined via visual inspection';...
              'All reductions applied before filtering/decimation';...
              'Gravity is not corrected for PCB-Temp effects (see Level1/AUX time series)';...
              'Calibration parameters: see IGETS-IGRAV-CAL-we006-20160700.cal';...
              'Processing scripts can be found at: https://github.com/emenems/igetsDataTools'};

% Set OUTPUT file format (more or less fixed)
channel_names = {'gravity','pressure','tides','atmosphere','polar motion','drift'};
channel_units = {'nm/s^2','hPa','nm/s^2','nm/s^2','nm/s^2','nm/s^2'};
header_offset = 21;
out_precision = {'%10.2f','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f'};
file_format = 'preterna';
nanval = 99999.999; % Flagged NaN values

%% Prepare for loading
% Convert the input starting time and ending time to matlab format 
% suitable for loading
j = 1;
for year = start_time(1):end_time(1)
    if j == 1
        mz = start_time(2);
    else
        mz = 1;
    end
    if year == end_time(1)
        mk = end_time(2);
    else
        mk = 12;
    end
    for m = mz:mk
        time_in(j,1) = year;
        time_in(j,2) = m;
        j = j + 1;
    end
end
time_in(:,3) = 1;
time_in(:,4) = datenum(time_in(:,1),time_in(:,2),time_in(:,3));

%% Load filter
% Load filter in ETERNA modified format (header must be commented using %)
Num = load(filter_file);     
% Stack the filter (ETERNA uses only one half of the repose = mirror the filter)              
Num = vertcat(Num(:,2),flipud(Num(1:end-1,2)));


%% Load data: Gravity
time = [];
data = [];
for m = 1:size(time_in,1)
    % Create output file path + name 
    file_input = fullfile(input_path,...
                    sprintf('%04d',time_in(m,1)),...
                    sprintf('%s%04d%02d%s',input_prefix,...
                    time_in(m,1),time_in(m,2),input_suffix));
    % Load data
    fprintf('Loading data %s\n',file_input);
    [timec,datac] = loadggp('file_in',file_input,'offset',0,...
                            'nanval',99999.999);    
    % Stack data for processing
    time = vertcat(time,timec);
    data = vertcat(data,datac);
    clc
end

%% Apply calibration parameters
% Phase delay
data(:,1) = interp1(time+calib_delay/86400,data(:,1),time); % re-interpolate
% Amplitude factor
data(:,1) = data(:,1)*calib_factor;
   
%% Prepare/load correction time series
fprintf('Loading/Computing corrections \n');
[temp_time,temp_data] = loadtsf(input_tides);
data(:,3) = interp1(temp_time,temp_data,time);

[~,data(:,4)] = getAtmacs(atmacs_url_link_loc,atmacs_url_link_glo,...
                            time,data(:,2),atmcas_admittance);
                        
[~,data(:,5)] = getEOPeffect(pol_lat,pol_lon,time,1.16);

data(:,6) = polyval(drift_fit,time);      
clc

%% Computer residuals + correct anaomalous intervals
temp = data(:,1) - sum(data(:,3:6),2);
data = horzcat(temp,data(:,2:end));clear temp
% Correct gaps/steps using given file
data = correctTimeInterval(time,data,corr_file,1);
% Load the correction file for logfile records
fileid = fopen(corr_file,'r');
in_cell = textscan(fileid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %s','CommentStyle','%'); 
% convert cell aray (standard textscan output) to matrix with double precision
in_corr = horzcat(double(cell2mat(in_cell(1:14))),double(cell2mat(in_cell(15:16))));
fclose(fileid);


% Automatically remove/interpolate missing data.This procedure is identical 
% to the one in 'igets_convert_tsf_to_1sec.m' script.
if fill_missing > 0
    [data,id_time,id_col] = fillnans('time',time,'data',data,...
                                'max_wind',fill_missing); 
end

%% Filter
[timeout,dataout] = mm_filt(time,data,Num,1/1440);

%% Write result
fprintf('Write results \n');
% Use only output time
time_in(time_in(:,4) > datenum(end_time_write),:) = [];
time_in(time_in(:,4) < datenum([start_time_write(1:2),1]),:) = [];
for m = 1:size(time_in,1)  
    % Create output file path + name 
    file_output1 = fullfile(output_path,...
                    sprintf('%04d',time_in(m,1)));
    file_output2 = sprintf('%s%04d%02d%s',output_prefix,...
                    time_in(m,1),time_in(m,2),output_suffix);
    % Check if output folder exist, if not, create it
    if exist(file_output1,'dir')~=7
        mkdir(file_output1)
    end
    % Use data for current month only
    if m == 1
        time1h = transpose(datenum(start_time_write):1/24:(time_in(m+1,4)-1/24));
    elseif m ~= 1 && m ~= size(time_in,1)
        time1h = transpose(time_in(m,4):1/24:(time_in(m+1,4)-1/24));
    else
        time1h = transpose(time_in(m,4):1/24:datenum(end_time_write));
    end
    data1h = interp1(timeout,dataout,time1h);
    
    % Store output file name in header
    header(1,2) = {file_output2};
    % Write monthly data
    writeggp('time',time1h,'data',data1h,'header_offset',header_offset,'header',header,...
          'header_add',header_add,'channels',channel_names,...
          'units',channel_units,'output_file',fullfile(file_output1,file_output2),...
          'out_precision',out_precision,'format',file_format,...
          'nanval',nanval);
    
    % Write monthly logfile
    file_log = sprintf('%s%04d%02d%s',logfile_prefix,...
                    time_in(m,1),time_in(m,2),logfile_suffix);
    % Open logfile+write header
    fid = fopen(fullfile(file_output1,file_log),'w'); 
    fprintf(fid,'Filename 			: %s\n',file_log);
    fprintf(fid,'Station 			: %s\n',header{2,2});
    fprintf(fid,'Instrument 			: %s\n',header{3,2});
    fprintf(fid,'Author 				: %s\n',header{7,2});
    fprintf(fid,'yyyymmdd hhmmss comments\nC*************************************************\n');
    fprintf(fid,'77777777\n');
    
    % Write user defined gaps and steps
    if ~isempty(in_corr)
        for r = 1:length(in_corr(:,1))
            if in_corr(r,1) == 2
                fstring = sprintf('%04d%02d%02d %02d%02d%02.0f to %04d%02d%02d %02d%02d%02.0f NaN data inserted (applied before filtering and decimation)',...
                            in_corr(r,3:14));
            elseif in_corr(r,1) == 3
                fstring = sprintf('%04d%02d%02d %02d%02d%02.0f to %04d%02d%02d %02d%02d%02.0f Data linearly interpolated (applied before filtering and decimation)',...
                            in_corr(r,3:14));
            elseif in_corr(r,1) == 1
                fstring = sprintf('%04d%02d%02d %02d%02d%02.0f Step of %.2f nm/s^2 removed (applied before filtering and decimation)',...
                            in_corr(r,3:8),in_corr(r,16)-in_corr(r,15));
            end
            % Chech if the correction has been applied to current month
            temp1 = datenum(in_corr(r,3:8));
            if (temp1 >= time1h(1) && ... % only current month not appended (head and tail)
                    temp1 <= time1h(end)) 
                fprintf(fid,'%s\n',fstring);
            end
            clear temp1 fstring
        end
    end
    
    % Write which data points have been automatically replaced by 
    % interpolation
    if ~isempty(id_col)
        for r = 1:length(id_time)
            % Create sting allowing arbitrary number of data columns
            fprintf_string = '%04d%02d%02d %02d%02d%02.0f NaN data was automatically (linearly) interpolated after filtering and decimation\n';
            % Write
            otime = datevec(id_time(r));
            if (id_time(r) >= time1h(1) && ... % only current month not appended (head and tail)
                    id_time(r) <= time1h(end)) 
                fprintf(fid,fprintf_string,otime(1),otime(2),otime(3),otime(4),...
                        otime(5),otime(6));  
            end
            clear j fprintf_string otime
        end
    end
    % Close logfile
    fprintf(fid,'99999999\n');
    fclose(fid);
    clear time1h data1h file_output1 file_output2 
end
clc
