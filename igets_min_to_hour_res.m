function igets_min_to_hour_res(varargin)
%IGETS_MIN_TO_HOUR_RES Convert 1 minute data to 1 hour residuals
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
% INPUTS:
%  'start' 		  ... starting time 
%                       Example:  [2015 03 05 14 00 00]
%  'stop'         ... end time 
%                       Example:  [2017 03 06 23 59 59];
%  'input_path'   ... input path for loading 1second ggp/igets data
%                       Example:  'f:\we006\Level1'
%  'input_prefix' ... input igets file name prefix
%                       Example:  'IGETS-IGRAV-MIN-we006-'
%  'input_suffix' ... input igets file name suffix
%                       Example:  '00.ggp'
%  'output_path'  ... output (monthly) data folder
%                       Example:  'f:\we006\Level3';
%  'output_prefix'... output (monthly) data prefix
%                       Example:  'IGETS-IGRAV-HOURES-we006-';
%  'output_suffix'... output (monthly) data suffix
%                       Example:  'r1.ggp'
%  'logfile_prefix'.. logfile prefix (same folder as 'output_path')
%                       Example:  'IGETS-IGRAV-RESHOUR-we006-'
%  'logfile_suffix'.. logfile suffix
%                       Example:  'r1.log';
%  'file_format'  ... 'preterna' or 'eterna'
%  'instrument'   ... name of the instrument for igets file header
%                       Example: 'iGrav006'
%  'station'      ... name of the station for igets file header
%                       Example:  'Wettzell'
%  'latitude'     ... station latitude + accuracy
%                       Example:  '49.1449    0.0001 measured'
%  'longitude'    ... station longitude + accuracy
%                       Example:  '12.8769    0.0001 measured'
%  'height'       ... station longitude + accuracy
%                       Example:  '609.76     0.3000 measured'
%  'author'       ... file author
%                       Example:  'M. Mikolaj (mikolaj@gfz-potsdam.de)'
%  'header_add'   ... add text to igets file header
%                       Example:  {'Sensor height 1.05 m (0.03 measured)'}
%  'out_precision'... output precision
%                       Example:  {'%10.6f','%10.4f'}
%  'nanval'       ... flagged NaN values
%                       Example:  99999.999 (=default)
%  'fill_missing' ... Set what longest time interval of NaNs or missing
%                       data should be filled with interpolated values
%                       Example: 10
%  'filter_file'  ... file used to filter data prior interpolation
%                       Example: fullfile('data','g1m1h.nlf')
%  'header_offset'... optional row offset in header
%                       Example: 21 (=default)
%  'input_tides'  ... input tides file in tsf format where the 'input_tides_chan'
%                       contains tidal signal. 
%                       Example: fullfile('data','tides_file_example.tsf')
%  'input_tides_chan' optional tsf channel containing tidal effect on gravity
%						Example: 1 (=default)
%  'tides_interp' ... optional interpolation method for tides re-sampling 
%						Example: 'linear' (=default)
%  'atmacs_loc'   ... urls to local Atmcas contribution
%                       Example: {'http://atmacs.bkg.bund.de/data/results/lm/we_lm2_12km_19deg.grav',...
%                                 'http://atmacs.bkg.bund.de/data/results/iconeu/we_iconeu_70km.grav'};
%  'atmacs_glo'   ... urls to local Atmcas contribution
%                       Example: {'http://atmacs.bkg.bund.de/data/results/icon4lm/we_icon384_19deg.grav',...
%                                 'http://atmacs.bkg.bund.de/data/results/icon/we_icon384_20deg.grav'};
%  'corr_file'    ... optional correction file with steps and gaps,
%                       correctTimeInterval.m will be used for this purpose
%                       Example: fullfile('data','corr_file_example.txt')
% Channel names, units, input columns and output precision are fixed
%
%                                                    M.Mikolaj
%                                                    mikolaj@gfz-potsdam.de

%% Read user input
% Default values
input_tides_chan = 1;
tides_interp = 'linear';
header_offset = 21;
fill_missing = 0;
nanval = 99999.999;
header_add = [];

% First check if correct number of input arguments
if nargin > 2 && mod(nargin,2) == 0
    % Count input parameters
    in = 1;
    % Try to find input parameters
    while in < nargin
        % Switch between function parameters
        switch varargin{in}
            case 'start'
                start_time = varargin{in+1};
            case 'stop'        
                end_time = varargin{in+1};
            case 'input_path'
                input_path = varargin{in+1};
            case 'input_prefix'
                input_prefix = varargin{in+1};
            case 'input_suffix'
                input_suffix = varargin{in+1};
            case 'output_path'
                output_path = varargin{in+1};
            case 'output_prefix'
                output_prefix = varargin{in+1};
            case 'output_suffix'
                output_suffix = varargin{in+1};
            case 'logfile_prefix'
                logfile_prefix = varargin{in+1};
            case 'logfile_suffix'
                logfile_suffix = varargin{in+1};
            case 'file_format'
                file_format = varargin{in+1};
            case 'instrument'
                instrument = varargin{in+1};
            case 'station'
                station = varargin{in+1};
            case 'header_add'
                header_add = varargin{in+1};
            case 'latitude'
                latitude = varargin{in+1};
            case 'longitude'
                longitude = varargin{in+1};
            case 'height'
                height = varargin{in+1};
            case 'author'
                author = varargin{in+1};
            case 'nanval'
                nanval = varargin{in+1};
            case 'fill_missing'
                fill_missing = varargin{in+1};
            case 'header_offset'
                header_offset = varargin{in+1};
            case 'input_tides'
                input_tides = varargin{in+1};
            case 'atmacs_loc'
                atmacs_url_link_loc = varargin{in+1};
            case 'atmacs_glo'
                atmacs_url_link_glo = varargin{in+1};
            case 'admittance'
                atmcas_admittance = varargin{in+1};
            case 'drift'
                drift_fit = varargin{in+1};
            case 'calib_delay'
                calib_delay = varargin{in+1};
            case 'calib_factor'
                calib_factor = varargin{in+1};
            case 'corr_file'
                corr_file = varargin{in+1};
            case 'filter_file'
                filter_file = varargin{in+1};
			case 'input_tides_chan'
                input_tides_chan = varargin{in+1};
			case 'tides_interp'
                tides_interp = varargin{in+1};
        end
        % Increase by 2 as parameters are in pairs!
        in = in + 2;
    end
elseif nargin > 0 && mod(nargin,2) ~= 0
    error('Set even number of input parameters')
end

%% Create header for output files
header =   {'Filename','';... % file name will be appended automatically
            'Station',station;...
            'Instrument',instrument;...
            'N. Latitude (deg)', latitude;...
            'E. Longitude (deg)',longitude;...
            'Elevation MSL (m)', height;...
            'Author',author};
header_add(end+1) = {'Tools for estimation of large-scale effects on gravity can be found at: https://github.com/emenems/mGlobe'};
header_add(end+1) = {'Processing scripts can be found at: https://github.com/emenems/igetsDataTools'};

% Polar motion effect (computed using getEOPeffect.m):
temp = strsplit(latitude);
pol_lat = str2double(temp{1});
temp = strsplit(longitude);
pol_lon = str2double(temp{1});

channel_names = {'gravity','pressure','tides','atmosphere','polar motion','drift'};
channel_units = {'nm/s^2','hPa','nm/s^2','nm/s^2','nm/s^2','nm/s^2'};
out_precision = {'%10.2f','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f'};

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
                            'nanval',nanval);    
    % Stack data for processing
    time = vertcat(time,timec);
    data = vertcat(data,datac);
    clc
end

%% Apply calibration parameters
% Phase delay
data(:,1) = interp1(time+calib_delay/86400,data(:,1),time,tides_interp); % re-interpolate
% Amplitude factor
data(:,1) = data(:,1)*calib_factor;
   
%% Prepare/load correction time series
fprintf('Loading/Computing corrections \n');
[temp_time,temp_data] = loadtsf(input_tides);
data(:,3) = interp1(temp_time,temp_data(:,input_tides_chan),time);

[~,data(:,4)] = getAtmacs(atmacs_url_link_loc,atmacs_url_link_glo,...
                            time,data(:,2),atmcas_admittance);
                        
[~,data(:,5)] = getEOPeffect(pol_lat,pol_lon,time,1.16);

data(:,6) = polyval(drift_fit,time);      
clc

%% Computer residuals + correct anaomalous intervals
temp = data(:,1) - sum(data(:,3:6),2);
data = horzcat(temp,data(:,2:end));clear temp
% Correct gaps/steps using given file
if ~isempty(corr_file)
	data = correctTimeInterval(time,data,corr_file,1);
	% Load the correction file for logfile records
	fileid = fopen(corr_file,'r');
	in_cell = textscan(fileid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %s','CommentStyle','%'); 
	% convert cell aray (standard textscan output) to matrix with double precision
	in_corr = horzcat(double(cell2mat(in_cell(1:14))),double(cell2mat(in_cell(15:16))));
	fclose(fileid);
else
	in_corr = [];
end


% Automatically remove/interpolate missing data.This procedure is identical 
% to the one in 'igets_raw_to_ggp.m' script.
if fill_missing > 0
    [data,id_time,id_col] = fillnans('time',time,'data',data,...
                                'max_wind',fill_missing); 
end

%% Filter
[timeout,dataout] = mm_filt(time,data,Num,1/1440);

%% Write result
fprintf('Write results \n');
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
        % Use only time interval not affected by filtering
        time1h = transpose(datenum(start_time):1/24:(time_in(m+1,4)-1/24));
        if length(time1h) > floor(length(Num)/120)+1
            time1h(1:floor(length(Num)/120)+1,:) = []; % /120 = 2*60 = half of filter * convert to hours
        else
            time1h = [];
        end
    elseif m ~= 1 && m ~= size(time_in,1)
        time1h = transpose(time_in(m,4):1/24:(time_in(m+1,4)-1/24));
    else
        time1h = transpose(time_in(m,4):1/24:datenum(end_time));
        if length(time1h) > floor(length(Num)/120)+1
            time1h(end-floor(length(Num)/120)-1:end,:) = [];
        else
            time1h = [];
        end
    end
    if ~isempty(time1h)
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
        fprintf(fid,'Station 			: %s\n',station);
        fprintf(fid,'Instrument 			: %s\n',instrument);
        fprintf(fid,'Author 				: %s\n',author);
        fprintf(fid,'All time stamps are approximate\n');
        fprintf(fid,'Processing scripts available at: https://github.com/emenems/igetsDataTools\n');
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
                fprintf_string = '%04d%02d%02d %02d%02d%02.0f NaN data automatically replaced by (linearly) interpolated value (after filtering and decimation)\n';
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
    end
    clear time1h data1h file_output1 file_output2 
end
clc
