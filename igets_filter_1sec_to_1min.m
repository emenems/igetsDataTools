%% Convert 1 second igets data to 1 minute data
% This script will convert one second monthly ggp/igets files to one minute 
% ggp/igets monthly files.
% Use the 'igets_convert_tsf_to_1sec.m' first to get the monthly 1 second 
% data.
% Following data structure of INPUT files is required
%     'input_path\YYYY\input_prefix+YYYYMM+input_suffix'
% Following OUTPUT structure will be created
%     'output_path\YYYY\output_prefix+YYYYMM+output_suffix'
% Following steps are carried out:
%   1. Load data for current month
%   2. Load data of the previous and next month to allow for filtering at 
%      the edges
%   3 If required, the input (1 second) time series will be automatically 
%      inspected for missing/NaN data and specified interval can be filled 
%      by interpolated values. This should be normally done while creating
%      1 second data.
%   4. Filter the loaded time series (if filter selected). 
%   5. Re-sample/decimate the time series to 1 min resolution
%   6. Write the result
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

%% Main settings
% Process time interval
start_time = [2015 03 05 14 00 00];% e.g., [2015 03 05 14 00 00];
end_time   = [2017 03 06 23 59 59];% e.g., [2017 03 06 23 59 59];
% INPUT File path/name settings
input_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1'; % year/month/day will be generated automatically
input_prefix = 'IGETS-IGRAV-SEC-we006-'; % file name prefix
input_suffix = '00.ggp';
% Set which channels should be loaded&exported (e.g., gravity and pressure))
input_channels = [1,2];
% Maximum time interval to be interpolated in case of missing or NaN data. Set 
% to 0 for not interpolation. Such interpolation would be done prior
% filtering (normally done within processing of 1 second data).
fill_missing = 0; % seconds (=in input units)
% Set filter input filter. This file must be modified: all header lines start 
% '%' and only half of the impulse response is given (second will be created via
% flipping). 
filter_file = fullfile('data','g1s1md.gwr');

% Set OUTPUT file naming
output_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1';
output_prefix = 'IGETS-IGRAV-MIN-we006-';
output_suffix = '00.ggp';
% Set output logfile. One logfile documenting the steps will be written (if not 
% set to []). This is just for your info, not for IGETS!
logfile = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1\IGETS-IGRAV-MIN-STATLOG-we006_ALL.log';

% Set header for OUTPUT file
header   = {'Filename','';... % file name will be appended automatically
            'Station','Wettzell';...
            'Instrument','iGrav006';...
            'N. Latitude (deg)', '49.1449    0.0001 measured';...
            'E. Longitude (deg)','12.8769    0.0001 measured';...
            'Elevation MSL (m)', '609.76     0.3000 measured';...
            'Author','M. Mikolaj (mikolaj@gfz-potsdam.de)'};      
% Header footer (OUTPUT)
header_add = {'To get (geoid) height of the sensor add 1.05 m (0.03 measured)';...
              'Always check the STATLOG files for details on time series quality';...
              'Processing scripts can be found at: https://github.com/emenems/igetsDataTools'};

% Set OUTPUT file format (more or less fixed)
channel_names = {'gravity','pressure'};
channel_units = {'V','hPa'};
header_offset = 21;
out_precision = {'%10.6f','%10.4f'};
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
if ~isempty(logfile)
    fid = fopen(logfile,'w');
end

%% Load filter
if ~isempty(filter_file)
    % Get the name for logfile and header of the ouput file_in_loadpath
    [~,filter_name,filter_ext] = fileparts(filter_file);  
    header_add(end+1) = {;...
                  sprintf('Data filtered using: %s.%s filter',filter_name,filter_ext)};
    % Load filter in ETERNA modified format (header must be commented using %)

    Num = load(filter_file);     
    % Stack the filter (ETERNA uses only one half of the repose = mirror the filter)              
    Num = vertcat(Num(:,2),flipud(Num(1:end-1,2)));
    if ~isempty(logfile)
        fprintf(fid,'Filter file loaded: %s',filter_file);
    end
else
    if ~isempty(logfile)
        fprintf(fid,'No filter loaded: %s',filter_file);
    end    
end

%% Load data
for m = 1:size(time_in,1)
    % Create input file name of the current month
    file_input_cur = fullfile(input_path,...
                        sprintf('%04d%02d',time_in(m,1)),...
                        sprintf('%s%04d%02d%s',input_prefix,...
                        time_in(m,1),time_in(m,2),output_suffix));
    
    % Load/get current month data. For the first month load the file, otherwise 
    % use the previously loaded for following month.
    if m == 1 
        fprintf('Loading data %s\n',file_input_cur);
        [timec,datac] = loadggp('file_in',file_input_cur,'offset',0,...
                            'nanval',99999.999);    
        datac = datac(:,input_channels);  
        % Previous month data 
        datap = [];
        timep = [];        
    else 
        % Previous
        timep = timec;
        datap = datac;
        % Current
        timec = timef;
        datac = dataf;
    end    
    % Load/get next month data (except for the last month)
    if m~=size(time_in,1)
        file_input_fol = fullfile(input_path,...
                    sprintf('%04d%02d',time_in(m+1,1)),...
                    sprintf('%s%04d%02d%s',input_prefix,...
                    time_in(m+1,1),time_in(m+1,2),output_suffix)); 
        fprintf('Loading data %s\n',file_input_fol);
        [timef,dataf] = loadggp('file_in',file_input_fol,'offset',0,...
                            'nanval',99999.999);   
        dataf = dataf(:,input_channels);                     
    else
        timef = [];
        dataf = [];
    end 
    
    % Create output file path + name 
    file_output1 = fullfile(output_path,...
                    sprintf('%04d',time_in(m,1)));
    file_output2 = sprintf('%s%04d%02d%s',output_prefix,...
                    time_in(m,1),time_in(m,2),output_suffix);
    % Check if output folder exist, if not, create it
    if exist(file_output1,'dir')~=7
        mkdir(file_output1)
    end
    
    % Stack data for processing
    fprintf('Preparing data for processing %s\n',file_output2);
    time = [];
    data = [];
    if ~isempty(timep)
        time = vertcat(timep(end-length(Num*4):end,:),...
                timec);
        data = vertcat(datap(end-length(Num*4):end,:),...
                datac);
    end
    if ~isempty(timef)
        time = vertcat(time,...
                timef(1:length(Num*4),:));
        data = vertcat(data,...
                dataf(1:length(Num*4),:));
    end
            
    % If required interpolate NaNs. This procedure is identical to the one
    % in 'igets_convert_tsf_to_1sec.m' script.
    if fill_missing > 0
        [data,id_time,id_col] = fillnans('time',time,'data',data,...
                                    'max_wind',fill_missing);
        % Write which data points have been replaced by interpolation
        if ~isempty(id_col)
            for r = 1:length(id_time)
                % Create sting allowing arbitrary number of data columns
                fprintf_string = '%04d%02d%02d %02d%02d%02.0f Missing or NaN data was linearly interpolated (';
                for j = 1:size(data,2)
                    if id_col(r,j)
                        if strcmp(fprintf_string(end),'s')
                            fprintf_string = [fprintf_string,' & %s'];
                        else
                            fprintf_string = [fprintf_string,'%s'];
                        end
                    end
                end
                fprintf_string = [fprintf_string,')\n'];
                % Write
                otime = datevec(id_time(r));
                if (id_time(r) >= time_in(index(1),7) && ... % only current month not appended (head and tail)
                        id_time(r) <= datenum([time_in(index(end),1:3),23,59,59])+1e-6) 
                    fprintf(fid,fprintf_string,otime(1),otime(2),otime(3),otime(4),...
                            otime(5),otime(6),channel_names{logical(id_col(r,:))});  
                end
                clear j fprintf_string otime
            end
        end
        clear r id_time id_cor  
    end
    
    %% Filter data + re-sample
    if ~isempty(filter_file)
        fprintf(fid,'Filtering data %s\n',file_output2);
        [timeout,dataout] = mm_filt(time,data,Num,1/86400);
    else
        timeout = time;
        dataout = data;
    end
    % Resample
    if m == 1
        time = transpose(datenum(start_time):1/1440:(time_in(m+1,4)-1/86400));
    elseif m ~= 1 && m ~= size(time_in,1)
        time = transpose(time_in(m,4):1/1440:(time_in(m+1,4)-1/86400));
    else
        time = transpose(time_in(m,4):1/1440:datenum(end_time));
    end
    data = interp1(timeout,dataout,time,'linear');
    
    %% Write
    % Store output file name in header
    header(1,2) = {file_output2};
    fprintf(fid,'Write data to %s\n',file_output2);
    writeggp('time',time,'data',data,'header_offset',header_offset,'header',header,...
          'header_add',header_add,'channels',channel_names,...
          'units',channel_units,'output_file',fullfile(file_output1,file_output2),...
          'out_precision',out_precision,'format',file_format,...
          'nanval',nanval);
    clear file_output1 file_output2 timeout dataout
    clc
end
    
%% End
% Close logfile
if ~isempty(logfile)
    fclose(fid);
end
rmpath('f:\mikolaj\code\libraries\matlab_octave_library')