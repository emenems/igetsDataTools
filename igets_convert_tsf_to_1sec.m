%% Convert 1 second iGrav tsoft files to ggp
% This script will convert daily tsoft files to ggp/igets monthly files
% Following data structure of INPUT files is required
%     'input_path\input_prefix_YYYY\MMDD\Data_input_prefix_MMDD.tsf
% Following OUTPUT structure will be created
%     'output_path\YYYY\output_prefix+YYYYMM+output_suffix'
% Following steps are carried out:
%   1. Load data
%   2. Stack all daily data within one month into one data matrix
%   3. Check for errors in file = sort according to time + remove
%   ambiguities
%   4. Force to equal sampling! This will insert NaNs into time series
%   where no data is available
%   5. Append (temporarily) data from previous and next month to allow for 
%   interpolation at the edges of time series
%   6. Replace NaNs if such interval is shorter than 10 seconds via linear 
%   interpolation. All replaced values will be noted in the monthly 
%   logfiles.
%   7. Save the results to output file (existence of sub-folder structure
%   will be checked)
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

%% Main SETTINGS
% Process time interval
start_time = [2015 03 05 14 00 00];% e.g., [2015 03 05 14 00 00];
end_time   = [2017 03 06 23 59 59];% e.g., [2017 03 06 23 59 59];
% INPUT File path/name settings (to tsf files)
input_path = 'y:\iGrav\iGrav006 Data'; % year/month/day will be generated automatically
input_prefix = 'iGrav006'; % file name prefix
input_suffix = '.tsf';
% Set which channels should be exported (e.g., gravity and pressure
% channels in INPUT tsf files)
input_channels = [1,2];% e.g, [15] for PCB temperature;
channel_names = {'gravity','pressure'}; % e.g., {'PCB-Temp'} ;
channel_units = {'V','hPa'};% e.g., {'degC'};
% Set what longest time interval of NaNs or missing data should be filled
% with interpolated values 
fill_missing = 10; % (seconds/input time step)

% Set OUTPUT file naming
output_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level1'; % Level1 data
output_prefix = 'IGETS-IGRAV-SEC-we006-'; % e.g. for auxiliary 'IGETS-IGRAV-AUX-we006-';
output_suffix = '00.ggp'; % e.g. 00.aux auxiliary
% Set output logfile prefix (will be created for each month). Same output 
% folder as for data will be used!
logfile_prefix = 'IGETS-IGRAV-STATLOG-we006-'; % e.g., 'IGETS-IGRAV-AUXLOG-we006-'
logfile_suffix = '00.log';

% Set header for OUTPUT file
header =   {'Filename','';... % file name will be appended automatically
            'Station','Wettzell';...
            'Instrument','iGrav006';...
            'N. Latitude (deg)', '49.1449    0.0001 measured';...
            'E. Longitude (deg)','12.8769    0.0001 measured';...
            'Elevation MSL (m)', '609.76     0.3000 measured';...
            'Author','M. Mikolaj (mikolaj@gfz-potsdam.de)'};
% Header footer (OUTPUT)
header_add = {'To get (geoid) height of the sensor add 1.05 m (0.03 measured)';... %. e.g, 'The PCB-Temp is measured inside the iGrav head'
              'Always check the STATLOG files for details on time series quality';... % e.g., 'Variations of PCB-Temp can affect the measured gravity'
              'Processing scripts can be found at: https://github.com/emenems/igetsDataTools'};

% Set OUTPUT file format (more or less fixed)
header_offset = 21; % fixed for IGETS
out_precision = {'%10.6f','%10.4f'}; % e.g., {'%10.4f'}
file_format = 'preterna';
nanval = 99999.999; % Flagged NaN values

%% Prepare for loading
% Convert the input starting time and ending time to matlab format and
% create a time vector. Time step is ONE day = one file per day written by
% iGrav or SG  
time_in(:,7) = [datenum(start_time(1:3)):1:datenum(end_time(1:3))]';      
time_in(:,1:6) = datevec(time_in(:,7));   
time_pattern = time_in(:,1)*100+time_in(:,2);

%% Load data
% Run loop for each month. First find number of months within given
% interval
months = unique(time_pattern);
count = 0; % count number of loaded files
for m = 1:length(months)
    % Create output file folder
    file_output1 = fullfile(output_path,...
                    sprintf('%04d',floor(months(m)/100)));
    % Check if output folder exist, if not, create it
    if exist(file_output1,'dir')~=7
        mkdir(file_output1)
    end
    % Create output logfile name
    file_log = sprintf('%s%06d%s',logfile_prefix,...
                    months(m),logfile_suffix);
    if ~isempty(logfile_prefix)
        % Open logfile+write header
        fid = fopen(fullfile(file_output1,file_log),'w'); 
        fprintf(fid,'Filename 			: %s\n',file_log);
        fprintf(fid,'Station 			: %s\n',header{2,2});
        fprintf(fid,'Instrument 			: %s\n',header{3,2});
        fprintf(fid,'Author 				: %s\n',header{7,2});
        fprintf(fid,'yyyymmdd hhmmss comments\nC*************************************************\n');
        fprintf(fid,'77777777\n');
    end
    % Declare output variables
    time = [];
    data = [];
    % For each month find days within this month and load corresponding
    % data. Following command will find indices in the 'time_in' variable
    index = find(time_pattern == months(m));
    for i = 1:length(index)                          
        try
            % create input (not zip) file name = file path + file prefix + date + .tsf
            file_name = fullfile(input_path,...
                        sprintf('%s_%04d',input_prefix,time_in(index(i),1)),...
                        sprintf('%02d%02d',time_in(index(i),2),time_in(index(i),3)),...
                        sprintf('Data_%s_%02d%02d%s',input_prefix,time_in(index(i),2),...
                        time_in(index(i),3),input_suffix));
            fprintf('Loading %s\n',file_name);
            % load file and store to temporary variables          
            [ttime,tdata] = loadtsf(file_name); 
            count = count + 1;
            % Remove data other than gravity and pressure
            tdata = tdata(:,input_channels);
        catch error_message
            fprintf('Error: %s',error_message.message);
            ttime = [];
            tdata = [];
            count = count + 1;
        end
        % stack the temporary variable on already loaded ones 
        time = vertcat(time,ttime);             
        data = vertcat(data,tdata);
        clear ttime tdata
        clc
    end
    
    %% Prepare data for processing
    % Create output name(2)
    file_output2 = sprintf('%s%06d%s',output_prefix,...
                    months(m),output_suffix);
    fprintf('Preparing %s\n',file_output2);
    % Check the sampling for errors or missing data 
    % Sort the time vector and data
    [time_temp,temp_index] = sort(time,1);
    data_temp = data(temp_index,:);
    % Use unique time values only
    [time,temp_index] = unique(time_temp);
    data = data_temp(temp_index,:);
    % To allow for interpolation at the 
    % end and start of the file, load data from first and last day of the 
    % next and previous day respectively.
    if length(data) ~= (length(index)*86400)
        % Load head and tail
        temp_time1 = [];temp_data1 = [];
        temp_time2 = [];temp_data2 = [];
        if m ~= 1
            temp_file = fullfile(input_path,...
                        sprintf('%s_%04d',input_prefix,time_in(index(1)-1,1)),...
                        sprintf('%02d%02d',time_in(index(1)-1,2),time_in(index(1)-1,3)),...
                        sprintf('Data_%s_%02d%02d%s',input_prefix,time_in(index(1)-1,2),...
                        time_in(index(1)-1,3),input_suffix));
            if exist(temp_file,'file') == 2
                [temp_time1,temp_data1] = loadtsf(temp_file); 
                % Use only first X seconds
                temp_data1 = temp_data1(end-fill_missing-1:end,input_channels);
                temp_time1 = temp_time1(end-fill_missing-1:end,:);
            end
            clear temp_file
        end
        if m ~= length(months)
            temp_file = fullfile(input_path,...
                        sprintf('%s_%04d',input_prefix,time_in(index(end)+1,1)),...
                        sprintf('%02d%02d',time_in(index(end)+1,2),time_in(index(end)+1,3)),...
                        sprintf('Data_%s_%02d%02d%s',input_prefix,time_in(index(end)+1,2),...
                        time_in(index(end)+1,3),input_suffix));
            if exist(temp_file,'file') == 2
                [temp_time2,temp_data2] = loadtsf(temp_file);
                % Use only first X seconds
                temp_data2 = temp_data2(1:fill_missing+1,input_channels);
                temp_time2 = temp_time2(1:fill_missing+1,:);
            end
        end
        % Re-sample to equal data. The appended values from both sides will
        % be removed later
        [time,data] = findTimeStep(vertcat(temp_time1,time,temp_time2),...
                                   vertcat(temp_data1,data,temp_data2),...
                                   1/86400);
    end
    
    %% Fill missing data
    fprintf('Filling missing/NaN data %s\n',file_output2);
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
            clear j fprintf_string
        end
    end
    clear r id_time id_cor  
    
    %% Write output
    % Remove data out of current month (appended for interpolation only)
    data(time < time_in(index(1),7),:) = [];
    time(time < time_in(index(1),7),:) = [];                        
    data(time > datenum([time_in(index(end),1:3),23,59,59])+1e-6,:) = []; 
    time(time > datenum([time_in(index(end),1:3),23,59,59])+1e-6,:) = []; 
    % Just in case, remove also data out of requested range
    data(time>datenum(end_time),:) = [];
    time(time>datenum(end_time),:) = [];
    data(time<datenum(start_time),:) = [];
    time(time<datenum(start_time),:) = [];
    
    % Store output file name in header
    header(1,2) = {file_output2};
    fprintf('Write data to %s\n',file_output2);
    writeggp('time',time,'data',data,'header_offset',header_offset,'header',header,...
          'header_add',header_add,'channels',channel_names,...
          'units',channel_units,'output_file',fullfile(file_output1,file_output2),...
          'out_precision',out_precision,'format',file_format,...
          'nanval',nanval);
    clear file_output1 file_output2
    % Close logfile
    if ~isempty(logfile_prefix)
        fprintf(fid,'99999999\n');
        fclose(fid);
    end
    clc
end

%% End
rmpath('f:\mikolaj\code\libraries\matlab_octave_library')
