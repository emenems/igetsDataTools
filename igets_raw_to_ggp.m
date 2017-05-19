function igets_raw_to_ggp(varargin)
%IGETS_RAW_TO_GGP Convert 1 second iGrav/SG tsoft files to ggp/igets
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
% Tested on Matlab R2015b (preferred) and Octave 4.2.1 (rather slow)
%
% INPUTS:
%  'start' 		  ... starting time 
%                       Example:  [2015 03 05 14 00 00]
%  'stop'         ... end time 
%                       Example:  [2017 03 06 23 59 59];
%  'input_path'   ... input path for loading raw tsoft files 
%                       Example:  'y:\iGrav\iGrav006 Data'
%  'input_prefix' ... input tsoft file name prefix
%                       Example:  'iGrav006'
%  'input_suffix' ... input tsoft file name suffix
%                       Example:  '.tsf'
%  'input_channels'.. channels/columns from input file
%                       Example:  [1,2]
%  'channel_names'... channel names used for output file
%                       Example:  {'gravity','pressure'}
%  'channel_units'... channel units used for output file
%                       Example:  {'V','hPa'}
%  'output_path'  ... output (monthly) data folder
%                       Example:  'f:\we006\Level1';
%  'output_prefix'... output (monthly) data prefix
%                       Example:  'IGETS-IGRAV-SEC-we006-'
%  'output_suffix'... output (monthly) data suffix
%                       Example:  '00.ggp'
%  'logfile_prefix'.. logfile prefix (same folder as 'output_path')
%                       Example:  'IGETS-IGRAV-AUXLOG-we006-'
%  'logfile_suffix'.. logfile suffix
%                       Example:  '00.log';
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
%                       Example:  99999.999
%  'fill_missing' ... Set what longest time interval of NaNs or missing
%                       data should be filled with interpolated values
%                       Example: 10
%  'header_offset'... row offset in header
%                       Example: 21
%
%                                                    M.Mikolaj
%                                                    mikolaj@gfz-potsdam.de

%% Read user input
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
            case 'out_precision'
                out_precision = varargin{in+1};
            case 'nanval'
                nanval = varargin{in+1};
            case 'fill_missing'
                fill_missing = varargin{in+1};
            case 'input_channels'
                input_channels = varargin{in+1};
            case 'channel_names'
                channel_names = varargin{in+1};
            case 'channel_units'
                channel_units = varargin{in+1};
            case 'header_offset'
                header_offset = varargin{in+1};
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
header_add(end+1) = {'Processing scripts can be found at: https://github.com/emenems/igetsDataTools'};

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
        fprintf(fid,'Station 			: %s\n',station);
        fprintf(fid,'Instrument 			: %s\n',instrument);
        fprintf(fid,'Author 				: %s\n',author);
        fprintf(fid,'All time stamps are approximate\n');
        fprintf(fid,'Processing scripts available at: https://github.com/emenems/igetsDataTools\n');
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

end % function
