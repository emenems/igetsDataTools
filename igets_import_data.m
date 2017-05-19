function igets_import_data(varargin)
%IGETS_IMPORT_DATA Convert data stored in one file to monthly IGETS format
% Use this scrip to convert tsf or mat (structure array) data to monthly
% ggp file as required by IGETS
% Following OUTPUT structure will be used
%     'output_path\YYYY\output_prefix+YYYYMM+output_suffix'
% Following steps are carried out:
%   1. load data
%   2. resample/decimate if selected
%   3. replace NaNs by inerpolated values if selected
%   4. (or 2.) apply filter if selected (if filter_after = 0)
%   5. Save monthly data
%
% INPUTS:
%  'start' 		  ... (optional) starting time of the output/ggp file
%                       Example:  [2015 03 05 14 00 00]
%  'stop'         ... (optional) end time time of the output/ggp file
%                       Example:  [2017 03 06 23 59 59];
%  'input_file'   ... full file name of the input
%                       Example:  'y:\iGrav\iGrav006\Whole_time.tsf'
%  'input_channels'.. channels/columns from input file to be exported
%                       Example:  [1,2]
%  'channel_names'... channel names used for output file
%                       Example:  {'gravity','pressure'}
%  'channel_units'... channel units used for output file
%                       Example:  {'V','hPa'}
%  'out_precision'... output precision
%                       Example:  {'%10.6f','%10.4f'}
%  'output_path'  ... output (monthly) data folder
%                       Example:  'f:\we006\Level1';
%  'output_prefix'... output (monthly) data prefix
%                       Example:  'IGETS-IGRAV-AUX-we006-'
%  'output_suffix'... output (monthly) data suffix
%                       Example:  '00.ggp'
%  'resample'     ... resample/decimate to 'sec','min', or 'hour'. Do not
%                       set (or []) for no resampling. Simple linear
%                       interpolation will be used if selected.
%                       Example: 'hour'
%  'fill_missing' ... Set what longest time interval of NaNs or missing
%                       data should be filled with interpolated values
%                       Use same time resolution as input or resampled time 
%                       series (if selected), e.g., if time(2)-time(1) = 
%                       10 seconds (10/86400 days)=> fill_missing must be 
%                       also 10
%                       Example: 10 (set to 0 for not filling)
%  'filter_file'  ... (optional) filter data if needed. The filter must be
%                       in the same time units as input/resampled time! 
%                       Use 'filter_after' to switch between filtering
%                       applied before (0) of after resampling and filling
%                       missing data (1)
%                       Example: fullfile('data','g1s1md.gwr')
%  'filter_after' ...  (optional) switch between filtering applied before 
%                       (0) of after resampling and filling missing data 
%                       (1 = default). Filling missing data will be always
%                       carried out after resampling!
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
%
%                                                    M.Mikolaj
%                                                    mikolaj@gfz-potsdam.de

%% Set default values
input_file = '';
nanval = 99999.999;
time = [];
data = [];
start_time = [];
end_time = [];
input_channels = 1;
filter_file = [];
fill_missing = 0;
resample = [];
filter_after = 1;
header_add = [];
latitude = '';
longitude = '';
height = '';
author = '';
station = '';
instrument = '';
out_precision = {'10.2f'};
channel_names = '';
channel_units = '';
header_offset = 21;
file_format = 'preterna';

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
            case 'input_file'
                input_file = varargin{in+1};
            case 'output_path'
                output_path = varargin{in+1};
            case 'output_prefix'
                output_prefix = varargin{in+1};
            case 'output_suffix'
                output_suffix = varargin{in+1};
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
            case 'resample'
                resample = varargin{in+1};
            case 'filter_file'
                filter_file = varargin{in+1};
            case 'filter_after'
                filter_after = varargin{in+1};
        end
        % Increase by 2 as parameters are in pairs!
        in = in + 2;
    end
elseif nargin > 0 && mod(nargin,2) ~= 0
    error('Set even number of input parameters')
end

%% Load input file
fprintf('Loading input file\n');
% Get file suffix to switch between input formats format
[~,~,file_suffix] = fileparts(input_file);
if strcmp(file_suffix,'.ggp')
    [time,data]=loadggp('file_in',input_file,...
                        'offset',0,'nanval',nanval);
elseif strcmp(file_suffix,'.tsf')
    [time,data] = loadtsf(input_file);
elseif strcmp(file_suffix,'.mat')
    temp = importdata(input_file);
    time = temp.time;
    data = temp.data;
    clear temp;
end
% Select input channels
data = data(:,input_channels);

% Load filter if on input
if ~isempty(filter_file)
    Num = load(filter_file);     
    % Stack the filter (ETERNA uses only one half of the repose = mirror the filter)              
    Num = vertcat(Num(:,2),flipud(Num(1:end-1,2)));
end

%% Prepare data
fprintf('Preparing input data\n');
% Cut time intervals out of requested range
if ~isempty(start_time)
    data(time<datenum(start_time),:) = [];
    time(time<datenum(start_time),:) = [];
else
    start_time = datevec(time(1));
end
if ~isempty(end_time)
    data(time>datenum(end_time),:) = [];
    time(time>datenum(end_time),:) = [];
else
    end_time = datevec(time(end));
end
% Input time resolution. Round to max 1 second precision
delta_t = round((time(2) - time(1))*86400)/86400;
% Filter if required
if ~isempty(filter_file) && filter_after == 0
    [timef,dataf] = mm_filt(time,data,Num,delta_t);
end
% Resample if needed
if ~isempty(resample)
    if strcmp(resample,'sec')
        delta_t = 1/86400;
    elseif strcmp(resample,'min')
        delta_t = 1/1440;
    elseif strcmp(resample,'hour')
        delta_t = 1/24;
    end
    if exist('timef','var') == 1
        data = interp1(timef,dataf,transpose(time(1):delta_t:time(end)));
    else
        data = interp1(time,data,transpose(time(1):delta_t:time(end)));
    end
    time = transpose(time(1):delta_t:time(end));
end
% If required interpolate/replace NaNs. 
if fill_missing > 0
    [data,~,~] = fillnans('time',time,'data',data,...
                                'max_wind',fill_missing);
end
% Filter if required
if ~isempty(filter_file) && filter_after == 1
    [time,data] = mm_filt(time,data,Num,delta_t);
end

%% Create header for output files
header =   {'Filename','';... % file name will be appended automatically
            'Station',station;...
            'Instrument',instrument;...
            'N. Latitude (deg)', latitude;...
            'E. Longitude (deg)',longitude;...
            'Elevation MSL (m)', height;...
            'Author',author};
header_add(end+1) = {'Processing scripts available at: https://github.com/emenems/igetsDataTools'};

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

%% Write 
% Continue only if some data was loaded
if ~isempty(time) 
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
        % Cut and write
        if m ~= size(time_in,1)
            timec = time(time>=time_in(m,4) & time<time_in(m+1,4));
            datac = data(time>=time_in(m,4) & time<time_in(m+1,4),:);
        else
            timec = time(time>=time_in(m,4) & time<=datenum(end_time));
            datac = data(time>=time_in(m,4) & time<=datenum(end_time),:);
        end
        fprintf('Writing %s\n',file_output2);
        header(1,2) = {file_output2};
        writeggp('time',timec,'data',datac,'header_offset',header_offset,'header',header,...
          'header_add',header_add,'channels',channel_names,...
          'units',channel_units,'output_file',fullfile(file_output1,file_output2),...
          'out_precision',out_precision,'format',file_format,...
          'nanval',nanval);
        clear timec datac file_output1 file_output2
        clc
    end  
else
    disp('No data loaded (or ready for processing)');
end

end % function