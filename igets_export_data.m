function [time,data] = igets_export_data(varargin)
%IGETS_EXPORT_DATA Load monthly second/minite/hour data +export to one file
% This script will load one second, minute or hour monthly ggp/igets files 
% to one data matrix that will be then written to a ggp/igets, tsf or mat
% file. 
% Following data structure of INPUT files is required
%     'input_path\YYYY\input_prefix+YYYYMM+input_suffix'
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
%  'output_file'  ... full file name of the output 
%                       Example:  'f:\we006\Level3\Whole_time.tsf';
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
%  'out_precision'... output precision
%                       Example:  '%10.2f'
%  'nanval_in'    ... flagged NaN values in input files
%                       Example:  99999.999
%  'nanval_out'   ... flagged NaN values in input files (for ggp only)
%                       Example:  99999.99
% Channel names, units, input columns and output precision are fixed
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
            case 'output_file'
                output_file = varargin{in+1};
            case 'file_format'
                file_format = varargin{in+1};
            case 'instrument'
                instrument = varargin{in+1};
            case 'station'
                station = varargin{in+1};
            case 'latitude'
                latitude = varargin{in+1};
            case 'longitude'
                longitude = varargin{in+1};
            case 'height'
                height = varargin{in+1};
            case 'author'
                author = varargin{in+1};
            case 'nanval_in'
                nanval_in = varargin{in+1};
            case 'nanval_out'
                nanval_out = varargin{in+1};
            case 'out_precision'
                out_precision = varargin{in+1};
        end
        % Increase by 2 as parameters are in pairs!
        in = in + 2;
    end
elseif nargin > 0 && mod(nargin,2) ~= 0
    error('Set even number of input parameters')
end

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
    [timec,datac,head_in] = loadggp('file_in',file_input,'offset',0,...
                            'nanval',nanval_in);   
    % Stack data for processing
    time = vertcat(time,timec);
    data = vertcat(data,datac);
    clc
end

%% Write result
% Use only required time interval
data(time>datenum(end_time),:) = [];
time(time>datenum(end_time),:) = [];
% Get channel names and units
try
    temp = strrep(head_in{end},')','');
    temp_all = strsplit(head_in{end},')');
    temp = strsplit(temp_all{1});
    temp = strsplit(temp{end},'(');
    channels(1) = {temp{1}};
    units(1) = {temp{2}};
    for i = 2:size(data,2)
        temp = strsplit(temp_all{i},'(');
        channels(i) = {temp{1}};
        units(i) = {temp{2}};
    end
catch
    channels = {''};
    units = {''};
    for i = 1:size(data,2)
        channels(i) = {sprintf('data%02d',i)};
        units(i) = {'?'};
    end
end
% Get file suffix to choose output format
[~,file_name,file_suffix] = fileparts(output_file);
if strcmp(file_suffix,'.ggp')
    header = {'Filename',[file_name,file_suffix];...
           'Station',station;...
           'Instrument',instrument,...
           'N. Latitude (deg)', latitude;...
           'E. Longitude (deg)',longitude;...
           'Elevation MSL (m)', height;...
           'Author',author};
    header_add = {'Processing scripts available at: https://github.com/emenems/igetsDataTools'};
    writeggp('time',time,'data',data,'header_offset',21,'header',header,...
          'channels',channels,'units',units,'output_file',output_file,...
          'out_precision',out_precision,'format',file_format,'header_add',header_add,...
          'nanval',nanval_out);
elseif strcmp(file_suffix,'.tsf')
    % Convert out_precision to writetsf format
    temp = strsplit(out_precision,'.');
    out_precision_tsf = str2double(temp{2}(1:end-1));
    header = {};
    % Create header for tsf
    for i = 1:size(data,2)
        header(i,1:4) = {station,instrument,channels{i},units{i}};
    end
    writetsf([datevec(time),data],header,output_file,out_precision_tsf,...
        {'Processing scripts available at: https://github.com/emenems/igetsDataTools'});
elseif strcmp(file_suffix,'.mat')
    data_all.time = time;
    data_all.data = data;
    data_all.channels = channels;
    data_all.units = units;
    save(output_file,'data_all');
end

