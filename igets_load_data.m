%% Load 1 second/minite/hour data + create one file with whole time series
% This script will load one second,minute or hour monthly ggp/igets files 
% to one data matrix that can be then written so a ggp/igets, tsf or mat
% file. 
% Following data structure of INPUT files is required
%     'input_path\YYYY\input_prefix+YYYYMM+input_suffix'
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

%% Main Settings
% Time interval to be loaded
start_time = [2015 03 05 14 00 00];
end_time   = [2017 02 27 15 00 00];
% INPUT File path/name settings
input_path = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3'; % year/month/day will be generated automatically
input_prefix = 'IGETS-IGRAV-RESHOUR-we006-'; % file name prefix
input_suffix = 'r1.ggp';
% Set which channels should be loaded (e.g., gravity and pressure)
input_channels = [1,2]; 
% Input file settings
file_format = 'preterna';
nanval = 99999.999; % Flagged NaN values

% Set OUTPUT file naming. Output file format will be created based on file
% suffix: tsf = tsoft, ggp = ggp/igets, mat = matlab array (plotGrav)
output_file = 'f:\mikolaj\data\wettzell\grav\sg\igrav006\igets\Wettzell\we006\Level3\IGETS-IGRAV-RESHOUR-we006-ALLr1.mat';
% Output header:
head_site = 'Wettzell';
head_instrument = 'iGrav006';
output_nanval = 99999.999; % only for ggp/igets
out_precision = '%10.2f';

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
                            'nanval',99999.999);   
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
           'Station',head_site;...
           'Instrument',head_instrument};
    writeggp('time',time,'data',data,'header_offset',21,'header',header,...
          'channels',channels,'units',units,'output_file',output_file,...
          'out_precision',out_precision,'format',file_format,...
          'nanval',output_nanval);
elseif strcmp(file_suffix,'.tsf')
    % Convert out_precision to writetsf format
    temp = strsplit(out_precision,'.');
    out_precision_tsf = str2double(temp{2}(1:end-1));
    % Create header for tsf
    for i = 1:size(data,2)
        header(i,1:4) = {head_site,head_instrument,channels{i},units{i}};
    end
    writetsf([datevec(time),data],header,output_file,out_precision_tsf);
elseif strcmp(file_suffix,'.mat')
    data_all.time = time;
    data_all.data = data;
    data_all.channels = channels;
    data_all.units = units;
    save(output_file,'data_all');
end

