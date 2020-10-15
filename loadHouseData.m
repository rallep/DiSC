% Script for loading consumption data of houses and saving it in .mat
% format. The data is sampled every 15 minutes.

clc; clear; close all;

% Setup
dataPath = 'HouseDataCSV/';
savePath = '/Users/rasmus/Documents/uni/benchmark_model_code/svn/matlab/lv_grid/consData/';

fromHouse = 1;
toHouse = 400;
numDays = 2;
% Date format: 'yyyy-mm-dd HH:MM:SS', starts from:  2011-01-01 00:15:00, 
% and is sampled every 15 min.
% To start at different date, set startSample.
% Example:
%   - If start time of wanted data is: 2011-03-01 00:00:00
%   - Then startSample = 5664, which is 59 days.
%   - However, be aware of difference in number of days in February
%   - Check the start date in the Data struct: Data.Time
startSample = 1;    % samplesDay*numDays*month (96 samples pr. day)

numSamples = numDays*24*4+1;
HouseP = zeros(numSamples,toHouse - fromHouse);
j=1;
for i=fromHouse:toHouse
    fileName = ['house' num2str(i) '.csv'];
    fid = fopen([dataPath fileName]);
    data = textscan(fid,'%s %f',numSamples,'delimiter',',','HeaderLines',startSample);
    Power = data{2};
    HouseP(:,j) = Power;
    fclose(fid);
    j = j+1;
end
time = cell2mat(data{1});
xdate = datenum(time,'yy-mm-dd HH:MM:SS');

%% Save data
saveFile = ['house' num2str(fromHouse) 'to' num2str(toHouse) 'days' num2str(numDays) 'startSample' num2str(startSample)];
Data.time = time;
Data.HouseP = HouseP;
curPath = pwd;
cd(savePath)
save(saveFile,'Data')
cd(curPath)