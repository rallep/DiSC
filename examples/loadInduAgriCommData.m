% Script for loading consumption data of industry, agriculture, commercial
% and residential. All data is on a hour sampling. There is one week of
% data.
clc; clear; close all;

% Setup
Range = 'B2:B26';   % B is winther and C is summer

%% Read Industry
% file to read
fileName = 'industry.xlsx';
savePath = '/Users/rasmus/Documents/uni/benchmark_model_code/svn/matlab/mv_grid/consData/';

season = Range; % B is winter and C is summer 
pIndu = xlsread(fileName,season);
p = pIndu*1000;
saveName = 'induPower';
curPath = pwd;
cd(savePath);
save(saveName,'p');
cd(curPath);

% season = Range; % B is winter and C is summer 
% p = xlsread(fileName,season);
% p = p*1000;
% saveName = 'induPowerSummer';
% curPath = pwd;
% cd(savePath);
% save(saveName,'p');
% cd(curPath);

%% Read Agriculture
% file to read
fileName = 'agriculture.xlsx';
%savePath = '/Users/rasmus/Dropbox/SmartC2Net-Control/benchmark_simulation_framework/matlab/mv_grid/consData/';

season = Range; % B is winter and C is summer 
p = xlsread(fileName,season);
p = p*1000;
saveName = 'agriPower';
curPath = pwd;
cd(savePath);
save(saveName,'p');
cd(curPath);

% season = 'C2:C169'; % B is winter and C is summer 
% p = xlsread(fileName,season);
% p = p*1000;
% saveName = 'agriPowerSummer';
% curPath = pwd;
% cd(savePath);
% save(saveName,'p');
% cd(curPath);

%% Read Commercial
% file to read
fileName = 'commercial.xlsx';
%savePath = '/Users/rasmus/Dropbox/SmartC2Net-Control/benchmark_simulation_framework/matlab/mv_grid/consData/';

season = Range; % B is winter and C is summer 
p = xlsread(fileName,season);
p = p*1000;
saveName = 'commPower';
curPath = pwd;
cd(savePath);
save(saveName,'p');
cd(curPath);

% season = 'C2:C169'; % B is winter and C is summer 
% p = xlsread(fileName,season);
% p = p*1000;
% saveName = 'commPowerSummer';
% curPath = pwd;
% cd(savePath);
% save(saveName,'p');
% cd(curPath);