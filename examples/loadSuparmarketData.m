% Script for loading and saving ESO2 - Supermarket Power data in mat format
clc; clear all; close all;

% Setup
month = 'Jan';
year = '2013';
dataPath = [year '/' month '/'];
savePath = '/Users/rasmus/Dropbox/SmartC2Net-Control/benchmark_simulation_framework/matlab/mv_grid/consData/';
fileNameLT = 'frost.csv';
fileNameMT = 'koel.csv';

% Import data form CSV file
dataLT = csvread([dataPath fileNameLT],2,1);
dataMT = csvread([dataPath fileNameMT],2,1);

power = (dataLT(:,3) + dataMT(:,3)).*1000;

% Save file
saveName = ['smPower' month year];
curPath = pwd;
cd(savePath);
save(saveName,'power');
cd(curPath);